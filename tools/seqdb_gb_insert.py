#! /usr/bin/env python
"""Script to use Entrez and SeqDB-WS APIs to load GenBank sequences into SeqDB
"""

from Bio import Entrez
import logging.config
import os
import re
import shelve
import sys
import urllib

from config import config_root
import httplib as http_client
import tools_helper
from api.BaseSeqdbApi import UnexpectedContent
from api.ConsensusSequenceApi import ConsensusSequenceApi
from api.FeatureTypeApi import FeatureTypeApi
from api.FeatureApi import FeatureApi
from api.DeterminationApi import DeterminationApi
from api.GeneRegionApi import GeneRegionApi
from api.RawSequenceApi import RawSequenceApi
from api.SeqSourceApi import SeqSourceApi
from api.SpecimenApi import SpecimenApi


def merge(a, b, path=None):
    "merges b into a"
    if path is None:
        path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass  # same leaf value
            else:
                raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a


def format_sequence_name(genbank_id, record):
    """Return the value to be added to the SeqDB Sequence Name field.
    Args:
        genbank_id (str): GenBank Id of the record being processed
        record (obj): GenBank record set retrieved from Entrez
    Kargs:
        None
    Returns:
        str. The formatted name
    Raises:
        None
    """
    logging.info("Formating Sequence Name for GI: {}".format(genbank_id))
#    name = "gi|" + str(genbank_id) + "|"
#    name = name + "gb|" + record[0]["GBSeq_accession-version"] + "|"
#    name = name + record[0]["GBSeq_definition"].replace(", ", "_").replace(
#        "; ", "_").replace(" ", "_").replace(";", "").replace(",", "")
#    return name
    name = record[0]["GBSeq_definition"]
    if len(name) > 255:
        name = name[:252] + "..."
    logging.info("Returning: {}".format(name))
    return name


def format_tracefile_name(**keywds):
    """Return the value to be added to the SeqDB Sequence TracefileName field.
    Args:
        Any number of key = value pairs
    Kargs:
        None
    Returns:
        The urllib urlencoded string representing the string = value pairs
    Raises:
        None
    Currently using the URL encoded query portion of the Entrez URL as this
    field will be updated in SeqDB.
    """
    logging.info("Formatting Tracefile Name")
    result = "?" + urllib.urlencode(keywds)
    logging.info("Returning: {}".format(result))
    return result


def extract_gene_names(record):
    """Return an array of  gene names from the Entrez record.
    Args:
        record (obj): Genbank record retrieved from Entrez
    Kargs:
        None
    Returns:
        str. A list of gene names contained in the record.
    Raises:
        None
    """
    logging.info("Extracting gene names from record: {}".format(record["GBSeq_accession-version"]))
    genes = {}
    for feature in record["GBSeq_feature-table"]:
        if feature["GBFeature_key"] == "gene":
            for qualifier in feature["GBFeature_quals"]:
                # TODO: Should only be looking at GBQualifier_name == "gene"
                if qualifier["GBQualifier_name"] == "gene":
                    genes[qualifier["GBQualifier_value"]] = 1
                    logging.debug("Gene name: {}".format(qualifier["GBQualifier_value"]))
    logging.debug("Found {} Gene Names".format(len(genes.keys())))
    return genes.keys()


def check_feature_type(api_key, url, ftn, create=False, lookup=None):
    """Check to see if SeqDB contains a Feature Type with desired name
    Args:
        ftn (str): name of feature type to lookup
    Kargs:
        create (bool): If True, create feature types which are not found
    Returns:
        int or None. SeqDB feature type id if found/created. None if not
                     found/created or there are multiple gene regions
                     with the same name.
    Raises:
        None
    """
    logging.debug("Checking SeqDB for feature type: {}.  Create == {}".format(ftn, create))
    
    featureTypeApi = FeatureTypeApi(api_key=api_key, base_url=url)

    feature_type_id = None

    if lookup is not None and ftn in lookup:
        feature_type_id = lookup[ftn]

    if feature_type_id is None:
        #feature_types = seqdb_ws.getFeatureTypesWithIds()
        feature_types = featureTypeApi.getFeatureTypesWithIds()
        if ftn in feature_types:
            feature_type_id = feature_types[ftn]
        elif create is True:
            #feature_type_id = seqdb_ws.createFeatureType(ftn, "Genbank Feature Type: %s" % (ftn))
            feature_type_id = featureTypeApi.create(ftn, "GenBank Feature Type: {}".format(ftn))

    if lookup is not None:
        lookup[ftn] = feature_type_id

    logging.debug("Returning Feature Type ID: {}".format(feature_type_id))

    return feature_type_id


# TODO filter by group
def check_region(api_key, url, gene, create=False):
    """Check to see if SeqDB contains a region with a name as per gene.
    Args:
        gene (str): name of region to lookup
    Kargs:
        create (bool): If True, create regions which are not found
    Returns:
        int or None. SeqDB Gene Region Id if found/created. None if not
                     found/created or there are multiple gene regions
                     with the same name.
    Raises:
        None
    """
    logging.debug("Checking SeqDB for region: {}. Create == {}".format(gene, create))
    region_id = None
    geneRegionApi = GeneRegionApi(api_key=api_key, base_url=url)
    geneRegionApi.nameFilter = gene
    #region_ids = seqdb_ws.getRegionIdsByName(gene)
    region_ids = geneRegionApi.getIds()
    if len(region_ids) == 0 and create is True:
        #region_id = seqdb_ws.createRegion(gene, "GenBank Gene: {}".format(gene))
        region_id = geneRegionApi.create(gene, "GenBank Gene: {}".format(gene))
        logging.debug("Created region: {} in SeqDB for {}".format(region_id, gene))
    elif len(region_ids) == 1:
        region_id = region_ids[0]
        logging.debug("Found region: {} in SeqDB for {}".format(region_id, gene))
    else:
        logging.warn("Found multiple regions for '{}' in SeqDB."
                     "Currently unable to assign sequences to regions."
                     "with non-unique names.".format(gene))

    logging.debug("Returning region ID: {}".format(region_id))
    return region_id


def entrez_search(query, retmax=1, retstart=0, database="nucleotide"):
    """Return record set retrieved via an Entrez search.
    Args:
        query (str): Entrez query to execute.
    Kargs:
        retmax (int): Limit number of results returned from Entrez (default=1).
        retstart (int): First result to return from Entrez (default=0).
        database (str): The Entrze database to query (default="nucleotide")
    Returns:
        obj. The entrez record set.
    Raises:
        None
    """
    handle = Entrez.esearch(db=database, retstart=retstart, retmax=retmax, term=query)
    record = Entrez.read(handle)
    handle.close()
    return record


def entrez_fetch(genbank_id, rettype="gb", database="nucleotide", retmode=None, cache=None):
    """Retrieve record retrieved via an Entrez fetch.
    Args:
        genbank_id (str): The genbank id of the record to retrieve.
    Kargs:
        rettype (str): The Entrez rettype (default="gb")
        retmode (str): The Entrez retmode (default="xml")
        database (str): The Entrze database to query (default="nucleotide")
    Returns:
        obj. The Entrez record.
    Raises:
        None
    """
    logging.info("Retrieving record for GI: {} from GenBank".format(genbank_id))

    add_to_cache = False
    record = None

    if cache is not None:
        if genbank_id in cache:
            logging.debug("Using cached record for GI: {}".format(genbank_id))
            record = cache[genbank_id]
        else:
            add_to_cache = True

    if record is None:
        handle_text = Entrez.efetch(db=database, id=genbank_id, rettype=rettype, retmode="text")
        first_line = handle_text.readline().strip()
        possible_error_keywords = ["Error", "Bad", "Cannot", "unavailable", "unable", "is empty"]

        while first_line and any(error_substring in first_line for error_substring in possible_error_keywords):
            print "Entrez Error:{}".format(first_line)
            handle_text = Entrez.efetch(db=database, id=genbank_id, rettype=rettype, retmode="text")
            first_line = handle_text.readline().strip()

        handle_xml = Entrez.efetch(db=database, id=genbank_id, rettype=rettype, retmode="xml")
                
        try:
            record = Entrez.read(handle_xml)
        except Entrez.Parser.ValidationError:
            logging.error("Failed to parse GenBank record for GI:{}".format(genbank_id))
            add_to_cache = False            

        handle_xml.close()
        handle_text.close()

    if add_to_cache:
        logging.debug("Adding record for GI: {} to cache".format(genbank_id))
        cache[genbank_id] = record

    tools_helper.pretty_log_json(record[0], level="debug", message="Retrieved GenBank Record: ")
    return record


def seqdb_ret_entrez_gene_region_id(api_key, url, record, products=None):
    """Retrieve the SeqDB gene region id corresponding to the gene name on this Entrez sequence.
    Args:
        genbank_id (str): The GenBank id of the record being processed
        record (obj): GenBank record retrieved from Entrez
    Kargs:
        None
    Returns:
        int or None. SeqDB Gene Region Id if found.
                     None if not found or if the Entrez record contains
                     multiple gene annotations.
    Raises:
        None
    Note, SeqDB only permits a single gene region to be associated with a
    sequence, hence we warn on and ignore Entrez sequences with more than one
    annotated gene.
    """
    seqdb_gene_region_id = None
    genes = extract_gene_names(record[0])
    if len(genes) > 1:
        #seqdb_gene_region_id = check_region(seqdb_ws, "Multigene Region", create=True)
        seqdb_gene_region_id = check_region(api_key, url, "Multigene Region", create=True)
        logging.debug("Adding multigene sequence, SeqDB region ID: {}".format(seqdb_gene_region_id))

    elif len(genes) == 1:
        #seqdb_gene_region_id = check_region(seqdb_ws, genes[0], create=True)
        seqdb_gene_region_id = check_region(api_key, url, genes[0], create=True)
        logging.debug("Found gene: {}, SeqDB region ID: {}".format(genes[0], seqdb_gene_region_id))

    elif len(genes) == 0 and products is not None:
        if "18S ribosomal RNA" or "internal transcribed spacer 1" or \
                "5.8S ribosomal RNA" or "internal transcribed spacer 2" or \
                "28S ribosomal RNA" in products:
            #seqdb_gene_region_id = check_region(seqdb_ws, "Ribosomal Cistron", create=True)
            seqdb_gene_region_id = check_region(api_key, url, "Ribosomal Cistron", create=True)
            logging.debug("Identified Ribosomal Cistron based on features in "
                          "0 gene region, SeqDB region ID: {}".format(seqdb_gene_region_id))
        else:
            logging.debug("No gene region for 0 gene region, SeqDB region "
                          "ID: {}".format(seqdb_gene_region_id))

    return seqdb_gene_region_id


def seqdb_insert_entrez_sequence(consensusSequenceEntity, genbank_id, record):
    """Insert the GenBank sequence into SeqDB.
    Args:
        genbank_id (str): The GenBank id of the record being inserted
        record (obj): GenBank record retrieved from Entrez
    Kargs:
        None
    Returns:
        int. SeqDB id for the inserted consensus sequence.
    Raises:
        None
    """
    logging.info("Adding sequence for GI: {} to SeqDB".format(genbank_id))
    seq_name = format_sequence_name(genbank_id, record)
    sequence = record[0]["GBSeq_sequence"]
    sequence = sequence.upper()
    tracefile_dir = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

    tracefile_name = format_tracefile_name(db="nucleotide", id=genbank_id, rettype="gb", retmode="xml")

    additional = {
        'fileDir': tracefile_dir,
        'fileName': tracefile_name,
        'genBankGI': genbank_id,
        'genBankAccession': record[0]["GBSeq_primary-accession"],
        'genBankVersion': record[0]["GBSeq_accession-version"],
        'submittedToInsdc': 'true'
    }

    dict_for_logging = {'name': seq_name, 'sequence': sequence}
    dict_for_logging.update(additional)
    tools_helper.pretty_log_json(dict_for_logging, level="debug", message="Creating consensus (non-default values): ")

    #seqdb_id, code, message = seqdb_ws.createConsensusSequence(seq_name, sequence, additional=additional)
    seqdb_id, code, message = consensusSequenceEntity.create(seq_name, sequence, additional=additional)

    logging.info(
        "Created Consensus Sequence (seqdbid: {}) "
        "for GI: {} ({}), Status: {}, Message: {}".format(seqdb_id, genbank_id, record[0]["GBSeq_accession-version"], code, message))

    return seqdb_id


def seqdb_link_to_specimen(api_key, url, seqdb_id, feature):
    """Associate the sequence with a Specimen in SeqDB based on
       supported GBQualifier values.  Currently supported values include:

        * culture collection
        * strain
        * specimen voucher
        * isolate

        Currently looks for entries beginning with "DAOM" or "CBS" in the
        above, or with the prefix "personal:". May need to look for additional
        prefixes and check values for additional qualifier keys.
    """
    logging.info("Linking sequence to available source Specimen")
    specimenApi = SpecimenApi(api_key=api_key, base_url=url)
    for supported in ['culture_collection', 'strain', 'specimen voucer', 'isolate']:

        if supported in feature['qualifiers']:

            for source_entry in feature['qualifiers'][supported]:

                code = None
                identifier = None
                specimenIds = None
                #jsn_resp = None

                # some source entries are a list split with a semi-colon
                sources = source_entry.split(";")
                for source in sources:

                    source = source.strip()

                    if source.startswith("DAOM") or \
                            source.startswith("CCFC") or \
                            source.startswith("CBS") or \
                            source.startswith("ATCC") or \
                            source.startswith("INVAM") or \
                            source.startswith("NISK") or \
                            source.startswith("BR"):

                        logging.info("\tPossible known source {}.".format(source))

                        matched = None

                        if source.startswith("DAOM"):

                            # TODO Also search / instead based on DAOM field in
                            # FungalInfo
                            matched = re.search(r"(?P<collection>\w+)[: ]?(?P<identifier>[\d.]+)", source)

                            if matched.groupdict()['identifier'].startswith("BR"):

                                matched = re.search(r"(?P<collection>BR)(?P<identifier>[\d-]+)", matched.groupdict()['identifier'])

                        elif source.startswith("CCFC"):

                            # TODO Also search / instead based on CCFC field in
                            # FungalInfo
                            matched = re.search(r"CCFC:DAOM (?P<collection>BR) (?P<identifier>[\d.]+)", source)

                        elif source.startswith("BR"):

                            matched = re.search(r"(?P<collection>\w+)[: ]?(?P<identifier>[\d.]+) \(DAOM\)", source)

                        elif source.startswith("INVAM"):

                            matched = re.search(r"(?P<collection>\w+)[ ]?(?P<identifier>[\w\d]+)", source)

                        elif source.startswith("CBS") or \
                                source.startswith("ATCC"):

                            matched = re.search(r"(?P<collection>\w+)[: ]?(?P<identifier>[\d.]+)", source)

                        elif source.startswith("NISK"):

                            matched = re.search(r"(?P<collection>NISK) (?P<identifier>\d+-\d+/\d+)", source)

                        if matched:

                            code = matched.groupdict()['collection']
                            identifier = matched.groupdict()['identifier']

                            logging.debug("\tChecking SeqDB for Specimen with OtherIds containing "
                                          "Code: {}, Identifier: {}".format(code, identifier))

                            try:
                                #jsn_resp = seqdb_ws.getJsonSpecimenIdsByOtherIds(code, identifier)
                                specimenApi.otherIdsFilter = code + identifier
                                specimenIds = specimenApi.getIds()
                                

                            except UnexpectedContent, e:
                                logging.error("Exception querying Specimen using "
                                              "Code: {}, Identifier: {}.  {}".format(code, identifier, e))

                    elif source.startswith("personal:"):

                        logging.info("\tUnevaluated personal identifier {}.".format(source))

                        prefix, code, identifier = source.split(":")

                        if isinstance(identifier, int):

                            logging.debug("\tChecking SeqDB for Specimen "
                                          "Code: {}, Identifier: {}".format(code, identifier))

                            #jsn_resp = seqdb_ws.getJsonSpecimenIdsBySpecimenId(code, identifier)
                            specimenApi.otherIdsFilter = code + identifier
                            specimenIds = specimenApi.getIds()
                            

                        else:

                            logging.warn("\tReview non-numeric identifier: {}".format(source))

                    # the following prefixes don't appear in our otherIds
                    # column could change in the future, but ignoring for now
                    elif source.startswith("RGR"):

                        logging.debug("\tSkiping unused source {}".format(source))

                    else:

                        logging.error("\tUnevaluated source {}".format(source))

                        continue

                    if not specimenIds:

                        logging.info("\tNo specimen found for code: {} and id: {}".format(code, identifier))

                    elif len(specimenIds) > 1:

                        logging.warn("\tUnable to link Sequence to Specimen using "
                                     "Code: {}, Identifier: {}. "
                                     "Found multiple Specimen with OtherIds containing".format(code, identifier))

                    else:  # jsn_resp['count'] == 1:

                        seqdb_update_seqsource_specimen(api_key, url, seqdb_id, specimenIds[0])


def seqdb_link_to_taxonomy(api_key, url, seqdb_id, taxon_id, organism, feature):
    """Associate the sequence with a Taxa in the SeqDB Taxonomy based on the
       GBQualifier "organism" value"""
    
    determinationApi = DeterminationApi(api_key=api_key, base_url=url)

    taxon_id_value = taxon_id.split(":")
    org_split = organism.split(" ")

    taxonomy = {
        "genus":org_split[0],
        "species":org_split[1]
    }

    determinationId = determinationApi.createSequenceDetermination(isAccepted="true", sequenceId=seqdb_id, taxonomy=taxonomy,ncbiTaxonId=taxon_id_value[1], notes="here are notes")
  
    logging.info("Sequence determination: {}".format(organism))


def seqdb_update_seqsource_region(api_key, url, seqdb_id, seqdb_region_id):
    """Associate the sequence with a gene region.
    Args:
        seqdb_id (int): SeqDB id of the sequence to associate with the Region
        seqdb_region_id (int):
                        SeqDB id of the Region to associate with the Sequence
    Kargs:
        None
    Returns:
        int. SeqDB id of updated Sequence (Strikes me as odd, but this
             is the result of the SeqDB WS API)
    Raises:
        None
    """
    logging.info("Linking Sequence (seqdbid: {}) to Region (seqdbid: {})".format(seqdb_id, seqdb_region_id))

    seqsource = {
        "seqSource": {
            "region": {
                "id": seqdb_region_id,
            }
        }
    }
    seqSourceApi = SeqSourceApi(api_key=api_key, base_url=url, sequence_id=seqdb_id)
    #existing = seqdb_ws.getJsonSeqSource(seqdb_id)
    existing = seqSourceApi.retrieveJson(seqSourceApi.request_url)
    region_id = None
    if 'result' in existing:
        # drop headers from response returned above by creating a new dict
        # containing only response 'result' portion and then add additional
        # properties to it
        existing = {'seqSource': existing['result']}

        tools_helper.pretty_log_json(seqsource, level="debug", message="Merging")
        tools_helper.pretty_log_json(existing, level="debug", message="Into")

        merge(existing, seqsource)

        region_id, code, message = seqSourceApi.update(seqSourceApi.request_url, existing)

        logging.debug("Updated SeqSource for sequence region linking "
                      "seqdbid: {}, Status: {}, Message: {}".format(seqdb_id, code, message))

    else:

        tools_helper.pretty_log_json(existing, level="error", message="Failed to retrieve seqSource for Sequence (seqdbid: {}):".format(seqdb_id))

    return region_id


def seqdb_update_seqsource_specimen(api_key, url, seqdb_id, seqdb_specimen_id):
    """Associate the sequence with a specimen.
    Args:
        seqdb_id (int): SeqDB id of the sequence to associate with the Region

        seqdb_specimen_id (int):
                        SeqDB id of the Specimen to associate with the Sequence
    Kargs:
        None
    Returns:
        int. SeqDB id of updated Sequence (Strikes me as odd, but this
             is the result of the SeqDB WS API)
    Raises:
        None
    """
    logging.info("Linking Sequence (seqdbid: {}) to Specimen (seqdbid: {})".format(seqdb_id, seqdb_specimen_id))
    
    seqSourceApi = SeqSourceApi(api_key=api_key, base_url=url, sequence_id=seqdb_id)


    seqsource = {
        "seqSource": {
            "specimen": {
                "id": seqdb_specimen_id,
            },
        }
    }

    #existing = seqdb_ws.getJsonSeqSource(seqdb_id)
    existing = seqSourceApi.retrieveJson(seqSourceApi.request_url)
    region_id = None
    if 'result' in existing:
        # drop headers from response returned above by
        # creating a new dict containing only response 'result' portion
        existing = {'seqSource': existing['result']}
        merge(existing, seqsource)
        #region_id, code, message = seqdb_ws.updateSeqSource(seqdb_id, existing)
        region_id, code, message = seqSourceApi.update(seqSourceApi.request_url, existing)
        logging.debug("Updated SeqSource for sequence specimen linking "
                      "seqdbid: {}, Status: {}, Message: {}".format(seqdb_id, code, message))

    else:
        tools_helper.pretty_log_json(existing, level="error", message="Failed to retrieve seqSource for Sequence (seqdbid: {}):".format(seqdb_id))

    return region_id


def parse_locations(gb_feature):
    """Parse the GBFeature_intervals block and create an array of locations in
       the local dictionary representation.
    Args:
        gb_feature: Reference to GBFeature block from Entrez XML
    Kargs:
        None
    Returns:
        Array of locations extracted from the Entrez GBFeature block.
    Raises:
        None
    """
    locations = []
    #print "******GBFeature_interval: {}\n".format(gb_feature['GBFeature_intervals'])
    for interval in gb_feature['GBFeature_intervals']:
        
        #According to Entrez documentation (example: https://github.com/biopython/biopython/blob/master/Tests/Entrez/nucleotide1.xml)
        #'GBInterval_from' and 'GBInterval_to' may not exist for a record.
        
        #if 'GBInterval_from' in interval and 'GBInterval_to' in interval:
            # TODO determined frame and strand, don't just default to 1
            # There is another spot where I adjust the frame
        
        #while 'GBInterval_from' not in interval and 'GBInterval_to' not in interval:
        #    print "=======Need to wait for GBInterval_from / to"
        #    time.sleep(3)
                
        locations.append(
            {
                "start": interval['GBInterval_from'],
                "end": interval['GBInterval_to'], "frame": 1, "strand": 1
            })
    return locations


def parse_qualifiers(gb_feature):
    """Parse the GBFeature_quals block of the GBFeature entry and create a
       dictionary of qualifiers using the GBQualifier_name as a key.
    Args:
        gb_feature: Reference to GBFeature block from Entrez XML
    Kargs:
        None
    Raises:
        None
    """
    qualifiers = {}
    #print "*******GBQualifier: {}\n".format(gb_feature['GBFeature_quals'])
    for qual in gb_feature['GBFeature_quals']:
        if qual['GBQualifier_name'] not in qualifiers:
            qualifiers[qual['GBQualifier_name']] = []
        qualifiers[qual['GBQualifier_name']].append(qual['GBQualifier_value'])
    return qualifiers


def parse_feature(gb_feature, api_key, url, lookup=None):
    """Parse the GBFeature to create a dict from the record.
    Args:
        gb_feature: Reference to GBFeature block from Entrez XML
    Kargs:
        lookup (obj): Reference to a dict to locally cache feature type entries
    Returns:
        Array of locations extracted from the Entrez GBFeature block.
    Raises:
        None
    """
    logging.debug("Parsing feature: \"{}\"".format(gb_feature['GBFeature_key']))
    gb_feature_record = {}
    gb_feature_record['location_description'] = gb_feature['GBFeature_location']
    gb_feature_record['feature_key'] = gb_feature['GBFeature_key']
    gb_feature_record['feature_type_id'] = check_feature_type(api_key, url, gb_feature['GBFeature_key'], create=True, lookup=lookup)
    gb_feature_record['locations'] = parse_locations(gb_feature)
    gb_feature_record['qualifiers'] = parse_qualifiers(gb_feature)
    
    return gb_feature_record


def adjust_feature_for_codon_start(feature):
    """Adjust feature start if feature contains a "codon_start" annotation.
    Args:
        feature (dict): Feature entry from internal representation
    Kargs:
        None
    Returns:
        The (possibly updated) feature.
    Raises:
        None
    """
    # adjust frame for feature based on codon_start qualifier
    if 'codon_start' in feature['qualifiers']:
        # doesn't make sense to me that there would be more than one
        # codon_start value
        if len(feature['qualifiers']['codon_start']) > 1:
            logging.warn(
                "Multiple codon_start feature qualifiers found.  Using first.")
        feature['locations'][0]['frame'] = feature[
            'qualifiers']['codon_start'][0]
    return feature


def process_features(seqdb_id, record, api_key, url, lookup=None):
    """Process features contained in Entrez GBSeq_feature-table and add them to
       the sequence in seqdb.
    Args:
        seqdb_id (int): SeqDB id of the sequence
        record (obj):   reference to the GBSeq block of the Entrez record

    Kargs:
        lookup (obj):   Reference to a dict to locally cache feature type entries
    Returns:
        dict containing all unique "product" GBFeature_quals  observed on this record.
    Raises:
        None
    """
    logging.info("Adding features from Entry: {} to Sequence (seqdbid: {})".format(record['GBSeq_accession-version'], seqdb_id))

    features = []
    for gb_feature in record['GBSeq_feature-table']:
        #while 'GBFeature_intervals' not in gb_feature:
        #    print "=======Need to wait for GBFeature_intervals"
        #    time.sleep(3)
        #while 'GBFeature_quals' not in gb_feature:
        #    print "======Need to wait for GBFeature_quals"
        #    time.sleep(5)
        features.append(parse_feature(gb_feature, api_key, url, lookup=lookup))

    products = {}
    gene_id = None
    mrna_id = None
    cds_id = None
    for feature in features:
        # create hash of unique product names for future use (returned from
        # function)
        if 'product' in feature['qualifiers']:
            for product in feature['qualifiers']['product']:
                products[product] = 1

        # default / initial name and description
        name = "unnamed " + feature['feature_key']
        description = feature['location_description']

        # update name, based on one of the following keys in order of
        # preference don't expect gene will ever have a 'product', but it
        # doesn't hurt anything
        name_keys = ['product', 'protein_id', 'gene', 'locus_tag']
        for key in name_keys:
            if key in feature['qualifiers']:
                # if the feature appeared multiple times, use the first entry
                # value for the name
                name = feature['qualifiers'][key][0]
                break

        # supplement description (prepend key value pairs to base description)

        taxon_id = None;
        organism_name = None;

        description_keys = ['db_xref', 'protein_id', 'locus_tag', 'note', 'rpt_type', 'gene', 'product', 'organism']
        for key in description_keys:
            if key in feature['qualifiers']:
                for value in feature['qualifiers'][key]:
                    description = key + ": " + value + "; " + description
                    #logging.warn("description: %s" % description)
                    if key =='db_xref':
                        taxon_id = value
                    if key =='organism':
                        organism_name = value

        # currently unsupported features that I've encountered
        #    1  STS
        #  276  assembly_gap
        #   46  exon
        #    1  gap
        #   39  intron
        #    1  mRNA
        #    4  misc_RNA
        #   98  misc_difference
        #    2  misc_feature
        # 2650  mobile_element
        #    1  source
        #   20  stem_loop
        #   39  tmRNA

        # Assumes I will encounter another gene before a feature that is not a
        # child of these gene.
        # TODO check range of parent and null gene/cds/mrna ids once we are
        # outside the range

        if feature['feature_key'] == 'gene':
            #gene_id = seqdb_ws.insertFeature(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description)
            featureApi = FeatureApi(api_key=api_key, base_url=url)
            gene_id = featureApi.create(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description)
            
        elif feature['feature_key'] == 'mRNA':
            feature = adjust_feature_for_codon_start(feature)
            #mrna_id = seqdb_ws.insertFeature(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description, parentId=gene_id)
            featureApi = FeatureApi(api_key=api_key, base_url=url)
            mrna_id = featureApi.create(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description, parentId=gene_id)

        elif feature['feature_key'] == 'CDS':
            parent_id = mrna_id
            if parent_id is None:
                parent_id = gene_id

            feature = adjust_feature_for_codon_start(feature)
            #cds_id = seqdb_ws.insertFeature(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description, parentId=parent_id)
            featureApi = FeatureApi(api_key=api_key, base_url=url)
            cds_id = featureApi.create(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description, parentId=parent_id)

        # TODO do these necessarily have a parent gene?
        elif feature['feature_key'] in ['tRNA', 'rRNA', 'misc_RNA']:
            #seqdb_ws.insertFeature(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description, parentId=gene_id)
            featureApi = FeatureApi(api_key=api_key, base_url=url)
            featureApi.create(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description, parentId=gene_id)


        elif feature['feature_key'] in ['repeat_region', 'misc_feature', 'misc_difference']:
            #seqdb_ws.insertFeature(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description)
            featureApi = FeatureApi(api_key=api_key, base_url=url)
            featureApi.create(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description)


        elif feature['feature_key'] == 'source':
            # seqdb_link_to_specimen(seqdb_ws, seqdb_id, feature)
            seqdb_link_to_taxonomy(api_key, url, seqdb_id, taxon_id, organism_name, feature)

        else:
            logging.warn("Unsupported feature type: {}".format(feature['feature_key']))

    return products


def process_entrez_entry(consensusSequenceEntity, api_key, url, genbank_id, cache=None, lookup=None, delete=False, update=False):
    """Process an Entrez entry.
    Args:
        genbank_id (str): The GenBank id of the record being processed
    Kargs:
        cache (obj):    a "stash" in which to cache Entrez results returned from GenBank.
        lookup (dict):  a dict to hold SeqDB Features and save queries.
        delete (bool):  Default False. Delete existing SeqDB records and recreate them.
    Returns:
        None
    Raises:
        None
    Overview:
        Check if the entry already exists in SeqDB.  If so, skip this entry.
        Otherwise:
            * Insert the sequence into SeqDB.
            * If this GenBank entry contains a single gene region, and it is
              present in SeqDB (or we have created the corresponding gene
              region), associate the Sequence with the Gene Region.
            * If this GenBank entry contains feature annotations, add them
              to SeqDB
    """
    logging.info("Processing GI: {}".format(genbank_id))

    # Ensure the record doesn't already exist in SeqDB
    # If it does, continue with the next record
    logging.debug("Checking for GI: {} in SeqDB".format(genbank_id))

    #result = seqdb_ws.getJsonConsensusSequenceIdsByGI(genbank_id)
    #seq_ids = seqdb_ws.getConsensusSequenceIds(genBankGI=genbank_id)
    consensusSequenceEntity.genBankGIFilter = genbank_id
    seq_ids = consensusSequenceEntity.getIds()
    
    if seq_ids and delete:
        logging.info("Deleting existing Sequence (seqdbid: {})".format(seq_ids[0]))
        #seqdb_ws.deleteConsensusSequence(seq_ids[0])
        consensusSequenceEntity.delete(seq_ids[0])
        
    record = None
    if not seq_ids or update:
        # retrieve GenBank record
        record = entrez_fetch(genbank_id, cache=cache)

    seqdb_id = None
    if seq_ids and update:
        seqdb_id = seq_ids[0]

    elif seq_ids:
        logging.info("Sequence for GI: {} already exists in SeqDB. Skipping.".format(genbank_id))

    if not seq_ids:

        if "GBSeq_sequence" in record[0]:
            seqdb_id = seqdb_insert_entrez_sequence(consensusSequenceEntity, genbank_id, record)

        else:
            logging.info("Skipping GI: {}, which does not contain a sequence.".format(genbank_id))

    if record is not None and seqdb_id is not None:
        features = process_features(seqdb_id, record[0], api_key, url, lookup=lookup)

        seqdb_gene_region_id = seqdb_ret_entrez_gene_region_id(api_key, url, record, features)

        if seqdb_gene_region_id is not None:
            
            seqdb_update_seqsource_region(api_key, url, seqdb_id, seqdb_gene_region_id)

        if logging.getLogger().isEnabledFor(logging.DEBUG):
            # retrieve inserted record and display to users for validation
            # purposes
            #result = seqdb_ws.getJsonSequence(seqdb_id)
            rawSequenceEntity = RawSequenceApi(api_key=api_key, base_url=url)
            result = rawSequenceEntity.retrieveJson(rawSequenceEntity.request_url)
            tools_helper.pretty_log_json(result, level="debug", message="Final Consensus Sequence:")


def main():
    """Load sequences matching Entrez query into SeqDB.
    Args:
        None
    Kargs:
        None
    Returns:
        None
    Raises:
        None
    """
    print("Loading configuration file: {}".format(config_root.path()) + '/config.yaml')
    print("Loading tools configuration file: {}".format(os.path.dirname(__file__)) + '/seqdb_gb_insert_config.yaml')
    
    main_config = tools_helper.load_config(config_root.path() + '/config.yaml')
    
    if not main_config:
        logging.error(tools_helper.log_msg_noConfig)
        sys.exit(tools_helper.log_msg_sysExit) 

    tool_config = tools_helper.load_config(os.path.dirname(__file__) + '/seqdb_gb_insert_config.yaml')

    if not tool_config:
        logging.error(tools_helper.log_msg_noConfig)
        sys.exit(tools_helper.log_msg_sysExit)

    url = main_config['seqdb']['url'] 
    api_key = tool_config['seqdb']['api_key']

    logging.config.dictConfig(main_config['logging'])
    
    http_client.HTTPConnection.debuglevel = main_config['http_connect']['debug_level']

    # caching the entrez records shaved 2 minutes off the time to load
    # ~740 sequences from query: "(*DAOM*[source] and levesque and not
    # 'unplaced genomics scaffold')"
    # real    11m40.754s
    # user    1m31.726s
    # sys     0m14.760s
    # - vs -
    # real    9m21.112s
    # user    1m27.726s
    # sys     0m13.619s
    entrez_cache = shelve.open(tool_config['entrez']['cache'])

    # however, caching the lookup shaved an additional ~7 minutes off the total
    # time to load above query
    # real    2m35.773s
    # user    0m16.539s
    # sys     0m2.486s
    # TODO May not be necessary any longer; instead use API lookup by feature
    # type name
    feature_type_lookup = {}

    logging.info("Script executed with the following command and arguments: {}".format(sys.argv))

    consensusSequenceEntity = ConsensusSequenceApi(api_key=api_key, base_url=url)

    Entrez.email = tool_config['entrez']['email']
    query = tool_config['entrez']['query']

    logging.info("Querying GenBank: \'{}\'".format(tool_config['entrez']['query']))

    # preliminary query to find out how many records there are
    record = entrez_search(query)

    # setup loop counters; retrieving records 50 at a time
    count = int(record["Count"])
    start = 0
    retrieve = 50
    logging.info("Query returned {} records. Retrieving them in batches of {}".format(count, retrieve))

    # repeat until we have all records
    while start < count:

        print 'Count:' + str(count)
        print 'Start:' + str(start)
        
        # retrieve block of records
        logging.debug("Retrieving {}..{}".format(start, start + retrieve))
        record = entrez_search(query, retstart=start, retmax=retrieve)

        # process each returned id in the batch
        for genbank_id in record["IdList"]:
            process_entrez_entry(
                consensusSequenceEntity,
                api_key,
                url,
                genbank_id,
                cache=entrez_cache,
                lookup=feature_type_lookup,
                delete=tool_config['gb_insert']['delete'],
                update=tool_config['gb_insert']['update'])
            print ("\n >Seqid: {}".format(genbank_id))

        start += retrieve
    
    
    print "***Done***"

if __name__ == "__main__":
    main()
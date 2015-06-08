#! /usr/bin/env python
"""Script to use Entrez and SeqDB-WS APIs to load GenBank sequences into SeqDB
"""

import sys
import json
import shelve
import urllib
import logging.config
import httplib as http_client
import tools_helper

from api.seqdbWebService import seqdbWebService
from Bio import Entrez



# Note: Not logging method call as method is a logging helper




def format_sequence_name(genbank_id, record):
    """Return the value to be added to the SeqDB Sequence Name field.

    Args:
        genbank_id (str): Genbank Id of the record being processed
        record (obj): Genbank record set retrieved from Entrez

    Kargs:
        None

    Returns:
        str. The formatted name

    Raises:
        None
    """
    logging.info("Formating sequence name for gi:%s" % (genbank_id))
#    name = "gi|" + str(genbank_id) + "|"
#    name = name + "gb|" + record[0]["GBSeq_accession-version"] + "|"
#    name = name + record[0]["GBSeq_definition"].replace(", ", "_").replace(
#        "; ", "_").replace(" ", "_").replace(";", "").replace(",", "")
#    return name
    name = record[0]["GBSeq_definition"]
    if len(name) > 255:
        name = name[:252] + "..."
    logging.info("Returning: %s" % (name))
    return name


def format_tracefile_name(**keywds):
    """Return the value to be added to the SeqDB Sequence TracefileName field.

    Args:
        Any number of key=value pairs

    Kargs:
        None

    Returns:
        The urllib urlencoded string representing the string=value pairs

    Raises:
        None

    Currently using the URL encoded query portion of the entrez URL as this
    field will be updated in SeqDB.
    """
    logging.info("Formatting tracefile name")
    result = "?" + urllib.urlencode(keywds)
    logging.info("Returning: %s" % (result))
    return result


def extract_gene_names(record):
    """Return an array of  gene names from the Entrez record.

    Args:
        record (obj): Genbank record retrieved from Entrez

    Kargs:
        None

    Returns:
        [str]. A list of gene names contained in the record.

    Raises:
        None
    """
    logging.info("Extracting gene names from record: %s" %
                 (record["GBSeq_accession-version"]))
    genes = {}
    for feature in record["GBSeq_feature-table"]:
        if feature["GBFeature_key"] == "gene":
            for qualifier in feature["GBFeature_quals"]:
                # TODO: Should only be looking at GBQualifier_name == "gene"
                if qualifier["GBQualifier_name"] == "gene":
                    genes[qualifier["GBQualifier_value"]] = 1
                    logging.debug("Gene name: %s" %
                                  (qualifier["GBQualifier_value"]))
    logging.info("Found %i gene names" % (len(genes.keys())))
    return genes.keys()


def check_feature_type(seqdb_ws, ftn, create=False, lookup=None):
    """Check to see if SeqDB contains a Feature Type with desired name

    Args:
        seqdb_ws (obj): reference to instance of api.seqdbWebService
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
    logging.info("Checking SeqDB for feature type: %s.  Create == %r" %
                 (ftn, create))

    feature_type_id = None

    if lookup is not None and ftn in lookup:
        feature_type_id = lookup[ftn]

    if feature_type_id is None:
        feature_types = seqdb_ws.getFeatureTypesWithIds()
        if ftn in feature_types:
            feature_type_id = feature_types[ftn]
        elif create is True:
            feature_type_id = seqdb_ws.createFeatureType(ftn, "Genbank Feature Type: %s" % (ftn))

    if lookup is not None:
        lookup[ftn] = feature_type_id

    logging.info("Returning feature type id: %r" % (feature_type_id))

    return feature_type_id

# TODO filter by group
def check_region(seqdb_ws, gene, create=False):
    """Check to see if SeqDB contains a region with a name as per gene.

    Args:
        seqdb_ws (obj): reference to instance of api.seqdbWebService
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
    logging.info("Checking SeqDB for region: %s.  Create == %r" %
                 (gene, create))
    region_id = None
    region_ids = seqdb_ws.getRegionIdsByName(gene)
    if len(region_ids) == 0 and create is True:
        region_id = seqdb_ws.createRegion(gene, "GenBank Gene: %s" % (gene))
        logging.debug("Created region: %i in seqdb for %s" % (region_id, gene))
    elif len(region_ids) == 1:
        region_id = region_ids[0]
        logging.debug("Found region: %i in seqdb for %s" % (region_id, gene))
    else:
        logging.warn("Found multiple regions for '%s' in SeqDB. "
                     "Currently unable to assign sequences to regions "
                     "with non-unique names." % (gene))

    logging.info("Returning region id: %r" % (region_id))
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
    handle = Entrez.esearch(
        db=database, retstart=retstart, retmax=retmax, term=query)
    record = Entrez.read(handle)
    handle.close()
    return record

def entrez_fetch(
        genbank_id, rettype="gb", retmode="xml", database="nucleotide", cache=None):
    """Retrieve record retrieved via an Entrez fetch.

    Args:
        genbank_id (str): The genbank id of the record to retrieve.

    Kargs:
        rettype (str): The Entrez rettype (default="gb")
        retmode (str): The Entrez retmode (default="xml")
        database (str): The Entrze database to query (default="nucleotide")

    Returns:
        obj. The Entreze record.

    Raises:
        None
    """
    logging.info("Retrieving record for gi:%s from GenBank" % (genbank_id))

    add_to_cache = False
    record = None

    if cache is not None:
        if cache.has_key(genbank_id):
            logging.debug("Using cached record for gi:%s" % (genbank_id))
            record = cache[genbank_id]
        else:
            add_to_cache = True
        
    if record is None:
        handle = Entrez.efetch(
            db=database, id=genbank_id, rettype=rettype, retmode=retmode)

        try:
            record = Entrez.read(handle)
        except Bio.Entrez.Parser.ValidationError:
            logging.error("Failed to parse genbank record for gi:%s" %(genbank_id))
            add_to_cache = False

        handle.close()

    if add_to_cache:
        logging.debug("Adding record for gi:%s to cache" % (genbank_id))
        cache[genbank_id] = record

    tools_helper.pretty_log_json(
        record[0], level="debug", message="Retrieved Genbank Record: ")
    return record


def seqdb_ret_entrez_gene_region_id(seqdb_ws, genbank_id, record, products=None):
    """Retrieve the seqdb gene region id corresponding to the gene name on this
    entrez sequence.

    Args:
        seqdb_ws (obj): reference to instance of api.seqdbWebService
        genbank_id (str): The genbank id of the record being processed
        record (obj): Genbank record retrieved from Entrez

    Kargs:
        None

    Returns:
        int or None. SeqDB Gene Region Id if found.
                     None if not found or if the Entrez record contains
                     multiple gene annotations.

    Raises:
        None

    Note, SeqDB only permits a single gene region to be associated with a
    sequence, hence we warn on and ignore entrez sequences with more than one
    annotated gene.
    """
    seqdb_gene_region_id = None
    genes = extract_gene_names(record[0])
    if len(genes) > 1:
        seqdb_gene_region_id = check_region(seqdb_ws, "Multigene Region", create=True)
        logging.debug("Adding multigene sequence, SeqDB region id: %i" % (seqdb_gene_region_id))
    elif len(genes) == 1:
        seqdb_gene_region_id = check_region(seqdb_ws, genes[0], create=True)
        logging.debug("Found gene: %s, SeqDB region id: %s" %
                      (genes[0], seqdb_gene_region_id))
    elif len(genes) == 0 and products is not None:
        if "18S ribosomal RNA" or "internal transcribed spacer 1" or "5.8S ribosomal RNA" or "internal transcribed spacer 2" or "28S ribosomal RNA" in products:
            seqdb_gene_region_id = check_region(seqdb_ws, "Ribosomal Cistron", create=True)
            logging.debug("Identified Ribosomal Cistron based on features in 0 gene region, SeqDB region id: %i" %
                        (seqdb_gene_region_id))
        else:
            logging.debug("No gene region for 0 gene region, SeqDB region id: %i" %
                        (seqdb_gene_region_id))
            

    return seqdb_gene_region_id


def seqdb_insert_entrez_sequence(seqdb_ws, genbank_id, record):
    """Insert the genbank sequence into SeqDB.

    Args:
        seqdb_ws (obj): reference to instance of api.seqdbWebService
        genbank_id (str): The genbank id of the record being inserted
        record (obj): Genbank record retrieved from Entrez

    Kargs:
        None

    Returns:
        int. SeqDB id for the inserted consensus sequence.

    Raises:
        None
    """
    logging.info("Adding sequence for gi:%s to SeqDB" % (genbank_id))
    seq_name = format_sequence_name(genbank_id, record)
    sequence = record[0]["GBSeq_sequence"]
    sequence = sequence.upper()
    tracefile_dir = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    tracefile_name = format_tracefile_name(
        db="nucleotide", id=genbank_id, rettype="gb", retmode="xml")
    additional = {
        'chromatDir': tracefile_dir,
        'chromatName': tracefile_name,
        'genBankGI': genbank_id,
        'genBankAccession': record[0]["GBSeq_primary-accession"],
        'genBankVersion': record[0]["GBSeq_accession-version"],
        'submittedToInsdc': 'true'
    }

    dict_for_logging = {'name': seq_name, 'sequence': sequence}
    dict_for_logging.update(additional)
    tools_helper.pretty_log_json(dict_for_logging, level="debug",
                    message="Creating consensus (non-default values): ")

    seqdb_id, code, message = seqdb_ws.createConsensusSequence(
        seq_name, sequence, additional=additional)

    logging.info(
        "Created consensus seqdbid:%i "
        "for gi:%s (%s), Status: %i, Message: %s" % (
            seqdb_id, genbank_id,
            record[0]["GBSeq_accession-version"],
            code, message))

    return seqdb_id

def seqdb_link_to_specimen(seqdb_ws, seqdb_id, feature):
    """Associate the sequence with a Specimen in SeqDB based on
       supported GBQualifier values.  Currently supported values include:

        * culture collection
        * strain
        * specimen voucher
        * isolate

        Currently looks for entries beginning with "DAOM" or "CBS" in the
        above, or with the prefix "personal:".
    """
    for qual in feature['GBFeature_quals']:

        if qualifier['GBQualifier_name'] == "culture collection":

            if qualifier['GBQualifer_value'].startswith("DAOM"):
                logging.warn(
                    "Found apparent DAOM identifer (%s), but Specimen "
                     "linking not yet implemented" % 
                     (qualifier['GBQualifer_value'])) 

            if qualifier['GBQualifer_value'].startswith("DAOM"):
                logging.warn(
                    "Found apparent DAOM identifer (%s), but Specimen "
                     "linking not yet implemented" % 
                     (qualifier['GBQualifer_value'])) 

            if qualifier['GBQualifer_value'].startswith("personal:"):
                logging.warn(
                    "Found candidate personal identifer (%s), but Specimen "
                    "linking not yet implemented" % 
                    (qualifier['GBQualifer_value']))

def seqdb_link_to_taxonomy(seqdb_ws, seqdb_id, feature):
    """Associate the sequence with a Taxa in the SeqDB Taxonomy based on the
       GBQualifier "organism" value"""

    for qual in gb_feature['GBFeature_quals']:
        if qualifier['GBQualifier_name'] == "organism":
            logging.warn("Taxonomy linking not yet implemented")
            break
    

def seqdb_update_seqsource(seqdb_ws, seqdb_id, seqdb_region_id):
    """Associate the sequence with a gene region.

    Args:
        seqdb_ws (obj): reference to instance of api.seqdbWebService
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
    logging.info(
        "Updating SeqSource for consensus sequence seqdbid:%i" % (seqdb_id))
    seqsource = {
        "seqSource": {
            "region": {
                "id": seqdb_region_id,
            }
        }
    }
    region_id, code, message = seqdb_ws.updateSeqSource(seqdb_id, seqsource)
    logging.info(
        "Updated SeqSource for consensus sequence"
        "seqdbid:%i, Status: %i, Message: %s" % (
            seqdb_id, code, message))

    return region_id

def parse_locations(gb_feature):
    locations = []
    for interval in gb_feature['GBFeature_intervals']:
        # TODO determined frame and strand, don't just default to 1
        locations.append({"start":interval['GBInterval_from'], "end":interval['GBInterval_to'], "frame":1, "strand":1})
    return locations

def parse_qualifiers(gb_feature):
    qualifiers = {}
    for qual in gb_feature['GBFeature_quals']:
        if qual['GBQualifier_name'] not in qualifiers:
            qualifiers[qual['GBQualifier_name']] = []
        qualifiers[qual['GBQualifier_name']].append(qual['GBQualifier_value'])
    return qualifiers

def parse_feature(gb_feature, seqdb_ws, lookup=None):
    logging.info("\tParsing feature: \"%s\"" % (gb_feature['GBFeature_key']))
    gb_feature_record = {}
    gb_feature_record['location_description'] = gb_feature['GBFeature_location']
    gb_feature_record['feature_key'] = gb_feature['GBFeature_key']
    gb_feature_record['feature_type_id'] = check_feature_type(seqdb_ws, gb_feature['GBFeature_key'], create=True, lookup=lookup)
    gb_feature_record['locations'] = parse_locations(gb_feature)
    gb_feature_record['qualifiers'] = parse_qualifiers(gb_feature)
    return gb_feature_record

def adjust_feature_for_codon_start (feature):
    # adjust frame for feature based on codon_start qualifer
    if 'codon_start' in feature['qualifiers']:
        # doesn't make sense to me that there would be more than one codon_start value
        if len(feature['qualifiers']['codon_start']) > 1:
            logging.warn("Multiple codon_start feature qualifiers found.  Using first.")
        feature['locations'][0]['frame'] = feature['qualifiers']['codon_start'][0]
    return feature
    
def process_features(seqdb_ws, seqdb_id, record, lookup=None):
    logging.info("Adding features from Entry: %s to Sequence seqdbid:%i" % (record['GBSeq_accession-version'], seqdb_id))

    features = []
    for gb_feature in record['GBSeq_feature-table']:
        features.append(parse_feature(gb_feature, seqdb_ws, lookup=lookup))

    products = {}
    gene_id = None
    mrna_id = None
    cds_id = None
    for feature in features:
        # create hash of unique product names for future use (returned from function)
        if 'product' in feature['qualifiers']:
            for product in feature['qualifiers']['product']:
                products[product] = 1;

        # default / initial name and description
        name = "unnamed " + feature['feature_key']
        description = feature['location_description']

        # update name, based on one of the following keys in order of preference
        # don't expect gene will ever have a 'product', but it doesn't hurt anything
        name_keys = ['product', 'protein_id', 'gene', 'locus_tag']
        for key in name_keys:
            if key in feature['qualifiers']:
                # if the feature appeared multiple times, use the first entry value for the name
                if (isinstance feature['qualifiers'][key], list):
                    name = feature['qualifiers'][key][0]
                else:
                    name = feature['qualifiers'][key]
                break
 
        # supplement description (prepend key value pairs to base description)
        description_keys = ['db_xref', 'protein_id', 'locus_tag', 'note', 'rpt_type', 'gene', 'product']
        for key in description_keys:
            if key in feature['qualifiers']:
                if (isinstance feature['qualifiers'][key], list):
                    for value in feature['qualifiers'][key]:
                        description = key + ": " + value + "; " + description
                else:
                    description = key + ": " + feature['qualifiers'][key] + "; " + description

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

        # Assumes I will encounter another gene before a feature that is not a child of these gene.
        # TODO check range of parent and null gene/cds/mrna ids once we are outside the range
        if feature['feature_key'] == 'gene':
            gene_id = seqdb_ws.insertFeature(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description)

        elif feature['feature_key'] == 'mRNA':
            feature = adjust_feature_for_codon_start(feature)
            mrna_id = seqdb_ws.insertFeature(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description, parentId=gene_id)

        elif feature['feature_key'] == 'CDS':
            parent_id = mrna_id
            if (parent_id is None):
                parent_id = gene_id
            
            feature = adjust_feature_for_codon_start(feature)
            cds_id = seqdb_ws.insertFeature(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description, parentId=parent_id)

        # TODO do these necessarily have a parent gene?
        elif feature['feature_key'] in ['tRNA', 'rRNA', 'misc_RNA']:
            seqdb_ws.insertFeature(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description, parentId=gene_id)

        elif feature['feature_key'] in ['repeat_region', 'misc_feature', 'misc_difference']:
            seqdb_ws.insertFeature(name, feature['feature_type_id'], feature['locations'], seqdb_id, description=description)

        elif feature['feature_key'] == 'source':
            seqdb_link_to_specimen(seqdb_ws, seqdb_id, feature)
            seqdb_link_to_taxon(seqdb_ws, seqdb_id, feature)

        else:
            logging.warn("Unsupported feature type: %s" %(feature['feature_key']))

    return products
        

def process_entrez_entry(seqdb_ws, genbank_id, cache=None, lookup=None):
    """Process an entrez entry.

    Args:
        seqdb_ws (obj): reference to instance of api.seqdbWebService
        genbank_id (str): The genbank id of the record being processed

    Kargs:
        None

    Returns:
        None

    Raises:
        None

    Overview:
        Check if the entry already exists in SeqDB.  If so, skip this entry.
        Otherwise:
            * Insert the sequence into SeqDB.
            * If this Genbank entry contains a single gene region, and it is
              present in SeqDB (or we have created the correponding gene
              region), associate the Sequence with the Gene Region.
            * If this Genbank entry contains feature annotations, add them
              to SeqDB
    """
    logging.info("Processing gi:%s" % (genbank_id))

    # Ensure the record doesn't already exist in SeqDB
    # If it does, continue with the next record
    logging.debug("Checking for gi:%s in SeqDB" % (genbank_id))

    result = seqdb_ws.getJsonConsensusSequenceIdsByGI(genbank_id)

    if result['count'] == 0:
        # retrieve genbank record
        record = entrez_fetch(genbank_id, cache=cache)

        if "GBSeq_sequence" in record[0]:
            seqdb_id = seqdb_insert_entrez_sequence(seqdb_ws, genbank_id, record)

            features = process_features(seqdb_ws, seqdb_id, record[0], lookup=lookup)

            seqdb_gene_region_id = seqdb_ret_entrez_gene_region_id(seqdb_ws, genbank_id, record, features)
            if seqdb_gene_region_id is not None:
                seqdb_update_seqsource(seqdb_ws, seqdb_id, seqdb_gene_region_id)

            if logging.getLogger().isEnabledFor(logging.DEBUG):
                # retrieve inserted record and display to users for validation
                # purposes
                result = seqdb_ws.getJsonSeq(seqdb_id)
                tools_helper.pretty_log_json(
                    result, level="debug", message="Record as inserted:")
        else:
            logging.info("Skipping gi:%s, which does not contain a sequence." %(genbank_id))

    elif result['count'] == 1:
        logging.info(
            "Sequence for gi:%s already exists in SeqDB. Skipping." % (
                genbank_id))

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
    main_conf = tools_helper.load_config('../config.yaml')
    tool_config = tools_helper.load_config('seqdb_gb_insert_config.yaml')
    logging.config.dictConfig(main_conf['logging'])
    http_client.HTTPConnection.debuglevel = main_conf[
        'http_connect']['debug_level']


    # caching the entrez records shaved 2 minutes off the time to load 
    # ~740 sequences from query: "(*DAOM*[source] and levesque and not 'unplaced genomics scaffold')"
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
    feature_type_lookup = {}

    logging.info(
        "Script executed with the following command and arguments: %s" % (
            sys.argv))

    seqdb_ws = seqdbWebService(
        tool_config['seqdb']['apikey'], main_conf['seqdb']['url'])

    Entrez.email = tool_config['entrez']['email']
    query = tool_config['entrez']['query']

    logging.info("Querying GenBank: \'%s\'" % (tool_config['entrez']['query']))

    # preliminary query to find out how many records there are
    record = entrez_search(query)

    # setup loop counters; retrieving records 50 at a time
    count = int(record["Count"])
    start = 0
    retrieve = 50
    logging.info(
        "Query returned %i records. Retrieving them in batches of %i" % (
            count, retrieve))

    # repeat until we have all records
    while start < count:

        # retrieve block of records
        logging.debug("Retrieving %i..%i" % (start, start + retrieve))
        record = entrez_search(query, retstart=start, retmax=retrieve)

        # process each returned id in the batch
        for genbank_id in record["IdList"]:
            process_entrez_entry(seqdb_ws, genbank_id, cache=entrez_cache, lookup=feature_type_lookup)

        start += retrieve

if __name__ == "__main__":
    main()

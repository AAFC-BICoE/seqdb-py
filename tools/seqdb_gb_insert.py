#! /usr/bin/env python
"""Script to use Entrez and SeqDB-WS APIs to load GenBank sequences into SeqDB
"""

import sys
import json
import yaml
import urllib
import logging.config
import httplib as http_client

from api.seqdbWebService import seqdbWebService
from Bio import Entrez

# Note: Can't log here - statements would occur before logger config


def load_config():
    """Return a config object populated from 'config.yaml'"""
    try:
        config = yaml.load(file('config.yaml', 'r'))
        return config
    except yaml.YAMLError, exc:
        print "Error in configuration file:", exc

# Note: Not logging method call as method is a logging helper


def pretty_log_json(obj, level="info", message=None):
    """Pretty print an object as a JSON formatted str to the log

    Args:
        obj: Object to pretty log

    Kargs:
        level (str): One of "debug", "info", "warn", "error".  Default="info".

    Returns:
        None

    Raises:
        None

    Example usage:
    >>> pretty_log_json({"a":1})
    >>> pretty_log_json({"a":1}, level="debug")
    >>> pretty_log_json({"a":1}, message="Contents of obj:")
    >>> pretty_log_json({"a":1}, level="debug", message="Contents of obj:")
    """
    display = ""
    if message is not None:
        display = display + message + "\n"
    display = display + json.dumps(obj, sort_keys=True, indent=4)
    getattr(logging, level)(display)


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
                genes[qualifier["GBQualifier_value"]] = 1
                logging.debug("Gene name: %s" %
                              (qualifier["GBQualifier_value"]))
    logging.info("Found %i gene names" % (len(genes.keys())))
    return genes.keys()


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
        genbank_id, rettype="gb", retmode="xml", database="nucleotide"):
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
    handle = Entrez.efetch(
        db=database, id=genbank_id, rettype=rettype, retmode=retmode)
    record = Entrez.read(handle)
    handle.close()
    pretty_log_json(
        record[0], level="debug", message="Retrieved Genbank Record: ")
    return record


def seqdb_ret_entrez_gene_region_id(seqdb_ws, genbank_id, record):
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
        logging.warn(
            "GenBank sequence 'gi:%s' contains multiple gene regions. "
            "Currently we only assign sequences to a single gene region. "
            "This sequence will not be assigned to a gene region." % (
                genbank_id))
    elif len(genes) == 1:
        seqdb_gene_region_id = check_region(seqdb_ws, genes[0], create=True)
        logging.debug("Found gene: %s, SeqDB region id: %i" %
                      (genes[0], seqdb_gene_region_id))

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
    pretty_log_json(dict_for_logging, level="debug",
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


def process_entrez_entry(seqdb_ws, genbank_id):
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
    """
    logging.info("Processing gi:%s" % (genbank_id))

    # Ensure the record doesn't already exist in SeqDB
    # If it does, continue with the next record
    logging.debug("Checking for gi:%s in SeqDB" % (genbank_id))

    result = seqdb_ws.getJsonConsensusSequenceIdsByGI(genbank_id)

    if result['count'] == 0:
        # retrieve genbank record
        record = entrez_fetch(genbank_id)

        seqdb_id = seqdb_insert_entrez_sequence(seqdb_ws, genbank_id, record)

        seqdb_gene_region_id = seqdb_ret_entrez_gene_region_id(
            seqdb_ws, genbank_id, record)

        if seqdb_gene_region_id is not None:
            seqdb_update_seqsource(seqdb_ws, seqdb_id, seqdb_gene_region_id)

        if logging.getLogger().isEnabledFor(logging.DEBUG):
            # retrieve inserted record and display to users for validation
            # purposes
            result = seqdb_ws.getJsonSeq(seqdb_id)
            pretty_log_json(
                result, level="debug", message="Record as inserted:")

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
    config = load_config()
    logging.config.dictConfig(config['logging'])
    http_client.HTTPConnection.debuglevel = config[
        'http_connect']['debug_level']

    logging.info(
        "Script executed with the following command and arguments: %s" % (
            sys.argv))

    seqdb_ws = seqdbWebService(
        config['seqdb']['apikey'], config['seqdb']['url'])

    Entrez.email = config['entrez']['email']
    query = config['entrez']['query']

    logging.info("Querying GenBank: \'%s\'" % (config['entrez']['query']))

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
            process_entrez_entry(seqdb_ws, genbank_id)

        start += retrieve

if __name__ == "__main__":
    main()

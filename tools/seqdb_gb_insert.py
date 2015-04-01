#! /usr/bin/env python

import sys, requests, json, yaml, urllib
import logging, logging.config
import httplib as http_client

from api.seqdbWebService import seqdbWebService
from Bio import Entrez

# Note: Can't log here - statements would occur before logger config
def load_config():
    try:
        config = yaml.load(file('config.yaml', 'r'))
        return config
    except yaml.YAMLError, exc:
        print "Error in configuration file:", exc

# Note: Not logging method call as method is a logging helper
def pretty_log_json(j, level="info", message=None):
    display = ""
    if message != None:
        display = display + message + "\n"
    display = display + json.dumps(j, sort_keys=True, indent=4)
    log = getattr(logging, level)(display)

def format_sequence_name(genbankId, record):
    logging.info("Formating sequence name for gi:%s" % (genbankId))
#    name = "gi|" + str(genbankId) + "|"
#    name = name + "gb|" + record[0]["GBSeq_accession-version"] + "|" 
#    name = name + record[0]["GBSeq_definition"].replace(", ", "_").replace("; ", "_").replace(" ", "_").replace(";","").replace(",","")
#    return name
    name = record[0]["GBSeq_definition"]
    logging.info("Returning: %s" % (name))
    return name

def format_tracefile_name(**keywds):
    logging.info("Formatting tracefile name")
    result = "?" + urllib.urlencode(keywds)
    logging.info("Returning: %s" % (result))
    return result
    
def extract_gene_names(record):
    logging.info("Extracting gene names from record: %s" % (record["GBSeq_accession-version"]))
    genes = {}
    for feature in record["GBSeq_feature-table"]:
        if feature["GBFeature_key"] == "gene":
            for qualifier in feature["GBFeature_quals"]:
                genes[qualifier["GBQualifier_value"]] = 1
                logging.debug("Gene name: %s" % (qualifier["GBQualifier_value"]))
    logging.info("Found %i gene names" % (len(genes.keys())))
    return genes.keys()

def check_region(seqdbWS, gene, create=False):
    logging.info("Checking SeqDB for region: %s.  Create == %r" %(gene, create))
    region_id = None
    region_ids = seqdbWS.getRegionIdsByName(gene)
    if len(region_ids) == 0 and create == True:
        region_id = seqdbWS.createRegion(gene, "GenBank Gene: %s" %(gene))
        logging.debug("Created region: %i in seqdb for %s" % (region_id, gene))
    elif len(region_ids) == 1:
        region_id = region_ids[0]
        logging.debug("Found region: %i in seqdb for %s" % (region_id, gene))
    else:
        logging.warn("Found multiple regions for '%s' in SeqDB. " \
                     "Currently unable to assign sequences to regions with non-unique names." % (gene))

    logging.info("Returning region id: %r" % (region_id))
    return region_id

def main(arv):

    config = load_config()
    logging.config.dictConfig(config['logging'])
    http_client.HTTPConnection.debuglevel = config['http_connect']['debug_level']

    logging.info("Script executed with the following command and arguments: %s" % sys.argv)

    seqdbWS = seqdbWebService(config['seqdb']['apikey'], config['seqdb']['url'])

    Entrez.email = config['entrez']['email']
    query = config['entrez']['query']

    logging.info("Querying GenBank: \'%s\'" % (query))
    # preliminary query to find out how many records there are
    handle = Entrez.esearch(db="nucleotide", retmax=10, term=query)
    record = Entrez.read(handle)
    handle.close()

    # setup loop counters; retrieving records 50 at a time
    count = int(record["Count"])
    start = 0;
    retrieve = 50;
    logging.info("Query returned %i records.  Retrieving them in batches of %i" % (count, retrieve))

    # repeat until we have all records
    while start < count:

        # retrieve block of records
        logging.debug("Retrieving %i..%i" % (start, start+retrieve))
        handle = Entrez.esearch(db="nucleotide", retstart = start, retmax=retrieve, term=query)
        record = Entrez.read(handle)
        handle.close()

        # for each returned id
        for genbankId in record["IdList"]:

            logging.debug("Processing gi:%s" % (genbankId))

            # Ensure the record doesn't already exist in SeqDB
            # If it does, continue with the next record
            logging.debug("Checking for gi:%s in SeqDB" % (genbankId))
            r = seqdbWS.getJSONConsensusSequenceIdsByGI(genbankId)
            if r['count'] == 1:
                logging.info("Sequence for gi:%s already exists in SeqDB. Skipping." % (genbankId))
                continue

            # retrieve genbank record
            logging.info("Retrieving record for gi:%s from GenBank" % (genbankId))
            handle = Entrez.efetch(db="nucleotide", id=genbankId, rettype="gb", retmode="xml")
            record = Entrez.read(handle)
            handle.close()
            pretty_log_json(record[0], level="debug", message="Retrieved Genbank Record: ")

            seqdb_gene_region = None
            genes = extract_gene_names(record[0])
            if len(genes) > 1:
                logging.warn("GenBank sequence 'gi:%i' contains multiple gene regions. "\
                             "Currently we only assign sequences to a single gene region. "\
                             "This sequence will not be assigned to a gene region." % (genbankId))
            elif len(genes) == 1:
                seqdb_gene_region = check_region(seqdbWS, genes[0], create=True)
                logging.debug("Found gene: %s, SeqDB region id: %i" % (genes[0], seqdb_gene_region))

            seq_name = format_sequence_name(genbankId, record)
            sequence = record[0]["GBSeq_sequence"]
            tracefileDir = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' 
            tracefileName = format_tracefile_name(db="nucleotide", id=genbankId, rettype="gb", retmode="xml")
            additional = {
                    'chromatDir': tracefileDir,
                    'chromatName': tracefileName,
                    'genBankGI': genbankId,
                    'genBankAccession': record[0]["GBSeq_primary-accession"], 
                    'genBankVersion': record[0]["GBSeq_accession-version"],
                    'submittedToInsdc': 'true'
                    }
            tmp = {'name': seq_name, 'sequence': sequence}
            tmp.update(additional)
            logging.info("Adding sequence for gi:%s to SeqDB" % (genbankId))
            pretty_log_json(tmp, level="debug", message="Creating consensus (non-default values): ")
            seqdb_id,code,message = seqdbWS.createConsensusSequence(seq_name, sequence, additional=additional)
            logging.info("Created consensus seqdbid:%i for gi:%s (%s), Status: %i, Message: %s" % (seqdb_id, genbankId, record[0]["GBSeq_accession-version"], code, message))

            if seqdb_gene_region != None:
                logging.info("Updating SeqSource for consensus sequence seqdbid:%i" % (seqdb_id))
                seqsource = {
                        "seqSource": {
                            "region": {
                                "id": seqdb_gene_region,
                                }
                            }
                        }
                region_id,code,message = seqdbWS.updateSeqSource(seqdb_id, seqsource)
                logging.info("Updated SeqSource for consensus sequence seqdbid:%i, Status: %i, Message: %s" % (seqdb_id, code, message))

            # retrieve inserted record and display to users for validation purposes
            # TODO only do this if in debug mode
            r = seqdbWS.getJSONSeq(seqdb_id)
            pretty_log_json(r, level="debug", message="Record as inserted:")

        start += retrieve

if __name__ == "__main__":
        main(sys.argv[1:])

#! /usr/bin/env python

import sys, requests, json, yaml, urllib
from api.seqdbWebService import seqdbWebService
from Bio import Entrez

def load_config():
    try:
        config = yaml.load(file('config.yaml', 'r'))
        return config
    except yaml.YAMLError, exc:
        print "Error in configuration file:", exc

def pretty_print_json(j, message=None):
    if message != None:
        print message
    print json.dumps(j, sort_keys=True, indent=4)

def format_sequence_name(genbankId, record):
#    name = "gi|" + str(genbankId) + "|"
#    name = name + "gb|" + record[0]["GBSeq_accession-version"] + "|" 
#    name = name + record[0]["GBSeq_definition"].replace(", ", "_").replace("; ", "_").replace(" ", "_").replace(";","").replace(",","")
#    return name
    return record[0]["GBSeq_definition"]

def format_tracefile_name(**keywds):
    return "?" + urllib.urlencode(keywds)
    
def main(arv):
    config = load_config()
    seqdbWS = seqdbWebService(config['seqdb']['apikey'], config['seqdb']['url'])

    Entrez.email = config['entrez']['email']
    query = config['entrez']['query']

    # preliminary query to find out how many records there are
    handle = Entrez.esearch(db="nucleotide", retmax=10, term=query)
    record = Entrez.read(handle)
    handle.close()

    # setup loop counters; retrieving records 50 at a time
    count = int(record["Count"])
    start = 0;
    retrieve = 50;

    # repeat until we have all records
    while start < count:

        # retrieve block of records
        handle = Entrez.esearch(db="nucleotide", retstart = start, retmax=retrieve, term=query)
        record = Entrez.read(handle)
        handle.close()

        # for each returned id
        for genbankId in record["IdList"]:

            # Ensure the record doesn't already exist in SeqDB
            # If it does, continue with the next record
            r = seqdbWS.getJSONConsensusSequenceIdsByGI(genbankId)
            if r['count'] == 1:
                print "Sequence for gi: %s already exists in SeqDB. Skipping." % (genbankId)
                continue

            # retrieve genbank record
            handle = Entrez.efetch(db="nucleotide", id=genbankId, rettype="gb", retmode="xml")
            record = Entrez.read(handle)
            handle.close()
            pretty_print_json(record[0], message="Retrieved Genbank Record: ")

            seq_name = format_sequence_name(genbankId, record)
            sequence = record[0]["GBSeq_sequence"]
            # Treating base URL as dir
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
            pretty_print_json(tmp, "Creating consensus (non-default values): ")
            seqdb_id,code,message = seqdbWS.createConsensusSequence(seq_name, sequence, additional=additional)
            print "Create consensus: Id: %i, Status: %i, Message: %s" % (seqdb_id, code, message)

            # retrieve inserted record and display to users for validation purposes
            r = seqdbWS.getJSONSeq(seqdb_id)
            pretty_print_json(r, message="Inserted record:")

        start += retrieve

if __name__ == "__main__":
        main(sys.argv[1:])

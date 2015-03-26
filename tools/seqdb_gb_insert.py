#! /usr/bin/env python

import sys, requests, json, yaml
from api.seqdbWebService import seqdbWebService
from Bio import Entrez

def load_config():
    try:
        config = yaml.load(file('config.yaml', 'r'))
        return config
    except yaml.YAMLError, exc:
        print "Error in configuration file:", exc

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

            # retrieve genank record in the returned list of ids
            handle = Entrez.efetch(db="nucleotide", id=genbankId, rettype="gb", retmode="xml")
            record = Entrez.read(handle)
            handle.close()
#           print(json.dumps(record[0], sort_keys=True, indent=4))

            # Try to extract the strain name if there is a feature table
            # The strainName seems to generally contain the collection code (e.g. DAOM_XXXXXX)
            strainName = ""
            if record[0]['GBSeq_feature-table'] == []:
                print "no feature-table \n"
            else:
                featuresTable = record[0]['GBSeq_feature-table']
                features = featuresTable[0]['GBFeature_quals']
                for feature in features:
                    qualName = feature['GBQualifier_name']
                    qualValue = feature['GBQualifier_value']
                    if qualName == 'strain':
                        strainName = qualValue.replace(" ", "_").replace(";","")                    

            # display record to users for validation purposes
            print "gi: " + genbankId + "\n" + \
                "accession: " + record[0]["GBSeq_primary-accession"] + "\n" + \
                "version: " + record[0]["GBSeq_accession-version"] + "\n" + \
                "description: " + record[0]["GBSeq_definition"] + "\n" + \
                "sequence: " + record[0]["GBSeq_sequence"] + "\n" + \
                "organism: " + record[0]["GBSeq_organism"].replace(" ", "_") + "\n" + \
                "strain: " + strainName + "\n"
            

            seq_name = "gi:" + str(genbankId) + "|" + record[0]["GBSeq_primary-accession"] + "|" + record[0]["GBSeq_organism"].replace(" ", "_") + "_" + strainName
            sequence = record[0]["GBSeq_sequence"]
            additional = {
                    'genBankGI': genbankId,
                    'genBankAccession': record[0]["GBSeq_primary-accession"], 
                    'genBankVersion': record[0]["GBSeq_accession-version"],
                    'submittedToInsdc': 'true'
                    }
            seqdb_id,code,message = seqdbWS.createConsensusSequence(seq_name, sequence, additional=additional)
            print "Create consensus: Id: %i, Status: %i, Message: %s" % (seqdb_id, code, message)

            # retrieve inserted record and display to users for validation purposes
            r = seqdbWS.getJSONSeq(seqdb_id)
            print json.dumps(r)

            # could also/instead retrive by id
            #    extract id of record and retrieve directly

            # delete the sequence after adding it to cleanup
            # not tested
            #   seqdbWS.deleteConsensusSequence(seqdb_id)

        start += retrieve

if __name__ == "__main__":
        main(sys.argv[1:])

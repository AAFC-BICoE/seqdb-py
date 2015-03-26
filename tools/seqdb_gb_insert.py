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

def form_url (config, obj):
    return config['seqdb']['url'] + "/" + obj

def main(arv):
    config = load_config()
    seqdbWS = seqdbWebService(config['seqdb']['apikey'], config['seqdb']['url'])
    Entrez.email = config['entrez']['email']
    query = config['entrez']['query']

    seqdb_headers_update = {
            'Content-type':'application/json',
            'apikey': config['seqdb']['apikey']
        }

    seqdb_headers_read = {
            'apikey': config['seqdb']['apikey']
        }

    handle = Entrez.esearch(db="nucleotide", retmax=10, term=query)
    record = Entrez.read(handle)
    handle.close()

    count = int(record["Count"])
    start = 0;
    retrieve = 50;

    while start < count:
        handle = Entrez.esearch(db="nucleotide", retstart = start, retmax=retrieve, term=query)
        record = Entrez.read(handle)
        handle.close()
        for id in record["IdList"]:
            handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="xml")
            record = Entrez.read(handle)
            handle.close()
#           print(json.dumps(record[0], sort_keys=True, indent=4))
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
            print "gi: " + id + "\n" + \
                "accession: " + record[0]["GBSeq_primary-accession"] + "\n" + \
                "version: " + record[0]["GBSeq_accession-version"] + "\n" + \
                "description: " + record[0]["GBSeq_definition"] + "\n" + \
                "sequence: " + record[0]["GBSeq_sequence"] + "\n" + \
                "organism: " + record[0]["GBSeq_organism"].replace(" ", "_") + "\n" + \
                "strain: " + strainName + "\n"
            

            name = "gi:" + str(id) + "|" + record[0]["GBSeq_primary-accession"] + "|" + record[0]["GBSeq_organism"].replace(" ", "_") + "_" + strainName
            sequence = record[0]["GBSeq_sequence"]
            additional = {
                    'genBankGI': id,
                    'genBankAccession': record[0]["GBSeq_primary-accession"], 
                    'genBankVersion': record[0]["GBSeq_accession-version"],
                    'submittedToInsdc': 'true'
                    }
            r = seqdbWS.createConsensusSequence( name, sequence, additional=additional )
            print r

            r = seqdbWS.getJSONConsensusSequencesByName(record[0]["GBSeq_primary-accession"])
            print json.dumps(r)

            # could also/instead retrive by id
            #    extract id of record and retrieve directly

            # delete the sequence after adding it to cleanup
            # not tested
            #   seqdbWS.deleteConsensusSequence(seqdb_id)

        start += retrieve

if __name__ == "__main__":
        main(sys.argv[1:])

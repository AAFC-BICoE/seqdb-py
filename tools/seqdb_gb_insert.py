#! /usr/bin/env python

import sys, requests, json, yaml
from Bio import Entrez

def load_config():
        try:
                config = yaml.load(file('config.yaml', 'r'))
                return config
        except yaml.YAMLError, exc:
                print "Error in configuration file:", exc

def retrieve_sequence(url, headers, name):
	params = {'filterName':'sequence.name',
		  'filterValue':name,
		  'filterOperator':'and',
		  'filterWildcard':'true'}
	r = requests.get(url, headers=headers, params=params)
	print r
	print r.text
	print json.dumps(json.loads(r.text))
	return 0

def delete_sequence(url, headers, seqdb_id):
	url = url + "/" + str(seqdb_id)
	print "delete: " + url
	#requests.delete(url, headers)

def submit_sequence(url, headers, gi, name, accession, version, sequence):
        payload = {'consensus': 
			{'name': name, 
			 'seq': sequence, 
			 'seqType': "N", 
			 'readNum': 0,
			 'genBankAccession':accession,
			 'genBankVersion':version,
			 'genBankGI':gi,
			 'submittedToInsdc': 'true'
			}
		}
        r = requests.post(url, data=json.dumps(payload), headers=headers)
        print r

def form_url (config, obj):
	return config['seqdb']['url'] + "/" + obj

def main(arv):
	config = load_config()
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
#			print(json.dumps(record[0], sort_keys=True, indent=4))
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
			submit_sequence(form_url(config, "consensus"), seqdb_headers_update, 
						id,
						"gi:" + str(id) + "|" + record[0]["GBSeq_primary-accession"] + "|" + record[0]["GBSeq_organism"].replace(" ", "_") + "_" + strainName, 
						record[0]["GBSeq_primary-accession"], 
						record[0]["GBSeq_accession-version"],
						record[0]["GBSeq_sequence"],
					)
			seqdb_id = retrieve_sequence(form_url(config, "consensus"), seqdb_headers_read, record[0]["GBSeq_primary-accession"])
			#if seqdb_id > 0:
			#	delete_sequence(form_url(config, "sequence"), seqdb_headers_update, seqdb_id)

		start += retrieve

if __name__ == "__main__":
        main(sys.argv[1:])

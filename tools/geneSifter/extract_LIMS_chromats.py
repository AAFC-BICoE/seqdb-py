'''
Created on Apr 14, 2015

@author: korolo

This module is designed to be used by the GeneSifter LIMS system.
The system has the ability to trigger actions based on events. This scipt is an action in response
to the completed run event. It extracts chromatogram blobls from the LIMS db, stores them in a 
directory, runs a script to import them to seqdb and cleans up.
'''

import sys, requests
import logging.config
import shutil
import psycopg2
import yaml
from api.seqdbWebService import seqdbWebService, UnexpectedContent

#log_file_name = "/tmp/extract_LIMS_chromats.log"
api_blob_file_name = 'blob_api.ab1.gz'
db_blob_file_name = 'blob_db.ab1.gz'

def write_simple_log(msg):
    with open("test.log", "a") as log_file:
        log_file.write(msg + "\n")
    log_file.close();
    


def load_config():
    """Return a config object populated from 'config.yaml'"""
    try:
        config = yaml.load(file('config.yaml', 'r'))
        return config
    except yaml.YAMLError, exc:
        print "Error in configuration file:", exc
        

def get_chromat_from_api(api_url, api_key, query_id):        
    params = {'key': api_key, 'id': query_id}
    resp = requests.get(api_url, params=params)
    
    logging.info("SQL API response heaeder: \n %s" %resp.headers)

    return resp.content

def write_blob_to_file(blob, file_name):
    output_file = open(file_name, 'wb')
    output_file.write(blob)
    output_file.close()

    logging.info("Data written to a file %s" % file_name)
    '''
    # Writing raw without decoding (same file as regular write)
    r = requests.get(api_url, params=params, stream=True)
    if r.status_code == 200:
        r.raw.decode_content = False
        with open("blob_raw.ab1.gz", 'wb') as f:
            shutil.copyfileobj(r.raw, f)  
            
    '''

def get_chromat_from_db(db_host, db_name, db_user, db_pass, db_query):
    try:
        db_connection_str = "host='%s' dbname='%s' user='%s' password='%s'" % (
             db_host, db_name, db_user, db_pass)
        
        conn = psycopg2.connect(db_connection_str)
    except:
        logging.error("Unable to connect to the database '%s' at '%s'" % (db_name, db_host))
        sys.exit(1)
    
    cursor = conn.cursor()
    logging.info("Executing LIMS db query: %s" % db_query)
    cursor.execute(db_query)
    blob = cursor.fetchone()
    logging.info("Received LIMS db query response of type '%s' and size %s" % (type(blob), len(blob)))

    return blob[0]    

def test_seqdb_api(base_url, api_key):
    seqdbWS = seqdbWebService(api_key, base_url)
    
    # Get all feature types from seqDB
    seqdb_feat_types = seqdbWS.getFeatureTypesWithIds()
    logging.info("Got SeqDB regions back! Number of regions: %s" % len(seqdb_feat_types))
    write_simple_log("Script: Got SeqDB regions back! Number of regions: %s" % len(seqdb_feat_types))


def main():
    """Load completed LIMS chromats into SeqDB."""
    
    write_simple_log("Script: in main().")
    
    config = load_config()
    logging.config.dictConfig(config['logging'])
    
    write_simple_log("Script: got yaml config.")
   
    logging.info("Script executed with the following command and arguments: %s" % (sys.argv))

    '''
    # Retrieving chromatograms through SQL API does not work (yet). The file comes back corrupted.
    blob = get_chromat_from_api(
                                api_key = config['lims']['sql_api']['key'], 
                                api_url = config['lims']['sql_api']['url'], 
                                query_id = config['lims']['query_ids']['get_sample_chromat'])

    write_blob_to_file(blob, config['chromat_files_dir'] + api_blob_file_name)
    '''

    blob = get_chromat_from_db(
                               db_host=config['lims']['db']['host'], 
                               db_name = config['lims']['db']['name'], 
                               db_user = config['lims']['db']['user'], 
                               db_pass = config['lims']['db']['password'], 
                               db_query = config['lims']['db']['sample_chromat_query'])
    
    write_blob_to_file(blob, config['chromat_files_dir'] + db_blob_file_name)
    
    write_simple_log("Script: Sample chromat was written to a file %s" % config['chromat_files_dir'] + db_blob_file_name)
    print "Sample chromat was written to a file %s" % db_blob_file_name
    
    test_seqdb_api(base_url=config['seqdb']['api']['url'], api_key=config['seqdb']['api']['key'])

    #open(db_blob_file_name, 'wb').write(blob[0])
      

if __name__ == '__main__':    
    write_simple_log("Script: begin.")
    main()
    write_simple_log("Script: end.")
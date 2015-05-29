'''
Created on Feb 12, 2015

@author: korolo

Usage:
pull_seqdb_its_seqs -c <Path to yaml config with SeqDB API info>
or
pull_seqdb_its_seqs -k <SeqDB API key> -u <SeqDB API URL> 
Other arguments:
   -h   help (prints this message)
'''
import sys 
import os
import getopt
import requests.exceptions
import logging.config
import tools_helper
from api.seqdbWebService import seqdbWebService, UnexpectedContent


usage_help_line = """Usage of the script: \npull_seqdb_its_seqs -c <Path to yaml config with SeqDB API info>
or
pull_seqdb_its_seqs -k <SeqDB API key> -u <SeqDB API URL> 
Other arguments:
   -h   help (prints this message)
"""

# File name where the pulled sequences will be stored
output_file_name = "seqdb_ITS_seqs.fasta"
# This log will provide users of Galaxy with extra information on the tool execution
# sysem statements should not go here, since full log is configured in yaml
user_log = tools_helper.SimpleLog("seqdb_pull.log")



# Parses command line arguments 
# Returns seqdb api_key and base url to use for web services requests
def parse_input_args(argv):
    ''' Parses command line arguments
    Returns:
        config_file: path to a config file with has api information, 
            or empty string if no such usage 
        seqdb api_key to use for web services requests
    '''
    config_file=''
    api_url=''
    api_key = ''
    
    try:
        opts, args = getopt.getopt(argv,"hc:k:u:",["config_file=", "seqdb_api_key=", "seqdb_ws_url="])
    except getopt.GetoptError:
        print usage_help_line
        logging.error("Invalid script arguments. Aborting.")
        user_log.error("Invalid script arguments. Aborting.")
        sys.exit(2)
        
    if len(opts)==0:
        print(usage_help_line)
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print(usage_help_line)
            sys.exit(2)
        elif opt in ("-c", "--config_file"):
            config_file = arg
        elif opt in ("-k", "--seqdb_api_key"):
            api_key = arg
        elif opt in ("-u", "--seqdb_ws_url"):
            api_url = arg
    
    if config_file and api_key and api_url:
        print("Supply either a configuration file or API url with key.")
        print(usage_help_line)
        logging.error("Invalid script arguments. Aborting.")
        user_log.error("Invalid script arguments. Aborting.")
        sys.exit(2)
    
    if bool(api_key) != bool(api_url):
        print("Both API url and API key have to be supplied.")
        print(usage_help_line)
        logging.error("Invalid script arguments. Aborting.")
        user_log.error("Invalid script arguments. Aborting.")
        sys.exit(2)
        
    return (config_file, api_url, api_key)




def pull_its_seqs(api_key,base_url):
    ''' Downloads all ITS sequences from SeqDB and writes them to a file '''
    
    
    seqdbWS = seqdbWebService(api_key, base_url)
    
    try:
        its_region_ids = seqdbWS.getItsRegionIds()
    except requests.exceptions.ConnectionError as e:
        user_log.error("Could not connect to Sequence DB. Contact your sysadmin for more details.\n")
        logging.error("Could not connect to Sequence DB. ")
        logging.error(e.message)
        sys.exit("Script error.")
    except requests.exceptions.ReadTimeout as e:
        user_log.error("Connection too slow for getting data from Sequence DB. \Contact your sysadmin for more details.")
        logging.error("Connection too slow for getting data from Sequence DB. ")
        logging.error(e.message)
        sys.exit("Script error.")
    except requests.exceptions.HTTPError as e:
        user_log.error("HTTP error when getting region ids from Sequence DB. Contact your sysadmin for more details.")
        logging.error("HTTP error when getting region ids from Sequence DB. ")
        logging.error(e.message)
        sys.exit("Script error.")
    except UnexpectedContent as e:
        user_log.error("Sequence DB API response did not match expected format. Contact your sysadmin for more details.")
        logging.error("Sequence DB API response did not match expected format.")
        logging.error(e.message)
        sys.exit("Script error.")
    except Exception as e:
        user_log.error("Script encountered an error. See your sysadmin for details.")
        logging.error(e.message)
        sys.exit("Script error.")
        
    logging.info("Number of ITS regions retrieved: %i " % len(its_region_ids))
    user_log.info("Number of ITS regions retrieved: %i " % len(its_region_ids))

   
    #Get sequence IDs for the ITS regions
    its_seq_ids = []
    for region_id in its_region_ids:
        try:
            curr_seq_ids = seqdbWS.getSeqIds(region_id)
            its_seq_ids.extend(curr_seq_ids)
        except requests.exceptions.ConnectionError as e:
            user_log.error("Could not connect to Sequence DB. Contact your sysadmin for more details.\n")
            logging.error("Could not connect to Sequence DB. ")
            logging.error(e.message)
            sys.exit("Script error.")
        except requests.exceptions.ReadTimeout as e:
            user_log.error("Connection too slow for getting data from Sequence DB. Contact your sysadmin for more details.")
            logging.error("Connection too slow for getting data from Sequence DB. ")
            logging.error(e.message)
            sys.exit("Script error.")
        except requests.exceptions.HTTPError as e:
            user_log.error("HTTP error when getting region ids from Sequence DB. Contact your sysadmin for more details.")
            logging.error("HTTP error when getting region ids from Sequence DB.")
            logging.error(e.message)
            sys.exit("Script error.")
        except UnexpectedContent as e:
            user_log.error("Sequence DB API response did not match expected format. Contact your sysadmin for more details.")
            logging.error("Sequence DB API response did not match expected format.")
            logging.error(e.message)
            sys.exit("Script error.")
                        

    logging.info("Number of ITS sequences retrieved: %i " % len(its_seq_ids))
    user_log.info("Number of ITS sequences retrieved: %i " % len(its_seq_ids))
     
    # Get fasta sequences based on ids and write to a file 
    output_file = open(output_file_name, 'w')
    
    success_ids = []
    for seq_id in its_seq_ids:
        try:
            # Request sequence in fasto format from SeqDB:
            fastaSequence = seqdbWS.getFastaSeq(seq_id)
            output_file.write(fastaSequence)
            success_ids.append(seq_id)
        except requests.exceptions.ConnectionError as e:
            user_log.error("Could not connect to Sequence Database. Contact your sysadmin for more details.")
            logging.error("Could not connect to Sequence DB.")
            logging.error(e.message)
            sys.exit("Script error.")
        except requests.exceptions.ReadTimeout as e:
            user_log.error("Connection too slow for getting data from Sequence DB. Contact your sysadmin for more details.")
            logging.error("Connection too slow for getting data from Sequence DB. ")
            logging.error(e.message)
            sys.exit("Script error.")
        except requests.exceptions.HTTPError as e:
            user_log.error("HTTP error when getting region ids from Sequence DB. Contact your sysadmin for more details.")
            logging.error("HTTP error when getting region ids from Sequence DB.")
            logging.error(e.message)
            sys.exit("Script error.")
        except UnexpectedContent as e:
            user_log.error("Sequence DB API response did not match expected format. Contact your sysadmin for more details.")
            logging.error("Sequence DB API response did not match expected format.")
            logging.error(e.message)
            sys.exit("Script error.")
     
    output_file.close()   

    logging.info("Number of ITS sequences successfully written to a file: %s" % len(success_ids) )
    user_log.info("Number of ITS sequences successfully written to a file: %s" % len(success_ids) )
    
    logging.info("ITS sequences written to a file: %s" % os.path.abspath(output_file.name))
    user_log.info("ITS sequences written to a file: %s" % output_file.name)

    return success_ids
    

def main():
    main_conf = tools_helper.load_config('../config.yaml')
    logging.config.dictConfig(main_conf['logging'])
    
    logging.info("Script executed with the following command and arguments: %s" % sys.argv)
    user_log.info("Script execution started")
    
    # Parse command line to get seqdb api key (necessary to request seqDB web services) and base url for web services requests
    config_file, api_url, api_key = parse_input_args(sys.argv[1:])
    
    if config_file:
        config = tools_helper.load_config(config_file)
        api_url = config['seqdb']['api_url'] 
        api_key = config['seqdb']['api_key'] 
    
    logging.info("Base URL for web services is: '%s'" % api_url)
    user_log.info("Base URL for web services is: '%s'" % api_url)
   
    success_seq_ids = pull_its_seqs(api_key, api_url)
    

    print("Number of sequences retrieved from Sequence Dababase:  %s" % len(success_seq_ids)) 
    print("Sequences are written to a file: '%s'" % output_file_name)
    print("Execution log is written to a file: '%s'" % user_log.getFileName())
    print("Execution complete.")

    user_log.info("Execution complete.")
    user_log.close()
    
    logging.info("Script execution ended.")


if __name__ == '__main__':
    main()
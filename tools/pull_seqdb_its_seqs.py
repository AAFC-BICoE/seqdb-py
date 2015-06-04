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
import logging.config
import requests.exceptions
import tools_helper
from api.seqdbWebService import seqdbWebService, UnexpectedContent
from config import config_root


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
        seqdb api_url to use for web services requests
    '''
    config_file=''
    api_url=''
    api_key = ''
    
    try:
        opts, args = getopt.getopt(argv,"hc:k:u:",["config_file=", "seqdb_api_key=", "seqdb_api_url="])
    except getopt.GetoptError:
        print usage_help_line
        logging.error(tools_helper.log_msg_argError)
        user_log.error(tools_helper.log_msg_argError)
        sys.exit(tools_helper.log_msg_sysExit)
        
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
        elif opt in ("-u", "--seqdb_api_url"):
            api_url = arg
    
    if config_file and api_key and api_url:
        print(tools_helper.log_msg_argErrorConfigFileUrl)
        print(usage_help_line)
        logging.error(tools_helper.log_msg_argError)
        user_log.error(tools_helper.log_msg_argError)
        sys.exit(2)
    
    if bool(api_key) != bool(api_url):
        print(tools_helper.log_msg_argErrorKeyUrl)
        print(usage_help_line)
        logging.error(tools_helper.log_msg_argError)
        user_log.error(tools_helper.log_msg_argError)
        sys.exit(2)
        
    return (config_file, api_url, api_key)




def pull_its_seqs(api_key,base_url):
    ''' Downloads all ITS sequences from SeqDB and writes them to a file '''
    
    
    seqdbWS = seqdbWebService(api_key, base_url)
    
    try:
        its_region_ids = seqdbWS.getItsRegionIds()
    except requests.exceptions.ConnectionError as e:
        user_log.error("%s %s" % (tools_helper.log_msg_noDbConnection, tools_helper.log_msg_sysAdmin))
        logging.error(tools_helper.log_msg_noDbConnection)
        logging.error(e.message)
        sys.exit(tools_helper.log_msg_sysExit)
    except requests.exceptions.ReadTimeout as e:
        user_log.error("%s %s" % (tools_helper.log_msg_slowConnection, tools_helper.log_msg_sysAdmin))
        logging.error(tools_helper.log_msg_slowConnection)
        logging.error(e.message)
        sys.exit(tools_helper.log_msg_sysExit)
    except requests.exceptions.HTTPError as e:
        user_log.error("%s %s" % (tools_helper.log_msg_httpError, tools_helper.log_msg_sysAdmin))
        logging.error(tools_helper.log_msg_httpError)
        logging.error(e.message)
        sys.exit(tools_helper.log_msg_sysExit)
    except UnexpectedContent as e:
        user_log.error("%s %s" % (tools_helper.log_msg_apiResponseFormat, tools_helper.log_msg_sysAdmin))
        logging.error(tools_helper.log_msg_apiResponseFormat)
        logging.error(e.message)
        sys.exit(tools_helper.log_msg_sysExit)
    except Exception as e:
        user_log.error("%s %s" % (tools_helper.log_msg_scriptError, tools_helper.log_msg_sysAdmin))
        logging.error(e.message)
        sys.exit(tools_helper.log_msg_sysExit)
        
    logging.info("Number of ITS regions retrieved: %i " % len(its_region_ids))
    user_log.info("Number of ITS regions retrieved: %i " % len(its_region_ids))

   
    #Get sequence IDs for the ITS regions
    its_seq_ids = []
    for region_id in its_region_ids:
        try:
            curr_seq_ids = seqdbWS.getSeqIds(region_id)
            its_seq_ids.extend(curr_seq_ids)
        except requests.exceptions.ConnectionError as e:
            user_log.error("%s %s" % (tools_helper.log_msg_noDbConnection, tools_helper.log_msg_sysAdmin))
            logging.error(tools_helper.log_msg_noDbConnection)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except requests.exceptions.ReadTimeout as e:
            user_log.error("%s %s" % (tools_helper.log_msg_slowConnection, tools_helper.log_msg_sysAdmin))
            logging.error(tools_helper.log_msg_slowConnection)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except requests.exceptions.HTTPError as e:
            user_log.error("%s %s" % (tools_helper.log_msg_httpError, tools_helper.log_msg_sysAdmin))
            logging.error(tools_helper.log_msg_httpError)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except UnexpectedContent as e:
            user_log.error("%s %s" % (tools_helper.log_msg_apiResponseFormat, tools_helper.log_msg_sysAdmin))
            logging.error(tools_helper.log_msg_apiResponseFormat)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except Exception as e:
            user_log.error("%s %s" % (tools_helper.log_msg_scriptError, tools_helper.log_msg_sysAdmin))
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
                        

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
    ''' Retrieves ITS sequenes from SeqDB '''
    main_conf = tools_helper.load_config(config_root.path() + '/config.yaml')

    if not main_conf:
        logging.error(tools_helper.log_msg_noConfig)
        sys.exit(tools_helper.log_msg_sysExit)
    
    logging.config.dictConfig(main_conf['logging'])
    
    logging.info("Script executed with the following command and arguments: %s" % sys.argv)
    user_log.info(tools_helper.log_msg_execStarted_simple)
    
    config_file, api_url, api_key = parse_input_args(sys.argv[1:])
    
    if config_file:
        tool_config = tools_helper.load_config(config_file)
        api_url = tool_config['seqdb']['api_url'] 
        api_key = tool_config['seqdb']['api_key'] 
    
    logging.info("%s '%s'" % (tools_helper.log_msg_apiUrl, api_url))
    user_log.info("%s '%s'" %  (tools_helper.log_msg_apiUrl, api_url))
   
    success_seq_ids = pull_its_seqs(api_key, api_url)
    

    print("Number of sequences retrieved from Sequence Dababase:  %s" % len(success_seq_ids)) 
    print("Sequences are written to a file: '%s'" % output_file_name)
    print("Execution log is written to a file: '%s'" % user_log.getFileName())
    print("Execution complete.")

    user_log.info(tools_helper.log_msg_execEnded)
    user_log.close()
    
    logging.info(tools_helper.log_msg_execEnded)


if __name__ == '__main__':
    main()
'''
Created on Feb 12, 2015

@author: korolo

Usage:
python pull_seqdb_its_seqs.py -k <SeqDB API key>
or
python pull_seqdb_its_seqs.py --seqdb_api_key <SeqDB API key>
-h is help
-p is to use production URL (uses ***REMOVED*** (Oksana's) machine by default, since we don't have UAT for SeqDB as of this time)
'''
import sys, getopt, requests, logging
import tools_helper
from api.seqdbWebService import seqdbWebService, UnexpectedContent
#from fileinput import filename
#from pip._vendor.distlib._backport.tarfile import filemode

usage_help_line = """Usage of the script: \npull_seqdb_its_seqs -c <Path to yaml config with SeqDB API info>
or
pull_seqdb_its_seqs -k <SeqDB API key> -u <SeqDB API URL> 
Other arguments:
   -h   help (prints this message)
"""

prod_url = "***REMOVED***/api/v1"
local_url = "***REMOVED***:2002/seqdb/api/v1"
# file name where the received sequences will be stores
output_file_name = "seqdb_ITS_seqs.fasta"
log_file_name = "seqdb_pull.log"



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
        opts, args = getopt.getopt(argv,"hk:u:",["seqdb_api_key=", "seqdb_ws_url="])
    except getopt.GetoptError:
        print usage_help_line
        sys.exit(2)
        
    if len(opts)==0:
        print(usage_help_line)
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print usage_help_line
            sys.exit()
        elif opt in ("-c", "--config_file"):
            config_file = arg
        elif opt in ("-k", "--seqdb_api_key"):
            api_key = arg
        elif opt in ("-u", "--seqdb_ws_url"):
            api_url = arg
            
    
    return (config_file, api_url, api_key)




def pull_its_seqs(api_key,base_url):
    ''' Downloads all ITS sequences from SeqDB and writes them to a file '''
    
    
    seqdbWS = seqdbWebService(api_key, base_url)
    
    try:
        its_region_ids = seqdbWS.getItsRegionIds()
    except requests.exceptions.ConnectionError as e:
        logging.error("Could not connect to Sequence DB. \n")
        logging.error(e.message)
        print("Could not connect to Sequence DB. See log file for details.")
        sys.exit(1)
    except requests.exceptions.ReadTimeout as e:
        logging.error("Connection too slow for getting data from Sequence DB. \n")
        logging.error(e.message)
        print "Connection too slow for getting data from Sequence DB. See log file for details."
        sys.exit(1)
    except requests.exceptions.HTTPError as e:
        logging.error("HTTP error when getting region ids from Sequence DB. \n")
        logging.error(e.message)
        print "HTTP error when getting region ids from Sequence DB. See log file for details."
        sys.exit(1)
    except UnexpectedContent as e:
        print e
        sys.exit(1)
        
    logging.info("Number of ITS regions retrieved: %i \n" % len(its_region_ids))

   
    #Get sequence IDs for the ITS regions
    its_seq_ids = []
    for region_id in its_region_ids:
        try:
            curr_seq_ids = seqdbWS.getSeqIds(region_id)
            its_seq_ids.extend(curr_seq_ids)
        except requests.exceptions.ConnectionError as e:
            logging.error("Could not connect to Sequence DB. \n")
            logging.error(e.message)
            print "Could not connect to Sequence DB. See log file for details."
            sys.exit(1)
        except requests.exceptions.ReadTimeout as e:
            logging.error("Connection too slow for getting data from Sequence DB. \n")
            logging.error(e.message)
            print "Connection too slow for getting data from Sequence DB. See log file for details."
            sys.exit(1)
        except requests.exceptions.HTTPError as e:
            logging.error("HTTP error when getting region ids from Sequence DB. \n")
            logging.error(e.message)
            print "HTTP error when getting region ids from Sequence DB. See log file for details."
            sys.exit(1)
        except UnexpectedContent as e:
            print e
            sys.exit(1)
            

    logging.info("Number of ITS sequences retrieved: %i \n" % len(its_seq_ids))
     
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
            logging.error("Could not connect to Sequence DB. \n")
            logging.error(e.message)
            print "Could not connect to Sequence DB. See log file for details."
            sys.exit(1)
        except requests.exceptions.ReadTimeout as e:
            logging.error("Connection too slow for getting data from Sequence DB. \n")
            logging.error(e.message)
            print "Connection too slow for getting data from Sequence DB. See log file for details."
            sys.exit(1)
        except requests.exceptions.HTTPError as e:
            logging.error("HTTP error when getting region ids from Sequence DB. \n")
            logging.error(e.message)
            print "HTTP error when getting region ids from Sequence DB. See log file for details."
            sys.exit(1)
        except UnexpectedContent as e:
            print e
            sys.exit(1)
     
    output_file.close()   
    
    logging.info("Number of ITS sequences written to a file (%s): %s" % (output_file.name, len(success_ids)) )
    
    return success_ids
    

def main():
    # Start a log file. filemode='w' overwrites the log for each program run
    logging.basicConfig(filename=log_file_name, filemode='w', level=logging.DEBUG)
    
    logging.info("Script executed with the following command and arguments: %s" % sys.argv)
    
    # Parse command line to get seqdb api key (necessary to request seqDB web services) and base url for web services requests
    config_file, api_url, api_key = parse_input_args(sys.argv[1:])
    
    if config_file:
        config = tools_helper.load_config(config_file)
        api_url = config['seqdb']['api_url'] 
        api_key = config['seqdb']['api_key'] 
    
    logging.info("Base URL for web services is: '%s'" % api_url)
   
    success_seq_ids = pull_its_seqs(api_key, api_url)
    
    print "Execution complete."
    print "Number of sequences loaded from Sequence Dababase:  %s" % len(success_seq_ids) 
    print "Sequences written to a file: '%s'" % output_file_name
    print "Execution log is written to a file: '%s'" % log_file_name



if __name__ == '__main__':
    main()
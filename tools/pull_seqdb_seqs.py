'''
Created on Feb 12, 2015

@author: korolo

Usage:
pull_seqdb_its_seqs -c <Path to yaml config with SeqDB API info> 
or
pull_seqdb_its_seqs -k <SeqDB API key> -u <SeqDB API URL> 
with one of the following options:
 --its
 --consensus
 --raw
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

import argparse

usage_help_line = """Usage of the script: \npull_seqdb_seqs -c <Path to yaml config with SeqDB API info>
or
pull_seqdb_its_seqs -k <SeqDB API key> -u <SeqDB API URL> 
  with one of the following options:
  --its
  --consensus
  --raw

Other arguments:
   -h   help (prints this message)
"""

### Values below are used in Galaxy wrappers, so make sure you know what 
### you're doing if you're changing any of them 

# File name where the pulled sequences will be stored. 
output_file_name = "seqdb_sequences.fasta"
# This log will provide users of Galaxy with extra information on the tool 
# execution sysem statements should not go here, since full log is configured
# in yaml
user_log = tools_helper.SimpleLog("seqdb_pull.log")
# Values for the types of sequences this script downloads. I.e. "its" loads 
# ITS sequences.
pull_types_dict = {"its":"its", "consensus":"consensus", "raw":"raw"}

###


# Parses command line arguments 
# Returns seqdb api_key and base url to use for web services requests
def parse_input_args(argv):
    ''' Parses command line arguments '''
    
    pull_types_set = frozenset(pull_types_dict.values())
    
    parser = argparse.ArgumentParser(description="Load sequences from SeqDB")
    parser.add_argument('seq_type', help="Type of sequences to load", type=str, choices=pull_types_set)
    parser.add_argument('-c', help="SeqDB config file", dest="config_file", required=False)
    parser.add_argument('-u', help="SeqDB API URL", dest="api_url", required=False)
    parser.add_argument('-k', help="SeqDB API key", dest="api_key", required=False)    
    parser.add_argument('-s', help="Specimen name", dest="specimen_id", required=False)    
    #parser.add_argument('-t', help="Type of sequences to load", dest="load_type", type=str, choices=set(("its","consensus")), required=True)
    
    args = parser.parse_args(argv)

    if not (args.config_file or (args.api_url and args.api_key)):
        parser.error('Either -c <configuration file>, or -u <api_url> -k <api_key> have to be specified')
    
    if args.seq_type == pull_types_dict["raw"] and not args.specimen_id:
        parser.error('To load raw sequences, please specify specimen id: -s <specimen_id>')
        
    return args
    

def __init__(self, api_url, api_key):
    self.api_url = api_url
    self.api_key = api_key 


def get_ITS_seq_ids(seqdbWS):
    ''' Get all sequence ids, which are associated with the ITS regions '''
    
    ### Get ITS regions
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
    logging.info("ITS region ids retrieved: %s " % its_region_ids)
   
    ### Get sequence IDs for the ITS regions
    its_seq_ids = []
    for region_id in its_region_ids:
        try:
            curr_seq_ids = seqdbWS.getSequenceIdsByRegion(region_id)
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
    
    msg_numITSseqs = "Number of ITS sequences retrieved:"
    logging.info("%s %i " % (msg_numITSseqs, len(its_seq_ids)))
    user_log.info("%s %i " % (msg_numITSseqs, len(its_seq_ids)))

    return its_seq_ids


def get_consensus_seq_ids(seqdbWS):
    ''' Get all SeqDB consensus sequence ids (accessible with this API key) '''
    try:
        consensus_seq_ids = seqdbWS.getConsensusSequenceIds()
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
        
    msg_numConsensusSeqs = "Number of consensus sequences retrieved:"
    logging.info("%s %i " % (msg_numConsensusSeqs, len(consensus_seq_ids)))
    user_log.info("%s %i " % (msg_numConsensusSeqs, len(consensus_seq_ids)))

    return consensus_seq_ids
    
    
def get_raw_seq_ids(seqdbWS, specimen_id):
    ''' Returns raw sequence ids for a specimen '''
    
    try:
        seq_ids = seqdbWS.getSequenceIdsBySpecimen(specimen_id)
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
        
    logging.info("Number of raw sequences retrieved for specimen '%s': %i " % (specimen_id, len(seq_ids)))
    user_log.info("Number of raw sequences retrieved for specimen '%s': %i " % (specimen_id, len(seq_ids)))
    
    return seq_ids
    
         
def write_fasta_file(seqdbWS, its_seq_ids, fasta_file_name):
    # Get fasta sequences based on ids and write to a file 
    output_file = open(fasta_file_name, 'w')
    
    success_ids = []
    for seq_id in its_seq_ids:
        try:
            # Request sequence in fasto format from SeqDB:
            fastaSequence = seqdbWS.getFastaSeq(seq_id)
            output_file.write(fastaSequence)
            success_ids.append(seq_id)
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
     
    output_file.close()   

    msg_fileName = "Sequences written to a file:"
    logging.info("%s %s" % (msg_fileName, os.path.abspath(output_file.name)))
    user_log.info("%s %s" % (msg_fileName, os.path.abspath(output_file.name)))

    msg_seqNum = "Number of sequences written:"
    logging.info("%s %s" % (msg_seqNum, len(success_ids)) )
    user_log.info("%s %s" % (msg_seqNum, len(success_ids)) )
    

    return success_ids
    

    
def main():
    ''' Retrieves ITS sequenes from SeqDB '''
    
    ### Load main configuration file and set up logging for the script
    
    main_conf = tools_helper.load_config(config_root.path() + '/config.yaml')

    if not main_conf:
        logging.error(tools_helper.log_msg_noConfig)
        sys.exit(tools_helper.log_msg_sysExit)
    
    logging.config.dictConfig(main_conf['logging'])
    
    logging.info("%s %s" % (tools_helper.log_msg_scriptExecutionWithParams, sys.argv))
    user_log.info(tools_helper.log_msg_execStarted_simple)
    
    ### Parse sript's input arguments
    
    parsed_args = parse_input_args(sys.argv[1:])
    
    if parsed_args.config_file:
        tool_config = tools_helper.load_config(parsed_args.config_file)
        api_url = tool_config['seqdb']['api_url'] 
        api_key = tool_config['seqdb']['api_key'] 
    
    logging.info("%s '%s'" % (tools_helper.log_msg_apiUrl, api_url))
    user_log.info("%s '%s'" %  (tools_helper.log_msg_apiUrl, api_url))
    
    ### Script execution
    
    seqdbWS = seqdbWebService(api_key, api_url)
     
    if pull_types_dict["its"] == parsed_args.seq_type:
        logging.info(tools_helper.log_msg_ITSLoad)
        user_log.info(tools_helper.log_msg_ITSLoad)
        
        seq_ids = get_ITS_seq_ids(seqdbWS)
        
    elif pull_types_dict["consensus"] == parsed_args.seq_type:
        logging.info(tools_helper.log_msg_ConsensusLoad)
        user_log.info(tools_helper.log_msg_ConsensusLoad)
        
        seq_ids = get_consensus_seq_ids(seqdbWS)
        
    elif pull_types_dict["raw"] == parsed_args.seq_type:
        logging.info("%s %s" % (tools_helper.log_msg_RawLoad, parsed_args.specimen_id))
        user_log.info("%s %s" % (tools_helper.log_msg_RawLoad, parsed_args.specimen_id))
        
        seq_ids = get_raw_seq_ids(seqdbWS, parsed_args.specimen_id)
        
    success_seq_ids = write_fasta_file(seqdbWS, seq_ids, output_file_name)

    
    ### Post-execution: messages and logging
    
    print("Number of sequences retrieved from Sequence Dababase:  %s" % len(success_seq_ids)) 
    print("Sequences are written to a file: '%s'" % output_file_name)
    print("Execution log is written to a file: '%s'" % user_log.getFileName())
    print("Execution complete.")

    user_log.info(tools_helper.log_msg_execEnded)
    user_log.close()
    
    logging.info(tools_helper.log_msg_execEnded)
    

if __name__ == '__main__':
    main()
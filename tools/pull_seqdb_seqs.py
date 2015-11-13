'''
Created on Feb 12, 2015

@author: korolo

Given a SeqDB URL and API key, uses SeqDB webservices to load sequences, filtered by parameters 
below, into a fasta file.

Usage (for most current usage run this script with no options: >python pull_seqdb_seqs.py):
pull_seqdb_seqs -c <Path to yaml config with SeqDB API info> 
or
pull_seqdb_its_seqs -k <SeqDB API key> -u <SeqDB API URL> 
with one of the following options:
 --its
 --consensus
 --raw
Other arguments:
   -h   help (prints this message)
'''
import argparse
import logging.config
import os
import sys 

import requests.exceptions

from api.seqdbWebService import seqdbWebService, UnexpectedContent
from config import config_root
import tools_helper



### Values below are used in Galaxy wrappers, so make sure you know what 
### you're doing if you're changing any of them 
# File name where the pulled sequences will be stored. 
output_file_name = "seqdb_sequences.fasta"
# File name where taxonomy for the sequences will be stored. Optional output. 
output_taxonomy_file_name = "seqdb_sequences_taxonomy.fasta"
# This log will provide users of Galaxy with extra information on the tool 
# execution sysem statements should not go here, since full log is configured
# in yaml
user_log = tools_helper.SimpleLog("seqdb_pull.log")
# Values for the types of sequences this script downloads. I.e. "its" loads 
# ITS sequences. Note that raw sequences are not implemented in SeqDB yet. "raw":"raw",
pull_types_dict = {"its":"its", "consensus":"consensus", "all":"all"}
# Taxonomy ranks that can be specified as a filter parameter with corresponding value
# These are used as a drop-down value in the wrapper
taxonomy_ranks = {"species", "genus", "family", "order", "class", "phylum"}

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
    parser.add_argument('-t', help="Output taxonomy file as well as fasta", dest="output_taxonomy_file", action='store_true', required=False)
    parser.add_argument('--specNum', help="Specimen number (identifier)", dest="specimen_num", required=False)    
    parser.add_argument('--seqName', help="Sequence name (keyword)", dest="sequence_name", required=False)   
    parser.add_argument('--geneRegion', help="Gene region name (keyword)", dest="gene_region_name", required=False)    
    parser.add_argument('--project', help="Project Name (keyword)", dest="project_name", required=False)    
    parser.add_argument('--collectionCode', help="Collection code (keyword)", dest="collection_code", required=False)    
    parser.add_argument('--pubRefSeqs', help="Public reference sequences", dest="pub_ref_seqs", action='store_true', required=False)    
    parser.add_argument('--taxRank', help="Taxonomic rank to filter sequences on (need to specify the value as well). Ex. --taxRank phylum --taxValue chordata ", dest="tax_rank", choices=taxonomy_ranks, required=False)
    parser.add_argument('--taxValue', help="Value for the taxonomic rank to be filtered on (need to specify the rank as well)", dest="tax_value", required=False)
    #parser.add_argument('-t', help="Type of sequences to load", dest="load_type", type=str, choices=set(("its","consensus")), required=True)
    
    args = parser.parse_args(argv)

    if not (args.config_file or (args.api_url and args.api_key)):
        parser.error('Either -c <configuration file>, or -u <api_url> -k <api_key> have to be specified')
    
    if args.seq_type == pull_types_dict["its"] and (args.specimen_num or args.sequence_name or args.pub_ref_seqs):
        parser.error('ITS sequences can not be restricted by filters at the moment. Please do not use --specNum, --seqName or any other additional options.')
    
    if bool(args.tax_rank) != bool(args.tax_value):
        parser.error('Either both --taxRank and --taxValue have to be specified, or none.')
    
        
    return args
    

def __init__(self, api_url, api_key):
    self.api_url = api_url
    self.api_key = api_key 


def get_ITS_seq_ids(seqdbWS):
    ''' Get all sequence ids, which are associated with the ITS regions '''
    
    ### Get sequence IDs for the ITS regions
    its_seq_ids = set()
    its_region_names = ["18s", "its", "28s",  "ssu", "16s", "5.8s", "lsu", "23s", "25s", "internal transcribed spacer"]
    for its_region_keyword in its_region_names:
        #TODO: parallelize; use locking when appending to its_seq_ids
        try:
            curr_seq_ids, offset = seqdbWS.getSequenceIdsWithOffset(regionName=its_region_keyword)
            its_seq_ids.update(curr_seq_ids)
            while offset:
                curr_seq_ids, offset = seqdbWS.getSequenceIdsWithOffset(regionName=its_region_keyword, offset=offset)
                its_seq_ids.update(curr_seq_ids)
                
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

    return list(its_seq_ids)

    
def get_seq_ids(seqdbWS, pull_type,
                specimen_num=None, 
                sequence_name=None, 
                pub_ref_seqs=None, 
                region_name=None, 
                project_name=None,
                collection_code=None,
                taxonomy_rank=None, taxonomy_value=None):
    ''' Gets sequence ids based on specified parameters 
    Agrs:
        pull_type: string of pre-determined values. Values should correspond to the values of pull_types_dict
        specimen_num: if specified, specimen number for which the sequence ids will be retrieved
    '''
    
    if pull_type not in pull_types_dict.values():
        msg = "Value for pull_type should be one of the following: %s" %pull_types_dict.values()
        logging.error(msg)
        sys.exit(tools_helper.log_msg_sysExit + msg)
    
    if pull_type == pull_types_dict["its"]:
        seq_ids = get_ITS_seq_ids(seqdbWS)
    else:
        try:
            if pull_type == pull_types_dict["consensus"]:
                
                seq_ids = seqdbWS.getAllConsensusSequenceIds(specimenNum=specimen_num, 
                                                          sequenceName=sequence_name, 
                                                          pubRefSeq=pub_ref_seqs,
                                                          regionName=region_name,
                                                          projectName=project_name,
                                                          collectionCode=collection_code,
                                                          taxonomyRank=taxonomy_rank, 
                                                          taxonomyValue=taxonomy_value)
                    
                log_msg = "Number of consensus sequences retrieved:"
            elif pull_type == pull_types_dict["all"]:
                
                seq_ids = seqdbWS.getAllSequenceIds(specimenNum=specimen_num, 
                                                sequenceName=sequence_name,
                                                pubRefSeq=pub_ref_seqs,
                                                regionName=region_name,
                                                collectionCode=collection_code)
                
                log_msg = "Number of sequences retrieved:"
            elif pull_type == pull_types_dict["raw"]:
                sys.exit("Raw sequence retrieval is not implemented yet.")
                #seq_ids = seqdbWS.getRawSequenceIds()
                #log_msg = "Number of raw sequences retrieved:"
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
        
        logging.info("%s %i " % (log_msg, len(seq_ids)))
        user_log.info("%s %i " % (log_msg, len(seq_ids)))
        
    return seq_ids
    
         
def write_fasta_file(seqdbWS, its_seq_ids, fasta_file_name):
    # Get fasta sequences based on ids and write to a file 
    output_file = open(fasta_file_name, 'w')
    
    success_ids = []
    for seq_id in its_seq_ids:
        #TODO: threads?
        '''
        import threading
        
        t = threading.Thread(target=<this try code>, args(seq_id,)
        '''
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
    

def write_taxonomy_file(seqdbWS, seq_ids, output_file_name):
    # Get fasta sequences based on ids and write to a file 
    output_file = open(output_file_name, 'w')
    
    success_ids = []
    for seq_id in seq_ids:
        
        try:
            #TODO get a taxonomy line from API
            taxonomy_line = "lala"
            output_file.write(taxonomy_line)
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
    else:
        api_url = parsed_args.api_url 
        api_key = parsed_args.api_key 
        
    logging.info("%s '%s'" % (tools_helper.log_msg_apiUrl, api_url))
    user_log.info("%s '%s'" %  (tools_helper.log_msg_apiUrl, api_url))
    
    
    ### Script execution
    
    seqdbWS = seqdbWebService(api_key, api_url)
     
    if pull_types_dict["its"] == parsed_args.seq_type:
        logging.info(tools_helper.log_msg_ITSLoad)
        user_log.info(tools_helper.log_msg_ITSLoad)
        
        seq_ids = get_ITS_seq_ids(seqdbWS)
    else:
        log_msg = "Loading %s sequences." %parsed_args.seq_type
        logging.info(log_msg)
        user_log.info(log_msg)

        seq_ids = get_seq_ids(seqdbWS=seqdbWS, 
                              pull_type=parsed_args.seq_type, 
                              specimen_num=parsed_args.specimen_num,
                              sequence_name=parsed_args.sequence_name,
                              region_name=parsed_args.gene_region_name,
                              project_name=parsed_args.project_name,
                              collection_code=parsed_args.collection_code,
                              pub_ref_seqs=parsed_args.pub_ref_seqs,
                              taxonomy_rank=parsed_args.tax_rank,
                              taxonomy_value=parsed_args.tax_value)

    success_seq_ids = write_fasta_file(seqdbWS, seq_ids, output_file_name)
    if (parsed_args.output_taxonomy_file):
        write_taxonomy_file(seqdbWS, seq_ids, output_taxonomy_file_name)

    
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
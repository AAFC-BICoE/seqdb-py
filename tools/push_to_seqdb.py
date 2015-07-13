'''
Created on July 8, 2015

@author: korolo

Usage:
push_stuff_to_seqdb -c <Path to yaml config with SeqDB API info> 
or
pull_seqdb_its_seqs -k <SeqDB API key> -u <SeqDB API URL> 
with one of the following options:
 -itsx_features
 -findLCA_taxonomy
Other arguments:
   -h   help (prints this message)
'''
import sys 
import os
import logging.config
import requests.exceptions
import tools_helper
from api.seqdbWebService import seqdbWebService, UnexpectedContent
from config import config_root
from TaxonomyLineage import TaxonomyLineage

import argparse
import time

usage_help_line = """Usage of the script: \npush_stuff_to_seqdb -c <Path to yaml config with SeqDB API info> 
or
pull_seqdb_its_seqs -k <SeqDB API key> -u <SeqDB API URL> 

with one of the following options:
 -itsx_features
 -findLCA_taxonomy

Other arguments:
   -h   help (prints this message)
"""

### Values below are used in Galaxy wrappers, so make sure you know what 
### you're doing if you're changing any of them 

# This log will provide users of Galaxy with extra information on the tool 
# execution sysem statements should not go here, since full log is configured
# in yaml
user_log = tools_helper.SimpleLog("seqdb_push.log")
# Values for the types of sequences this script downloads. I.e. "its" loads 
# ITS sequences. Note that raw sequences are not implemented in SeqDB yet. "raw":"raw",
push_types_dict = {"its":"itsx_features", "taxonomy":"findLCA_taxonomy"}

###



def parse_input_args(argv):
    ''' Parses command line arguments '''
    
    push_types_set = frozenset(push_types_dict.values())
    
    parser = argparse.ArgumentParser(description="Push information to SeqDB")
    parser.add_argument('push_type', help="Type of information to write", type=str, choices=push_types_set)
    parser.add_argument('-c', help="SeqDB config file", dest="config_file", required=False)
    parser.add_argument('-u', help="SeqDB API URL", dest="api_url", required=False)
    parser.add_argument('-k', help="SeqDB API key", dest="api_key", required=False)    
    parser.add_argument('-ifile', help="File with information to write to SeqDB", dest="info_file", required=True)    
    
    args = parser.parse_args(argv)

    if not (args.config_file or (args.api_url and args.api_key)):
        parser.error('Either -c <configuration file>, or -u <api_url> -k <api_key> have to be specified')
        
    return args
    

def open_file(file_name):
    ''' Opens a file and writes error/log messages if not successful 
    Raises:
        IOError
    '''
    file_handler = None
    try:
        file_handler = open(file_name,"r")
    except IOError as e:
        if e.errno == 2:
            error_msg = "Could not open file: %s." % file_name
            user_log.error(error_msg)
            logging.error(error_msg)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        else:
            raise

    return file_handler

def log_error(error_msg):
    ''' Logs the message to user and developer logs'''
    user_log.error(error_msg)
    logging.error(error_msg)
    
def report_log_error(error_msg):
    ''' Logs the message to user and developer logs and prints the message to stdout'''
    log_error(error_msg)
    print error_msg

def log_info(msg):
    ''' Logs the message to user and developer logs'''
    user_log.info(msg)
    logging.info(msg)
    
def report_log_info(msg):
    ''' Logs the message to user and developer logs and prints the message to stdout'''
    log_info(msg)
    print msg


def get_lieage_taxids(tax_parent_ids, taxon_id, lineage=None):
    ''' Finds a taxonomic lineage for taxon_id.
    Return: list of taxid in lineage order, from given tax id up. 
    '''
    

def push_taxonomy_data(seqdbWS, info_file_handler):
    ''' Extracts taxon id from the findLCA output file, finds lineage for the taxon id
        and pushes taxonomic identification to SeqDB
    '''
    
    input_file_line_example = """ Query: seqdb|6726    LCA: 129355  Name:Phytophthora lateralis    Rank: species    Matches: 1    Names: Phytophthora lateralis:1,  Evalue: 0.0,    Pid: 98.62,    Qcover: 0.867647058823529,    SourceGI: 320336337, 
    or
    Query: seqdb|6802    No suitable matches found."""
    
    """
    Query: seqdb|6802    No suitable matches found.
    """
    """
    Query: seqdb|28954    LCA: 164328    Name:Phytophthora ramorum    Rank: species    Matches: 1    Names: Phytophthora ramorum:55,    Evalue: 0.0,    Pid: 96.97,    Qcover: 0.929705215419501,    SourceGI: 46812520,74476198,78192394,81176726,81176727,81176728,408690088,378750406,378750408,643192938,643192939,44194092,47934204,49617503,60308861,89275889,108885436,114146744,55586083,643192936,55586084,46812521,42521138,32490530,315419600,643192937,55586085,643192935,227810703,227810722,46399123,323301609,227810723,295409544,284432200,308745533,551031928,227810464,227810545,227810570,372125558,284432218,284432396,308745534,478739030,227810622,227810639,42529489,227810623,284432190,227810547,417357152,417357145,38348758,171262896,
    """
    """
    Query: seqdb|64071    LCA: 4783    Name:Phytophthora    Rank: genus    Matches: 4    Names: Phytophthora ramorum:2,Phytophthora mengei:1,Phytophthora citricola:1,Phytophthora capsici:2,    Evalue: 0.0,    Pid: 95.00,95.15,    Qcover: 0.928864569083447,    SourceGI: 320336471,308913225,308913137,308912963,308912935,320336197,
    """
    # parse the file and extract seqdb sequence ID - ncbi taxon id pairs
    seqdb_sequence_taxon = {}
    for line in info_file_handler:
        error_msg = "%s Found line: \n%s\n Example of an expected line: \n%s" %(tools_helper.log_msg_wrongFileFormat, line, input_file_line_example)
        
        line_tokens = line.split('\t')
        
        # Do basic file format verification:
        if len(line_tokens) != 10 and len(line_tokens) != 2:
            report_log_error(tools_helper.log_msg_wrongFileFormat)
            sys.exit(tools_helper.log_msg_sysExitFile)

        try:
            token1 = line_tokens[0].split(': ')
            if token1[0] != "Query":
                report_log_error(error_msg)
                sys.exit(tools_helper.log_msg_sysExitFile)

            sequenceId = token1[1].split("|")[1]
            sequenceId = int(sequenceId)

            if ": " in line_tokens[1]:
                lca = line_tokens[1].split(": ")[1]
                lca = int(lca)
                seqdb_sequence_taxon[sequenceId] = lca
            else:
                # There were no matched found for this sequence id, so save the message from a file
                seqdb_sequence_taxon[sequenceId] = line_tokens[1].rstrip()
                
            
        except:
            report_log_error(error_msg)
            sys.exit(tools_helper.log_msg_sysExitFile)
        
    
    print seqdb_sequence_taxon
    
    # Find lineage for each seqdb id
  
    start_time = time.time()

    tax_lineage = TaxonomyLineage("./data/ncbi_taxonomy/")
    lineage_names = tax_lineage.findLineage(129355)
        
    print lineage_names  
    print "Time to get lineage new way: %s" %(time.time()-start_time)
    
    # Write each lineage to SeqDB
    
    # Return taxonomy table id / feature id for each identification written


# B15_17_SH817_ITS_ITS5    622 bp.    SSU: Not found    ITS1: 1-241    5.8S: 242-399    ITS2: 400-557    LSU: 558-622
def push_its_features(seqdbWS, features_file):
    ''' Extracts found ITS features from ITSx tools results .positions.txt file
        and writes these features to SeqDB
        
    Args:
        seqdbWS: api.seqdbWebService object with accessor methods to Seq
        features_file: file handler with opened .positions.txt file
        
    Return:
        list of SeqDB feature_ids for those features that were successfully written to SeqDB
    '''
    
    # The name of the feature type for each itsx feature in SeqDB
    its_feature_type_name = "misc_RNA"  # agreen on convention
    rRna_feature_type_name = "rRNA"     # agreed on convention
    feature_description = "Generated by ITSx tool."
    output_file_name = 'seqdb_feature_ids.txt'
    input_file_line_example = "seqdb|1685    321 bp.    SSU: 1-46    ITS1: 47-283    5.8S: No end    ITS2: Not found    LSU: Not found    Broken or partial sequence, only partial 5.8S!"

    
    created_feature_ids = []
    
   
    # Get all feature types from seqDB
    seqdb_feat_types = seqdbWS.getFeatureTypesWithIds()
    
    its_feature_type_id = ''
    rRna_feature_type_id = ''
    # Create feature type in SeqDB, if don't exist already
    if its_feature_type_name not in seqdb_feat_types:
        its_feature_type_id = seqdbWS.createFeatureType(its_feature_type_name)
    else:
        its_feature_type_id = seqdb_feat_types[its_feature_type_name]
        
    if rRna_feature_type_name not in seqdb_feat_types:
        rRna_feature_type_id = seqdbWS.createFeatureType(rRna_feature_type_name)
    else:
        rRna_feature_type_id = seqdb_feat_types[rRna_feature_type_name]
    
            
    ### TODO: Finalize what to use for the frame: and strand:
    ###        (using default 1 and 1 for now)
    
    # parse the file and insert each found feature to SeqDB
    for line in features_file:
        line_tokens = line.split('\t')

        try:
            sequenceId = line_tokens[0].split('|')[1]
            sequenceId = int(sequenceId)
        except:
            error_msg = "Input file error. Could not extract sequence id from the input file. Example of an expected line:\n%s\n" % input_file_line_example
            report_log_error(error_msg)
            sys.exit(tools_helper.log_msg_sysExit)
        
        try:
            itsx_features = line_tokens[2:7]
        except IndexError as e:
            error_msg = "Could not extract ITS features from the input file. Example of an expected line:\n%s\n" % input_file_line_example
            report_log_error(error_msg)
            sys.exit(tools_helper.log_msg_sysExit)
        
        if not sequenceId or not itsx_features:
            error_msg = "Input file not in the correct format. Example of an expected line:\n%s\n" % input_file_line_example
            report_log_error(error_msg)
            sys.exit(tools_helper.log_msg_sysExit)
            
                
        for itsx_feature_token in itsx_features:
            feature_location_pair = itsx_feature_token.split(": ")
            try:
                location = feature_location_pair[1].split("-")
                location = map(int, location)
                location = [{"start":location[0],"end":location[1],"frame":1,"strand":1}]
                feature_type_id = its_feature_type_id
                if feature_location_pair[0] in {"SSU", "5.8S", "LSU"}:
                    feature_type_id = rRna_feature_type_id
                fid = seqdbWS.insertFeature(feature_location_pair[0], feature_type_id, location, sequenceId, description=feature_description)
                created_feature_ids.append(fid)
            except requests.exceptions.ConnectionError as e:
                user_log.error("%s %s" % (tools_helper.log_msg_noDbConnection, tools_helper.log_msg_sysAdmin))
                logging.error(tools_helper.log_msg_noDbConnection)
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
            except:
                warning_msg = "File token '%s' is not in the expected format of <feature name>:<position>. Ignoring." % itsx_feature_token
                user_log.info(warning_msg)
                logging.info(warning_msg)


    # Write ids of the inserted features into a file
    output_file = open(output_file_name, 'w')
    for fid in created_feature_ids:
        output_file.write(str(fid) + '\n')
    output_file.close()
  
    log_msg1 = "Number of features written to Sequence Dababase:   %i " % len(created_feature_ids)
    log_msg2 = "Created feature ids are written to a file: '%s'" % output_file_name
    report_log_info(log_msg1)
    report_log_info(log_msg2)
        
    return created_feature_ids
        

    
def main():
    ''' Write provided information to SeqDB '''
    
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
    
    log_info("%s '%s'" % (tools_helper.log_msg_apiUrl, api_url))
    
    msg_inputFile = "File with information to write to SeqDB: %s" % parsed_args.info_file
    log_info(msg_inputFile)
    
    ### Script execution
    
    seqdbWS = seqdbWebService(api_key, api_url)
    
    # Open input file
    info_file_handler = ''
    try:
        info_file_handler = open(parsed_args.info_file,"r")
    except IOError as e:
        if e.errno == 2:
            error_msg = "Could not open input positions file '%s'." % parsed_args.info_file
            log_error(error_msg)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        else:
            raise
             
    if push_types_dict["its"] == parsed_args.push_type:
        #sys.exit("Not yet implemented")
        log_msg = "Writing ITS features to SeqDB."
        log_info(log_msg) 
        success_feat_ids = push_its_features(seqdbWS, info_file_handler)
            
        
    elif push_types_dict["taxonomy"] == parsed_args.push_type:
        log_msg = "Writing taxonomy lineage information to SeqDB."
        log_info(log_msg) 
        
        push_taxonomy_data(seqdbWS, info_file_handler)


    
    ### Post-execution: messages and logging
    
    
    print(tools_helper.log_msg_execEnded)
    print("Execution log is written to a file: '%s'" % user_log.getFileName())

    log_info(tools_helper.log_msg_execEnded)
    user_log.close()
    
    

if __name__ == '__main__':
    main()
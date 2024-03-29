'''
Created on July 8, 2015

@author: korolo

Usage (for the most updated usage, run "python push_to_seqdb.py -h"):
push_stuff_to_seqdb -c <Path to yaml config with SeqDB API info> 
or
pull_seqdb_its_seqs -k <SeqDB API key> -u <SeqDB API URL> 
with one of the following options:
 -itsx_features
 -findLCA_taxonomy
Other arguments:
   -h   help (prints this message)
'''

import argparse
import logging.config
import sys 

import requests.exceptions

from TaxonomyLineage import TaxonomyLineage
from api.BaseSeqdbApi import UnexpectedContent
from config import config_root
import tools_helper
from api.FeatureTypeApi import FeatureTypeApi
from api.FeatureApi import FeatureApi
from api.DeterminationApi import DeterminationApi


### Values below are used in Galaxy wrappers, so make sure you know what you're doing if you're changing any of them 
# Values for the types of sequences this script downloads. I.e. "its" loads 
# ITS sequences. Note that raw sequences are not implemented in SeqDB yet. "raw":"raw",
push_types_dict = {"its":"itsx_features", "taxonomy":"findLCA_taxonomy"}


###


def parse_input_args(argv):
    ''' Parses command line arguments '''
    
    push_types_set = frozenset(push_types_dict.values())
    
    parser = argparse.ArgumentParser(description="Write information to SeqDB")
    parser.add_argument('push_type', help="Type of information to write", type=str, choices=push_types_set)
    parser.add_argument('-c', help="SeqDB config file", dest="config_file", required=False)
    parser.add_argument('-u', help="SeqDB API URL", dest="base_url", required=False)
    parser.add_argument('-k', help="SeqDB API key", dest="api_key", required=False)    
    parser.add_argument('--lca_results_file', help="FindLCA results file", dest="lca_results_file", required=False)    
    parser.add_argument('--itsx_positions_file', help="ITSx feature positions file (.positions.txt)", dest="itsx_positions_file", required=False)    
    parser.add_argument('--itsx_extraction_file', help="ITSx extraction results file (.extraction.results)", dest="itsx_extraction_file", required=False)    
    
    args = parser.parse_args(argv)

    if not (args.config_file or (args.base_url and args.api_key)):
        parser.error('Either -c <configuration file>, or -u <base_url> -k <api_key> have to be specified')
    
    if args.push_type == push_types_dict["its"] and not (args.itsx_positions_file and args.itsx_extraction_file):
        parser.error('To write its features to SeqDB, ITSx positions (.positions.txt) \
        and extraction results (extraction.results) files need to be supplied. \
        Use --itsx_positions_file and --itsx_extraction_file options respectively.')
        
    if args.push_type == push_types_dict["taxonomy"] and not args.lca_results_file:
        parser.error('To write taxonomy information to SeqDB, FindLCA results file \
        needs to be supplied. Use --lca_results_file option.')
              
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
            error_msg = "Could not open file: {}.".format(file_name)
            logging.error(error_msg)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        else:
            raise

    return file_handler


def get_lieage_taxids(tax_parent_ids, taxon_id, lineage=None):
    ''' Finds a taxonomic lineage for taxon_id.
    Return: list of taxid in lineage order, from given tax id up. 
    '''
    

def push_taxonomy_data(determinationApi, info_file_name, taxonomy_dir):
    ''' Extracts taxon id from the findLCA output file, finds lineage for the taxon id
        and pushes taxonomic identification to SeqDB
    Args:
        info_file_name: FindLCA output file. Example of 3 lines in a file:
        Query: seqdb|6802    No suitable matches found.
        Query: seqdb|28954    LCA: 164328    Name:Phytophthora ramorum    Rank: species    Matches: 1    Names: Phytophthora ramorum:55,    Evalue: 0.0,    Pid: 96.97,    Qcover: 0.929705215419501,    SourceGI: 46812520,74476198,78192394,81176726,81176727,81176728,408690088,378750406,378750408,643192938,643192939,44194092,47934204,49617503,60308861,89275889,108885436,114146744,55586083,643192936,55586084,46812521,42521138,32490530,315419600,643192937,55586085,643192935,227810703,227810722,46399123,323301609,227810723,295409544,284432200,308745533,551031928,227810464,227810545,227810570,372125558,284432218,284432396,308745534,478739030,227810622,227810639,42529489,227810623,284432190,227810547,417357152,417357145,38348758,171262896,
        Query: seqdb|64071    LCA: 4783    Name:Phytophthora    Rank: genus    Matches: 4    Names: Phytophthora ramorum:2,Phytophthora mengei:1,Phytophthora citricola:1,Phytophthora capsici:2,    Evalue: 0.0,    Pid: 95.00,95.15,    Qcover: 0.928864569083447,    SourceGI: 320336471,308913225,308913137,308912963,308912935,320336197,
    
    Return:
        list of determination ids, written to seqdb
    '''
    
    # Output file name for this execution. Is used in Galaxy, so don't change UYKWYD
    output_file_name = 'seqdb_taxon_ids.txt'
    
    input_file_line_example = """ Query: seqdb|6726    LCA: 129355  Name:Phytophthora lateralis    Rank: species    Matches: 1    Names: Phytophthora lateralis:1,  Evalue: 0.0,    Pid: 98.62,    Qcover: 0.867647058823529,    SourceGI: 320336337, 
    or
    Query: seqdb|6802    No suitable matches found."""
    
    # Parse input file and extract (SeqDB sequence ID, NCBI Taxon ID) pairs
    seqdb_sequence_taxon = {}
    try:
        with open(info_file_name, "r") as info_file_handler:
            
            msg_inputFile = "Find LCA results file: {}".format(info_file_name)
            logging.info(msg_inputFile)
    
            for line in info_file_handler:
                error_msg = "{} Found line: \n{}\n Example of an expected line: \n{}".format(tools_helper.log_msg_wrongFileFormat, line, input_file_line_example)
                
                line_tokens = line.split('\t')
                
                # Do basic file format verification:
                if len(line_tokens) != 10 and len(line_tokens) != 2:
                    logging.error("{} Example of an expected line:\n{}".format(tools_helper.log_msg_wrongFileFormat, input_file_line_example))
                    sys.exit(tools_helper.log_msg_sysExitFile)
        
                try:
                    token1 = line_tokens[0].split(': ')
                    if token1[0] != "Query":
                        logging.error(error_msg)
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
                    logging.error(error_msg)
                    sys.exit(tools_helper.log_msg_sysExitFile)
        
    except IOError as e:
        if e.errno == 2:
            error_msg = "Could not open input file '{}'.".format(info_file_name)
            logging.error(error_msg)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        else:
            raise


    # For each (sequence_id, taxonomy_id) pair, find full lineage (from NCBI) and write to SeqDB
    tax_lineage = TaxonomyLineage(taxonomy_dir)
    
    determinationIds = list()
    for sequenceId in seqdb_sequence_taxon:
        taxonomyId = seqdb_sequence_taxon[sequenceId]
        notes = "Taxonomy derived from BLAST / LCA pipeline in Galaxy. "
        lineage_names = {}
        
        if isinstance( taxonomyId, ( int, long ) ):
            #start_time = time.time()
            lineage_names = tax_lineage.findLineage(taxonomyId)   
            #print "Time to get lineage new way: %s" %(time.time()-start_time)
        else:
            notes = taxonomyId + " " + notes
        
        try:
            determinationId = determinationApi.create_sequence_determination(sequence_id=sequenceId,
                                                                             taxonomy=lineage_names, is_accepted=False,
                                                                             notes=notes)
            determinationIds.append(determinationId)
        except requests.exceptions.ConnectionError as e:
            logging.error(tools_helper.log_msg_noDbConnection)
            logging.error(e.message)
            if repr(e.args[0].args[1]) == tools_helper.seqdbConnection_error_msg:
                logging.error(tools_helper.log_msg_seqDBConnection)
            if repr(e.args[0].args[1]) == tools_helper.invalidURL_error_msg:
                logging.error(tools_helper.log_msg_invalidURL)
            sys.exit(tools_helper.log_msg_sysExit)
        except requests.exceptions.HTTPError as e:
            logging.error(tools_helper.log_msg_httpError)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except UnexpectedContent as e:
            logging.error(tools_helper.log_msg_apiResponseFormat)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        except:
            warning_msg = "Unexpected error writing determination"
            logging.info(warning_msg)

    
    # Write determination IDs, which were written to SeqDB, to a file (for Galaxy purposes)
    output_file = open(output_file_name, 'w')
    for detId in determinationIds:
        output_file.write(str(detId) + '\n')
    output_file.close()
  
    log_msg1 = "Number of determinations written to Sequence Database: {}".format(len(determinationIds))
    log_msg2 = "Created determination IDs are written to a file: '{}'".format(output_file_name)
    logging.info(log_msg1)
    logging.info(log_msg2)
    
    logging.info("Determination IDs, written to SeqDB: {}".format(determinationIds))
   
    return determinationIds


# B15_17_SH817_ITS_ITS5    622 bp.    SSU: Not found    ITS1: 1-241    5.8S: 242-399    ITS2: 400-557    LSU: 558-622
def push_its_features(featureTypeApi, featureApi, features_file_name, extraction_results_file_name=None):
    ''' Extracts found ITS features from ITSx tools results .positions.txt file and writes these features to SeqDB     
    
    Return:
        list of SeqDB feature_ids for those features that were successfully written to SeqDB
    '''
    # Output file name for this execution. Is used in Galaxy, so don't change UYKWYD
    output_file_name = 'seqdb_feature_ids.txt'
   
    # The name of the feature type for each ITSx feature in SeqDB
    its_feature_type_name = "misc_RNA"  # agreed on convention
    rRna_feature_type_name = "rRNA"     # agreed on convention
    feature_description = "Generated by ITSx tool."
    input_file_line_example = "seqdb|1685    321 bp.    SSU: 1-46    ITS1: 47-283    5.8S: No end    ITS2: Not found    LSU: Not found    Broken or partial sequence, only partial 5.8S!"

    created_feature_ids = []
    
    # Get all feature types from seqDB
    seqdb_feat_types = featureTypeApi.getFeatureTypesWithIds()
    
    its_feature_type_id = ''
    rRna_feature_type_id = ''
    # Create feature type in SeqDB, if doesn't exist already
    if its_feature_type_name not in seqdb_feat_types:
        its_feature_type_id = featureTypeApi.create(its_feature_type_name)
    else:
        its_feature_type_id = seqdb_feat_types[its_feature_type_name]
        
    if rRna_feature_type_name not in seqdb_feat_types:
        rRna_feature_type_id = featureTypeApi.create(rRna_feature_type_name)
    else:
        rRna_feature_type_id = seqdb_feat_types[rRna_feature_type_name]
    
            
    ## Read strand information from the .extraction.results file
    match_strand = {}
    if extraction_results_file_name:
        msg_inputFile = "ITSx extraction results file (.extraction.results): {}".format(extraction_results_file_name)
        logging.info(msg_inputFile)
    
        try:
            with open(extraction_results_file_name, "r") as detailed_results_file:
                for line in detailed_results_file:
                    line = line.strip()
                    if line:
                        line_tokens = line.split('\t')
                        sequenceId = int(line_tokens[0].split('|')[1])
                        strand = int(line_tokens[3])
                        match_strand[sequenceId] = strand
        except IOError as e:
            if e.errno == 2:
                error_msg = "Could not open input file '{}'.".format(extraction_results_file_name)
                logging.error(error_msg)
                logging.error(e.message)
                sys.exit(tools_helper.log_msg_sysExit)
            else:
                raise
    
    ## Parse the .positions.txt file and insert each found feature to SeqDB
    
    msg_inputFile = "ITSx feature positions file (.positions.txt): {}".format(features_file_name)
    logging.info(msg_inputFile)
    
    
    info_file_handler = ''
    try:
        info_file_handler = open(features_file_name,"r")
    except IOError as e:
        if e.errno == 2:
            error_msg = "Could not open input file '{}'.".format(features_file_name)
            logging.error(error_msg)
            logging.error(e.message)
            sys.exit(tools_helper.log_msg_sysExit)
        else:
            raise
        
    for line in info_file_handler:
        line_tokens = line.split('\t')

        try:
            sequenceId = line_tokens[0].split('|')[1]
            sequenceId = int(sequenceId)
        except:
            error_msg = "Input file error. Could not extract sequence id from the input file. Example of an expected line:\n{}\n".format(input_file_line_example)
            logging.error(error_msg)
            sys.exit(tools_helper.log_msg_sysExit)
        
        try:
            itsx_features = line_tokens[2:7]
        except IndexError as e:
            error_msg = "Could not extract ITS features from the input file. Example of an expected line:\n{}\n".format(input_file_line_example)
            logging.error(error_msg)
            sys.exit(tools_helper.log_msg_sysExit)
        
        if not sequenceId or not itsx_features:
            error_msg = "Input file not in the correct format. Example of an expected line:\n{}\n".format(input_file_line_example)
            logging.error(error_msg)
            sys.exit(tools_helper.log_msg_sysExit)
            
        ## Strand information is read from .extraction.results file (above)
        ## 0 is main strand, 1 is complementary (same convention as in ITSx tool)
        if match_strand and sequenceId in match_strand:
            strand = match_strand[sequenceId]
        else:
            strand = 0  
            
                
        for itsx_feature_token in itsx_features:
            feature_location_pair = itsx_feature_token.split(": ")
            try:
                location = feature_location_pair[1].split("-")
                location = map(int, location)
                
                ### TODO: Finalize what to use for the frame: 
                
                location = [{"start":location[0],"end":location[1],"frame":1,"strand":strand}]
                feature_type_id = its_feature_type_id
                if feature_location_pair[0] in {"SSU", "5.8S", "LSU"}:
                    feature_type_id = rRna_feature_type_id  
                fid = featureApi.create(feature_location_pair[0], feature_type_id, location, sequenceId, description=feature_description)
                created_feature_ids.append(fid)
            except requests.exceptions.ConnectionError as e:
                logging.error(tools_helper.log_msg_noDbConnection)
                logging.error(e.message)
                if repr(e.args[0].args[1]) == tools_helper.seqdbConnection_error_msg:
                    logging.error(tools_helper.log_msg_seqDBConnection)
                if repr(e.args[0].args[1]) == tools_helper.invalidURL_error_msg:
                    logging.error(tools_helper.log_msg_invalidURL)
                sys.exit(tools_helper.log_msg_sysExit)
            except requests.exceptions.HTTPError as e:
                logging.error(tools_helper.log_msg_httpError)
                logging.error(e.message)
                sys.exit(tools_helper.log_msg_sysExit)
            except UnexpectedContent as e:
                logging.error(tools_helper.log_msg_apiResponseFormat)
                logging.error(e.message)
                sys.exit(tools_helper.log_msg_sysExit)
            except:
                warning_msg = "File token '{}' is not in the expected format of <feature name>:<position>. Ignoring.".format(itsx_feature_token)
                logging.info(warning_msg)

    info_file_handler.close()

    # Write IDs of the inserted features into a file
    output_file = open(output_file_name, 'w')
    for fid in created_feature_ids:
        output_file.write(str(fid) + '\n')
    output_file.close()
  
    log_msg1 = "Number of features written to Sequence Database: {}".format(len(created_feature_ids))
    log_msg2 = "Created feature IDs are written to a file: '{}'".format(output_file_name)
    logging.info(log_msg1)
    logging.info(log_msg2)
        
    return created_feature_ids
        
    
def main():
    ''' Write provided information to SeqDB '''
    
    
    ### Load main configuration file and set up logging for the script
    main_conf = tools_helper.load_config(config_root.path() + '/config.yaml')

    if not main_conf:
        logging.error(tools_helper.log_msg_noConfig)
        sys.exit(tools_helper.log_msg_sysExit)
    
    logging.config.dictConfig(main_conf['logging'])
    
    logging.info("{} '{}'".format(tools_helper.log_msg_scriptExecutionWithParams, sys.argv))
    
    
    ### Parse sript's input arguments 
    parsed_args = parse_input_args(sys.argv[1:])
    
    if parsed_args.config_file:
        tool_config = tools_helper.load_config(parsed_args.config_file)
        base_url = tool_config['seqdb']['base_url'] 
        api_key = tool_config['seqdb']['api_key'] 
    else:
        base_url = parsed_args.base_url 
        api_key = parsed_args.api_key 
    
    logging.info("{} '{}'".format(tools_helper.log_msg_apiUrl, base_url))
    
    
    ### Script execution             
    if push_types_dict["its"] == parsed_args.push_type:
        #sys.exit("Not yet implemented")
        log_msg = "Writing ITS features to SeqDB."
        logging.info(log_msg)
        featureTypeApi = FeatureTypeApi(api_key=api_key, base_url=base_url)
        featureApi = FeatureApi(api_key=api_key, base_url=base_url)
        success_feat_ids = push_its_features(featureTypeApi, featureApi, parsed_args.itsx_positions_file, parsed_args.itsx_extraction_file)
        print success_feat_ids
            
        
    elif push_types_dict["taxonomy"] == parsed_args.push_type:
        log_msg = "Writing taxonomy lineage information to SeqDB."
        logging.info(log_msg) 
        determinationApi = DeterminationApi(api_key=api_key, base_url=base_url)
        push_taxonomy_data(determinationApi, parsed_args.lca_results_file, main_conf['galaxy']['ncbi_taxonomy_dir'])

    
    ### Post-execution: messages and logging
    print(tools_helper.log_msg_execEnded)
    logging.info(tools_helper.log_msg_execEnded)
    
    
if __name__ == '__main__':
# Usage exmaple for ITS features:
#-c user_config.yaml itsx_features --itsx_positions_file ./test/data/test.positions.txt --itsx_extraction_file ./test/data/test.extraction.results

# Usage exmaple for LCA:
#-c user_config.yaml findLCA_taxonomy --lca_results_file ./test/data/LCA_Results.tabular
    main()
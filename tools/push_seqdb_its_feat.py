'''
Created on Mar 5, 2015

@author: korolo
'''
import sys, getopt, logging
from api.seqdbWebService import seqdbWebService, UnexpectedContent
#from unittest.test.support import LoggingResult


usage_help_line = """Usage of the script: \npush_seqdb_its_feat.py -k <SeqDB API key> -f <features file name>
Other arguments:
   -h       help (prints this message)
   -p       use production url for web service requests
   -u <url> specify base url for web services requests"""


prod_url = "***REMOVED***/api/v1"
local_url = "http://localhost:8181/seqdb.web-2.5/api/v1"

# The name of the feature type for each itsx feature
its_feature_type_name = "misc_RNA"
rRna_feature_type_name = "rRNA"
feature_description = "Generated by ITSx tool."
output_file_name = 'seqdb_feature_ids.txt'
log_file_name = "seqdb_push.log"
input_file_line_example = "seqdb|1685    321 bp.    SSU: 1-46    ITS1: 47-283    5.8S: No end    ITS2: Not found    LSU: Not found    Broken or partial sequence, only partial 5.8S!"


# Parses command line arguments 
# Returns:
#    seqdb api_key
#    name of the ITSx poisitons file (containing feature to import)
#    base url to use for web services requests
def parse_input_args(argv):
    seqdb_api_key = ''
    base_url=''
    features_file_name = ''
    
    prod=False
    user_url=''
    
    try:
        opts, args = getopt.getopt(argv,"hpk:f:u:",["seqdb_api_key="])
    except getopt.GetoptError:
        print usage_help_line
        sys.exit(2)
        
    if len(opts)==0:
        print usage_help_line
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print usage_help_line
            sys.exit()
        elif opt in ("-k", "--seqdb_api_key"):
            seqdb_api_key = arg
        elif opt in ("-f", "--features_file"):
            features_file_name = arg
        elif opt == "-p":
            prod=True
        elif opt == "-u":
            user_url = arg
            
    if user_url:
        base_url = user_url
    elif prod:
        base_url = prod_url
    else:
        base_url = local_url
    
    return seqdb_api_key, features_file_name, base_url




# B15_17_SH817_ITS_ITS5    622 bp.    SSU: Not found    ITS1: 1-241    5.8S: 242-399    ITS2: 400-557    LSU: 558-622
def main(api_key, features_file_name, base_url):
    print_error_line = "Execution was not successful. See log for details: '%s'" %log_file_name
    
    created_feature_ids = []
    
    # Open an ITSx positions file
    feat_file = ''
    try:
        feat_file = open(features_file_name,"r")
    except IOError as e:
        if e.errno == 2:
            logging.error("Could not open ITSx feature positions file '%s'." % features_file_name)
            logging.error(e.message)
            print "Could not open ITSx feature positions file. See log file for details."
            sys.exit(1)
        else:
            raise
    
    seqdbWS = seqdbWebService(api_key, base_url)
    
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
    for line in feat_file:
        line_tokens = line.split('\t')

        try:
            sequenceId = line_tokens[0].split('|')[1]
            sequenceId = int(sequenceId)
        except:
            logging.error("Could not extract sequence id from the input file. Example of an expected line:\n%s\n" % input_file_line_example)
            print print_error_line
            sys.exit(1)
        
        try:
            itsx_features = line_tokens[2:7]
        except IndexError as e:
            logging.error("Could not extract ITS features from the input file. Example of an expected line:\n%s\n" % input_file_line_example)
            logging.error(e.message)
            print print_error_line
            sys.exit(1)
        
        if not sequenceId or not itsx_features:
            logging.error("Input file not in the correct format. Example of an expected line:\n%s\n" % input_file_line_example)
            print print_error_line
            sys.exit(1)
            
                
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
            except UnexpectedContent as e:
                print e
                sys.exit(1)
            except:
                logging.warning("File token '%s' is not in the expected format of <feature name>:<position>. Ignoring." % itsx_feature_token)


    # Write ids of the inserted features into a file
    output_file = open(output_file_name, 'w')
    for fid in created_feature_ids:
        output_file.write(str(fid) + '\n')
    output_file.close()
  
    logging.info("Number of features written to Sequence Dababase:   %i " % len(created_feature_ids))  
    logging.info("Created feature ids are written to a file: '%s'" % output_file_name)
        
    return created_feature_ids

    

if __name__ == '__main__':
    # Start a log file. filemode='w' overwrites the log for each program run
    logging.basicConfig(filename=log_file_name, filemode='w', level=logging.DEBUG)
    
    logging.info("Script executed with the following command and arguments: %s" % sys.argv)
    
    seqdb_api_key, features_file_name, base_url = parse_input_args(sys.argv[1:])
    
    logging.info("Base URL for web services is: '%s'" % base_url)
    logging.info("File name with features: %s" % features_file_name)
    
    ok_feat_ids = main(seqdb_api_key, features_file_name, base_url)
    
    print "Execution complete."
    print "Number of features written to Sequence Dababase:   %i " % len(ok_feat_ids)  
    print "Created feature ids are written to a file: '%s'" % output_file_name
    print "Execution log is written to a file: '%s'" % log_file_name


'''
Created on Mar 5, 2015

Deletes features

@author: korolo
'''
import sys, getopt, logging
from api.seqdbWebService import seqdbWebService, UnexpectedContent


usage_help_line = """Usage of the script: \ndelete_seqdb_feature.py -k <SeqDB API key> -f <feature ids file name>
Other arguments:
   -h       help (prints this message)
   -p       use production url for web service requests
   -u <url> specify base url for web services requests"""


prod_url = "***REMOVED***/api/v1"
local_url = "http://localhost:8181/seqdb.web-2.5/api/v1"
output_file_name = "delete_failed_feature_ids.txt"
log_file_name = "seqdb_delete.log"



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
        opts, args = getopt.getopt(argv,"hpk:f:u:",["seqdb_api_key=", "feature_ids_file="])
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
        elif opt in ("-f", "--feature_ids_file"):
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




def delete_features(api_key, feature_ids_file_name, base_url):
    ''' Deletes features with the ids, specified in the file. 
    File should contain feature ids either:
        1) on a separate line; 
        2) one line, separated by comma; 
        3) separated by tab
    '''
    #print "Sending request to: " + base_url
    
    # Open an ITSx positions file
    feat_file = ''
    try:
        feat_file = open(feature_ids_file_name,"r")
    except IOError as e:
        if e.errno == 2:
            logging.error("Could not open feature ids file <%s>." % features_file_name)
            logging.error(e.message)
            print "Could not open feature ids file. See log file for details."
            sys.exit(1)
        else:
            raise
    
    seqdbWS = seqdbWebService(api_key, base_url)
    
    feature_ids = []
    line_number = 0
    # parse the file delete features
    for line in feat_file:
        line_number = line_number +1
        current_feature_ids = []
        
        if line.__contains__(','):
            current_feature_ids = line.replace(' ', '').split(',')
        elif line.__contains__('\t'):
            current_feature_ids = line.split('\t')
        else:
            current_feature_ids = [line]
        
        try:
            current_feature_ids = map(int, current_feature_ids)
            feature_ids.extend(current_feature_ids)
        except:
            logging.warning("Line number: %i. Could not parse line '%s' for feature ids. Ignoring." % (line_number,line.replace('\n', '')))
    
    logging.info("Identified %i features to be deleted." % len(feature_ids))
    
    success_ids = []
    fail_ids = []        
    for feature_id in feature_ids:  
        try:
            seqdbWS.deleteFeature(feature_id)
            success_ids.append(feature_id)
        except UnexpectedContent as e:
            print e
            sys.exit(1)
        except:
            fail_ids.append(feature_id)    
    
        
    output_file = open(output_file_name, 'w')
    for fid in fail_ids:
        output_file.write(str(fid) + '\n')        
    output_file.close()        
    
    logging.info("Number of features deleted from Sequence Dababase:   %i " % len(success_ids))  
    logging.info("Number of features which failed to be deleted:   %i " % len(fail_ids))
    logging.info("Feature ids, which failed to be deleted are written to a file: '%s'" % output_file_name)
   
    return success_ids, fail_ids
    
    
def main():
    ''' Deletes features from SeqDB with the feature ids, specified in the input file. '''
    
    # Start a log file. filemode='w' overwrites the log for each program run
    logging.basicConfig(filename=log_file_name, filemode='w', level=logging.DEBUG)
    
    logging.info("Script executed with the following command and arguments: %s" % sys.argv)
    
    seqdb_api_key, features_file_name, base_url = parse_input_args(sys.argv[1:])
    
    logging.info("Base URL for web services is: '%s'" % base_url)
    logging.info("File name with feature ids to delete: %s" % features_file_name)
    
    success_ids,fail_ids = delete_features(seqdb_api_key, features_file_name, base_url)
    
    
    print "Execution complete."
    print "Number of features deleted from Sequence Dababase:   %i " % len(success_ids)  
    print "Number of features which failed to be deleted:   %i " % len(fail_ids)  
    print "Feature ids, which failed to be deleted are written to a file: '%s'" % output_file_name
    print "Execution log is written to a file: '%s'" % log_file_name


if __name__ == '__main__':
    main()
    
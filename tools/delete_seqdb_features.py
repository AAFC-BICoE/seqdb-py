'''
Created on Mar 5, 2015

Deletes features

@author: korolo
'''
import sys, getopt
from api.seqdbWebService import seqdbWebService


usage_help_line = """Usage of the script: \ndelete_seqdb_feature.py -k <SeqDB API key> -f <feature ids file name>
Other arguments:
   -h       help (prints this message)
   -p       use production url for web service requests
   -u <url> specify base url for web services requests"""


prod_url = "***REMOVED***/api/v1"
local_url = "***REMOVED***:8181/seqdb.web-2.5/api/v1"



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




# File should contain feature ids either: 1) on a separate line; 2) one line, separated by comma; 3) separated by tab
def main(api_key, feature_ids_file_name, base_url):
    #print "Sending request to: " + base_url
    
    # Open an ITSx positions file
    feat_file = ''
    try:
        feat_file = open(feature_ids_file_name,"r")
    except IOError as e:
        if e.errno == 2:
            print "Could not open feature ids file <%s>." % feature_ids_file_name
            sys.exit(1)
        else:
            raise
    
    seqdbWS = seqdbWebService(api_key, base_url)
    
    feature_ids = []
    
    # parse the file delete features
    for line in feat_file:
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
            pass
    
    success_ids = []
    fail_ids = []        
    for feature_id in feature_ids:  
        try:
            seqdbWS.deleteFeature(feature_id)
            success_ids.append(feature_id)
        except:
            fail_ids.append(feature_id)    
    
    
    return success_ids, fail_ids
    

if __name__ == '__main__':
    
    sys.path.append("/Users/korolo/git/seqdb-py")
    sys.path.append("/Users/korolo/git/seqdb-py/api")
    
    seqdb_api_key, features_file_name, base_url = parse_input_args(sys.argv[1:])
    main(seqdb_api_key, features_file_name, base_url)
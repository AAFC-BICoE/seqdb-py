'''
Created on Mar 5, 2015

Deletes features

@author: korolo
'''
import sys
import getopt
import logging.config
import tools_helper
from api.seqdbWebService import seqdbWebService, UnexpectedContent
from config import config_root

usage_help_line = """Usage of the script: \ndelete_seqdb_feature -c <Path to yaml config with SeqDB API info> -f <feature ids file name>
or
delete_seqdb_feature -k <SeqDB API key> -u <SeqDB API URL> -f <feature ids file name>
Other arguments:
   -h   help (prints this message)
"""


output_file_name = "delete_failed_feature_ids.txt"
log_file_name = "seqdb_delete.log"
user_log = tools_helper.SimpleLog(log_file_name)
log_fail_msg = "%s '%s'" %(tools_helper.log_msg_errorSeeLog, log_file_name)



def parse_input_args(argv):
    ''' Parses command line arguments
    Returns:
        features_file_name: name of the file that contains ITS features to be pushed to SeqDB
        config_file: path to a config file with has api information, 
            or empty string if no such usage 
        seqdb api_key to use for web services requests
        seqdb api_url to use for web services requests
    '''
    config_file=''
    api_url=''
    api_key = ''
    features_file_name = ''
    
    
    try:
        opts, args = getopt.getopt(argv,"hc:f:k:u:",["config_file=", "seqdb_api_key=", "seqdb_api_url=", "feature_ids_file="])
    except getopt.GetoptError:
        print usage_help_line
        logging.error(tools_helper.log_msg_argError)
        user_log.error(tools_helper.log_msg_argError)
        sys.exit(tools_helper.log_msg_sysExit)
        
        
    if len(opts)==0:
        print usage_help_line
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print usage_help_line
            sys.exit()
        elif opt in ("-f", "--features_file"):
            features_file_name = arg
        elif opt in ("-c", "--config_file="):
            config_file = arg
        elif opt in ("-k", "--seqdb_api_key"):
            api_key = arg
        elif opt in ("-u", "--seqdb_api_url="):
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
    
            
            
    return (features_file_name, config_file, api_url, api_key)



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
            logging.error("Could not open feature ids file <%s>." % feature_ids_file_name)
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
    
    main_conf = tools_helper.load_config(config_root.path() + '/config.yaml')

    if not main_conf:
        logging.error(tools_helper.log_msg_noConfig)
        sys.exit(tools_helper.log_msg_sysExit)
    
    logging.config.dictConfig(main_conf['logging'])

    logging.info("Script executed with the following command and arguments: %s" % sys.argv)
    user_log.info(tools_helper.log_msg_execStarted_simple)
    
    features_file_name, config_file, api_url, api_key = parse_input_args(sys.argv[1:])
    
    if config_file:
        tool_config = tools_helper.load_config(config_file)
        api_url = tool_config['seqdb']['api_url'] 
        api_key = tool_config['seqdb']['api_key'] 
    
    
    logging.info("%s '%s'" %  (tools_helper.log_msg_apiUrl, api_url))
    user_log.info("%s '%s'" %  (tools_helper.log_msg_apiUrl, api_url))

    log_msg = "File name with feature ids: %s" % features_file_name
    logging.info(log_msg)    
    user_log.info(log_msg)
    
    
    success_ids,fail_ids = delete_features(api_key, features_file_name, api_url)

    
    print(tools_helper.log_msg_execEnded)
    print("Number of features deleted from Sequence Dababase:   %i " % len(success_ids))  
    print("Number of features which failed to be deleted:   %i " % len(fail_ids))
    print("Feature ids, which failed to be deleted are written to a file: '%s'" % output_file_name)
    print("Execution log is written to a file: '%s'" % log_file_name)
    
    user_log.info(tools_helper.log_msg_execEnded)
    user_log.close()
    
    logging.info(tools_helper.log_msg_execEnded)

    
    

if __name__ == '__main__':
    main()
    
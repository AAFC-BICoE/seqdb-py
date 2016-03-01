'''
Created on Mar 5, 2015

Deletes features

@author: korolo
'''
import argparse
import logging.config
import sys

from api.seqdbWebService import seqdbWebService, UnexpectedContent
from config import config_root
import tools_helper


output_file_name = "delete_failed_feature_ids.txt"
delete_types_dict = {"feature":"features", "determination":"taxonomy"}


def parse_input_args(argv):
    ''' Parses command line arguments '''
    delete_types_set = frozenset(delete_types_dict.values())
    
    parser = argparse.ArgumentParser(description="Delete items from SeqDB")
    parser.add_argument('delete_type', help="Type of information to be deleted", type=str, choices=delete_types_set)
    parser.add_argument('-c', help="SeqDB config file", dest="config_file", required=False)
    parser.add_argument('-u', help="SeqDB API URL", dest="api_url", required=False)
    parser.add_argument('-k', help="SeqDB API key", dest="api_key", required=False)    
    parser.add_argument('-f', help="File with feature ids to be deleted from SeqDB", dest="features_file_name", required=True)    
    
    args = parser.parse_args(argv)

    if not (args.config_file or (args.api_url and args.api_key)):
        parser.error('Either -c <configuration file>, or -u <api_url> -k <api_key> have to be specified')
        
    return args


def delete_from_seqdb(seqdbWS, seqdb_ids_file_name, delete_type):
    ''' Deletes items from seqDB. Need to specify SeqDB API connection details (seqdbWS),
        type of items to delete (delete_type), and ids of items to be deleted (seqdb_ids_file_name)
    Agrs:
        seqdbWS: api.seqdbWebService object with accessor methods to SeqDB
        delete_type: type of items to delete from SeqDB. Specified in the delete_types_dict
        seqdb_ids_file_name: name of the file with the seqdb ids to be deleted
            File should contain ids either:
                1) on a separate line; 
                2) one line, separated by comma; 
                3) separated by tab
    Return:
        success_ids: list of ids that where successfully deleted from SeqDB
        fail_ids: list of ids that failed to be deleted from SeqDB
    '''
    
    ids_file = ''
    try:
        ids_file = open(seqdb_ids_file_name,"r")
    except IOError as e:
        if e.errno == 2:
            logging.error("Could not open file <%s>." % seqdb_ids_file_name)
            logging.error(e.message)
            print "Could not open seqdb ids file. See log file for details."
            sys.exit(1)
        else:
            raise
    
    delete_ids = []
    line_number = 0
    
    for line in ids_file:
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
            delete_ids.extend(current_feature_ids)
        except:
            logging.warning("Line number: %i. Could not parse line '%s'. Ignoring." % (line_number,line.replace('\n', '')))
    
    logging.info("Identified %i items to be deleted." % len(delete_ids))
    
    success_ids = []
    fail_ids = []        
    for delete_item_id in delete_ids:  
        try:
            if delete_types_dict["feature"] == delete_type:
                seqdbWS.deleteFeature(delete_item_id)
            elif delete_types_dict["determination"] == delete_type:
                seqdbWS.deleteDetermination(delete_item_id)
            success_ids.append(delete_item_id)
        except UnexpectedContent as e:
            print e
            sys.exit(1)
        except:
            fail_ids.append(delete_item_id)    
    
        
    output_file = open(output_file_name, 'w')
    for fid in fail_ids:
        output_file.write(str(fid) + '\n')        
    output_file.close()        
    
    return success_ids, fail_ids
    
    
def main():
    ''' Deletes items from SeqDB. '''
    
    main_conf = tools_helper.load_config(config_root.path() + '/config.yaml')

    if not main_conf:
        logging.error(tools_helper.log_msg_noConfig)
        sys.exit(tools_helper.log_msg_sysExit)
    
    logging.config.dictConfig(main_conf['logging'])

    logging.info("Script executed with the following command and arguments: %s" % sys.argv)
    
    parsed_args = parse_input_args(sys.argv[1:])
    
    if parsed_args.config_file:
        tool_config = tools_helper.load_config(parsed_args.config_file)
        api_url = tool_config['seqdb']['api_url'] 
        api_key = tool_config['seqdb']['api_key'] 
    else:
        api_url = parsed_args.api_url 
        api_key = parsed_args.api_key 
    
    logging.info("%s '%s'" %  (tools_helper.log_msg_apiUrl, api_url))
    
    seqdbWS = seqdbWebService(api_key, api_url)
             
    log_msg = "Deleting items from SeqDB. File name with SeqDB ids to be deleted: %s" % parsed_args.features_file_name
    logging.info(log_msg)    
        
        
    success_ids,fail_ids = delete_from_seqdb(seqdbWS=seqdbWS, 
                                             delete_type=parsed_args.delete_type,
                                             seqdb_ids_file_name=parsed_args.features_file_name) 
    
        
    
    ### Post-execution: messages and logging
    logging.info("Number of items deleted from Sequence Dababase:   %i " % len(success_ids))  
    print("Number of items deleted from Sequence Dababase:   %i " % len(success_ids))  
    logging.info("Number of items which failed to be deleted:   %i " % len(fail_ids))
    print("Number of items which failed to be deleted:   %i " % len(fail_ids))
    logging.info("IDs, which failed to be deleted are written to a file: '%s'" % output_file_name)
   
    print(tools_helper.log_msg_execEnded)
    
    
    logging.info(tools_helper.log_msg_execEnded)

    
    

if __name__ == '__main__':
    main()
    
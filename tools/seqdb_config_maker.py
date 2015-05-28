'''
Created on May 28, 2015

@author: korolo

Usage:
python pull_seqdb_its_seqs.py -k <SeqDB API key>
or
python pull_seqdb_its_seqs.py --seqdb_api_key <SeqDB API key>
-h is help
-p is to use production URL (uses ***REMOVED*** (Oksana's) machine by default, since we don't have UAT for SeqDB as of this time)
'''
import sys
import os
import getopt
import logging
import yaml
from api.seqdbWebService import seqdbWebService, UnexpectedContent


usage_help_line = """Usage of the script: \nseqdb_config_maker -k <SeqDB API key>
Other arguments:
   -h   help (prints this message)
 """

prod_url = "***REMOVED***/api/v1"
local_url = "***REMOVED***:2002/seqdb/api/v1"
# file name where the received sequences will be stores
output_file_name = "seqdb_ITS_seqs.fasta"
log_file_name = "seqdb_pull.log"


class seqdbConfigMaker:
    """ Creates a user-specific yaml config file to be used in SeqDB Web Services calls. 
        Config file name, seqdb api url and other setup parameters are retrieved from config.yaml
    """
    
    def __init__(self, api_url='', config_file_name=''):

        if api_url and config_file_name:
            self.apiUrl = api_url
            self.configFileName = config_file_name
        else:
            main_config = self.loadConfig('../config.yaml')
            self.apiUrl = main_config['seqdb']['url']
            self.configFileName = main_config['user_config_file']
            
            
        
    def loadConfig(self, config_file):
        ''' Return a config object populated from "config.yaml" 

        Raises: 
            yaml.YAMLError
        '''
        try:
            config = yaml.load(file(config_file, 'r'))
            return config
        except yaml.YAMLError, exc:
            print("Error in configuration file:", exc)


    def createConfigFile(self, api_key):
        ''' Writes a yaml file with seqdb api configuration info
        Kwargs:
            api_key: seqdb api key
        Returns:
            absolute path to the config file that was created
        '''
        config_file = open(self.configFileName, 'w')
        config_file.write("seqdb:\n")
        config_file.write('    api_key: "%s"\n' % api_key)
        config_file.write('    api_url:"%s"\n' % self.apiUrl)
        config_file.close()
        
        return os.path.abspath(self.configFileName)

''' Class End '''


def parse_input_args(argv):
    ''' Parses command line arguments
    Returns:
        seqdb api_key to use for web services requests
    '''
    seqdb_api_key = ''
    
    try:
        opts, args = getopt.getopt(argv,"hk:",["seqdb_api_key=", ])
    except getopt.GetoptError:
        print(usage_help_line)
        sys.exit(2)
        
    if len(opts)==0:
        print(usage_help_line)
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print(usage_help_line)
            sys.exit()
        elif opt in ("-k", "--seqdb_api_key"):
            seqdb_api_key = arg
            
    
    return seqdb_api_key
    

def main():
    # Start a log file. filemode='w' overwrites the log for each program run
    logging.basicConfig(filename=log_file_name, filemode='w', level=logging.DEBUG)
    
    logging.info("Script executed with the following command and arguments: %s" % sys.argv)
    
    # Parse command line arguments
    api_key = parse_input_args(sys.argv[1:])
    
    #logging.info("Base URL for web services is: '%s'" % base_url)
   
    config_file = seqdbConfigMaker().createConfigFile(api_key)
    
    print("Configuration is written to a file: '%s'" % config_file)
    #print "Execution log is written to a file: '%s'" % log_file_name



if __name__ == '__main__':
    main()
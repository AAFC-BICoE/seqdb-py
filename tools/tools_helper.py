'''
Created on May 28, 2015

Collection of methods, used in tools

@author: korolo
'''
import yaml
import os
import datetime
import json
import logging


log_msg_sysAdmin = "Contact your sysadmin for more details."
log_msg_sysExit = "Script error."
log_msg_sysExitFile = "Input file error."
log_msg_noDbConnection = "Could not connect to Sequence DB."
log_msg_slowConnection = "Connection too slow for getting data from Sequence DB."
log_msg_httpError = "HTTP error when connecting to Sequence DB."
log_msg_apiResponseFormat = "Sequence DB API response did not match expected format."
log_msg_scriptError = "Script encountered an error."
log_msg_argError = "Invalid script arguments. Aborting."
log_msg_argErrorKeyUrl = "Both API url and API key have to be supplied."
log_msg_argErrorConfigFileUrl = "Supply either a configuration file or API url with key."
log_msg_noConfig = "Could not load configuration file. Exiting..."
log_msg_scriptExecutionWithParams = "Script executed with the following command and arguments:"
log_msg_execStarted_simple = "Script execution started."
log_msg_execEnded = "Execution complete."
log_msg_apiUrl = "Base URL for web services is:"
log_msg_errorSeeLog = "Execution was not successful. See log for details."

log_msg_ITSLoad="Loading ITS sequences."
log_msg_ConsensusLoad="Loading consensus sequences."
log_msg_RawLoad="Loading raw sequences for specimen: "

log_msg_wrongFileFormat="Supplied file is not in the required format. "


def load_config(config_file):
    ''' Return a config object populated from "config.yaml" 

    Raises: 
        yaml.YAMLError
    '''
    try:
        config = yaml.load(file(config_file, 'r'))
        return config
    except yaml.YAMLError, exc:
        print("Error in configuration file:", exc)
    except Exception, msg:
        print("Could not load configuration file: %s" % os.path.abspath(config_file))
        print("Current working script is: %s" % os.path.abspath(__file__))



def pretty_log_json(obj, level="info", message=None):
    """Pretty print an object as a JSON formatted str to the log

    Args:
        obj: Object to pretty log

    Kargs:
        level (str): One of "debug", "info", "warn", "error".  Default="info".

    Returns:
        None

    Raises:
        None

    Example usage:
    >>> pretty_log_json({"a":1})
    >>> pretty_log_json({"a":1}, level="debug")
    >>> pretty_log_json({"a":1}, message="Contents of obj:")
    >>> pretty_log_json({"a":1}, level="debug", message="Contents of obj:")
    """
    display = ""
    if message is not None:
        display = display + message + "\n"
    display = display + json.dumps(obj, sort_keys=True, indent=4)
    getattr(logging, level)(display)



class SimpleLog():
    def __init__(self, file_name):
        self.log_file = open(file_name, 'w')
        self.log_file_name = file_name
        
    
    def writeLogEntry(self, log_prefix, msg):
        self.log_file.write("%s - %s: %s \n" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), log_prefix, msg))
        
    def info(self, msg):
        self.writeLogEntry("INFO", msg)
        
    def warn(self, msg):
        self.writeLogEntry("WARN", msg)
        
    def error(self, msg):
        self.writeLogEntry("ERROR", msg)
        
    def getFullPath(self):
        return os.path.abspath(self.log_file_name)

    def getFileName(self):
        return os.path.basename(self.log_file_name)
    
    def close(self):
        self.log_file.close()
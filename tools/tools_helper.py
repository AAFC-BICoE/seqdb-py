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
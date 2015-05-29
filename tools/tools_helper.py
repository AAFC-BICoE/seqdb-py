'''
Created on May 28, 2015

Collection of methods, used in tools

@author: korolo
'''
import yaml
import os
import datetime

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
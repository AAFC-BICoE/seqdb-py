'''
This is a placeholder to be able to access the location of project root directory.
This is needed to be able to read config.yaml from different packages.
Created on Jun 1, 2015

@author: korolo
'''

import os

def path():
    return os.path.dirname(__file__)


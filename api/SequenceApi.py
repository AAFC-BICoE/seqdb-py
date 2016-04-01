'''
Created on Apr 1, 2016

@author: Oksana Korol
'''
from api.BaseSequenceEntity import BaseSequenceEntity

class SequenceApi(BaseSequenceEntity):
    '''
    classdocs
    '''


    def __init__(self, api_key, base_url):
        super(SequenceApi, self).__init__(api_key=api_key, base_url=base_url, request_url="sequence")
        
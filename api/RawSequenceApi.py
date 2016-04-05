'''
Created on Apr 1, 2016

@author: Oksana Korol
'''
from api.BaseSequenceEntity import BaseSequenceEntity

class RawSequenceApi(BaseSequenceEntity):
    '''
    classdocs
    '''


    def __init__(self, api_key, base_url):
        super(RawSequenceApi, self).__init__(api_key=api_key, base_url=base_url, request_url="sequence")
        
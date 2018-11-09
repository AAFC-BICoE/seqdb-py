'''
Created on April 06, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
'''

from api.BaseApiEntity import BaseApiEntity

class SeqSourceApi(BaseApiEntity):

    def __init__(self, api_key, base_url, sequence_id):
        # self.clearAllFilters()
        super(SeqSourceApi, self).__init__(api_key=api_key, base_url=base_url, request_url="sequence/{}/seqSource".format(sequence_id))
        
        
    def getParamsStr(self):
        return ''
'''
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
'''

from api.BaseApiEntity import BaseApiEntity
from api.BaseSeqdbApi import BaseSeqdbApi
from api.BaseSeqdbApi import UnexpectedContent


class SpecimenApi(BaseApiEntity):

    def __init__(self, api_key, base_url):
        super(SpecimenApi, self).__init__(api_key=api_key, base_url=base_url, request_url="specimen")
        
        
    @property
    def otherIdsFilter(self):
        return self.otherIdsFilter
    
    @otherIdsFilter.setter
    def otherIdsFilter(self, otherIdsFilter):
        ''' otherIdsFilter = code + identifier
        '''
        self.otherIdsFilter = otherIdsFilter
        
    def getParamsStr(self):
        params = ''
        if self.otherIdsFilter:
            params = "filterName=otherIds&filterValue={}&filterWildcard=true&".format(self.otherIdsFilter)
        
        return params
    
    def clearAllFilters(self):
        self.otherIdsFilter = None
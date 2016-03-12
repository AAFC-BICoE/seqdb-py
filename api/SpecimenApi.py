'''
Created on Feb 16, 2016

@author: korolo

Class that extracts common functionality for all SeqDB API entities
'''

from api.BaseApiEntity import BaseApiEntity
from api.BaseSeqdbApi import BaseSeqdbApi
from api.BaseSeqdbApi import UnexpectedContent


class SpecimenApi(BaseApiEntity):

    def __init__(self, api_key, base_url, request_url):
        super(BaseApiEntity, self).__init__(api_key, base_url, request_url="specimen")
        
        
    @property
    def otherIdsFilter(self):
        return self.__nameFilter
    
    @otherIdsFilter.setter
    def otherIdsFilter(self, code, identifier):
        self.__otherIdsFilter = code + identifier
        
    def getParamsStr(self):
        params = "filterName=otherIds&filterValue={}&filterWildcard=false&".format(self.__otherIdsFilter)
        
        return params
    
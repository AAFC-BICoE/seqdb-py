'''
Created on Feb 16, 2016

@author: korolo

Class that extracts common functionality for all SeqDB API entities
'''

from api.BaseApiEntity import BaseApiEntity
from api.BaseSeqdbApi import BaseSeqdbApi
from api.BaseSeqdbApi import UnexpectedContent


class ProjectTagApi(BaseApiEntity):

    def __init__(self, api_key, base_url, request_url):
        super(BaseApiEntity, self).__init__(api_key, base_url, request_url="projectTag")
        
        
    @property
    def nameFilter(self):
        return self.__nameFilter
    
    @nameFilter.setter
    def nameFilter(self, name):
        self.__nameFilter = name
        
    def getParamsStr(self):
        params = "filterName=name&filterValue=%s&filterWildcard=false&" %self.__nameFilter
        
        return params
    

'''
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
'''

from api.BaseApiEntity import BaseApiEntity
from api.BaseSeqdbApi import BaseSeqdbApi
from api.BaseSeqdbApi import UnexpectedContent


class ProjectTagApi(BaseApiEntity):

    def __init__(self, api_key, base_url):
        super(ProjectTagApi, self).__init__(api_key=api_key, base_url=base_url, request_url="projectTag")
        
        
    @property
    def nameFilter(self):
        return self.__nameFilter
    
    @nameFilter.setter
    def nameFilter(self, name):
        self.__nameFilter = name
        
    def getParamsStr(self):
        params = ''
        if self.__nameFilter:
            params = "filterName=name&filterValue={}&filterWildcard=true&".format(self.__nameFilter)
        
        return params
    
    def clearAllFilters(self):
        self.__nameFilter = None

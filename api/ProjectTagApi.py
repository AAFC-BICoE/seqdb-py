'''
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
'''

from api.BaseApiEntity import BaseApiEntity

class ProjectTagApi(BaseApiEntity):

    def __init__(self, api_key, base_url):
        self.clearAllFilters()
        super(ProjectTagApi, self).__init__(api_key=api_key, base_url=base_url, request_url="projectTag")
        
        
    @property
    def nameFilter(self):
        return self.__nameFilter
    
    @nameFilter.setter
    def nameFilter(self, name):
        self.__nameFilter = name
        
    def getParamsStr(self):
        params = ''
        if self.nameFilter:
            params = "filterName=name&filterValue={}&filterWildcard=true&".format(self.nameFilter)
        
        return params
    
    def clearAllFilters(self):
        self.nameFilter = None
    
    '''    
    def getIds(self):
        return super(ProjectTagApi, self).getIds(self.getParamsStr())
    '''

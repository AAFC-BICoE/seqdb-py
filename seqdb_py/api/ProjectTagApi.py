"""
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
"""

from api.BaseApiEntity import BaseApiEntity


class ProjectTagApi(BaseApiEntity):

    def __init__(self, api_key, base_url):
        self.clear_all_filters()
        super(ProjectTagApi, self).__init__(api_key=api_key,
                                            base_url=base_url,
                                            request_url='projectTag'
                                            )
        
        
    @property
    def name_filter(self):
        return self.__name_filter
    
    @name_filter.setter
    def name_filter(self, name):
        self.__name_filter = name
        
    def get_param_str(self):
        params = ''
        if self.name_filter:
            params = 'filterName=name&filterValue={}&filterWildcard=true&'\
                .format(self.name_filter)
        
        return params
    
    def clear_all_filters(self):
        self.name_filter = None
    
    """    
    def get_ids(self):
        return super(ProjectTagApi, self).get_ids(self.get_param_str())
    """

"""
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
"""

from api.BaseApiEntity import BaseApiEntity
# from api.BaseSeqdbApi import BaseSeqdbApi
# from api.BaseSeqdbApi import UnexpectedContent


class SpecimenApi(BaseApiEntity):

    def __init__(self, api_key, base_url, specimen_request_url='specimen'):
        self.clear_all_filters()
        super(SpecimenApi, self).__init__(api_key=api_key, base_url=base_url,
                                          request_url=specimen_request_url)

    @property
    def other_ids_filter(self):
        return self.__other_ids_filter
    
    @other_ids_filter.setter
    def other_ids_filter(self, other_ids_filter):
        """
        other_ids_filter = code + identifier
        """
        self.__other_ids_filter = other_ids_filter
        
    def get_param_str(self):
        params = ''
        if self.other_ids_filter:
            params = 'filterName=otherIds&filterValue={}&filterWildcard=true&'\
                .format(self.other_ids_filter)
        
        return params
    
    def clear_all_filters(self):
        self.other_ids_filter = None

    '''    
    def get_ids(self):
        return super(SpecimenApi, self).get_ids(self.get_param_str())
    '''
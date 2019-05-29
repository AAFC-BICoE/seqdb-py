"""
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
"""

import json

from api.BaseApiEntity import BaseApiEntity
from api.BaseSeqdbApi import UnexpectedContent


class GeneRegionApi(BaseApiEntity):

    def __init__(self, api_key, base_url):
        self.clear_all_filters()
        super(GeneRegionApi, self).__init__(api_key=api_key, 
                                            base_url=base_url, 
                                            request_url='region')
    
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
    
    def get_its_region_ids(self, offset=0):
        """ Get region IDs of ITS sequences
        Args:
            offset: 0 if it is a first query, then number of records from
            which to load the next set of ids
        Returns:
            a list of seqdb its region ids 
            offset of results. If 0 then all/last set of results have been
            retrieved, if > 0, then the function has to be called again with
            this offset to retrieve more results
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        """
        
        its_region_names = ['ssu', '16s', '18s', 'its', '5.8s', 'lsu',
                            '23s', '25s', '28s', 'internal transcribed spacer']
        its_region_ids = set()
        
        for its_name in its_region_names:
            self.name_filter = its_name
            current_ids, offset = self.get_ids_with_offset(offset=offset)
            its_region_ids.update(current_ids)
            while offset:
                current_ids, offset = self.get_ids_with_offset(offset=offset)
                its_region_ids.update(current_ids)       
        
        self.clear_all_filters()
            
        return list(its_region_ids)
    
    def create(self, name, description, group_id=1):
        """ Creates a region
        Args:
            name: region name
            description: region description
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        """
        # TODO Question requirement for gene region to be associated with a
        # group
        # TODO Don't hard code group
        post_data = {
            'region': {
                'description': description,
                'group': {'id': group_id},
                'name': name}
            }

        resp = super(GeneRegionApi, self).create(
            '{}{}'.format(self.base_url, self.request_url),
            json.dumps(post_data))
        jsn_resp = resp.json()

        if 'result' and 'metadata' not in jsn_resp:
            raise UnexpectedContent(response=jsn_resp)
        
        if 'statusCode' and 'message' not in jsn_resp['metadata']:
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result']

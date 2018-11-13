"""
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
"""

import json

from api.BaseApiEntity import BaseApiEntity
from api.BaseSeqdbApi import UnexpectedContent


class FeatureApi(BaseApiEntity):

    def __init__(self, api_key, base_url):
        super(FeatureApi, self).__init__(api_key=api_key, base_url=base_url,
                                         request_url='feature')
        
    def get_param_str(self):
        return ''
    
    def create(self, name, feature_type_id, feature_locations, sequence_id,
               description='', feature_default=False, parent_id=None):
        """ Creates a Feature
        Args:
            name: name of the feature
        Returns:
            'result' from the json response
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.HTTPError
            UnexpectedContent
            :param description:
        """
        post_data = {
            'feature': {
                'name': name,
                'featureType': {'id': feature_type_id},
                'featureLocations': feature_locations,
                'description': description,
                'featureDefault': feature_default
            }
        }
        if parent_id is not None:
            post_data['feature']['parentFeature'] = {'id': parent_id}

        resp = super(FeatureApi, self)\
            .create('{}sequence/{}/{}'
                    .format(self.base_url, sequence_id, self.request_url),
                    json.dumps(post_data))
        jsn_resp = resp.json()

        if 'result' not in jsn_resp and \
                ('statusCode' and 'message' not in jsn_resp['metadata']):
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result']

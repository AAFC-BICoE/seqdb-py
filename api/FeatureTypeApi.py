'''
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
'''

import json

from api.BaseApiEntity import BaseApiEntity
from api.BaseSeqdbApi import UnexpectedContent


class FeatureTypeApi(BaseApiEntity):

    def __init__(self, api_key, base_url):
        super(FeatureTypeApi, self).__init__(api_key=api_key, base_url=base_url, request_url="featureType")
        
        
    def getParamsStr(self):
        return ''
    
    def getFeatureTypesWithIds(self):
        ''' 
        Returns:
            a dictionary of Feature types with feature name as keys and featureIds as values
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        feature_types = ''
        jsn_resp, result_offset = self.retrieveJsonWithOffset(request_url=self.request_url)
        

        if jsn_resp:
  
            feature_type_ids = jsn_resp['result']
            # Get the rest of the results, if all where not returned with the first query
            while result_offset:
                jsn_resp, result_offset = self.retrieveJsonWithOffset(request_url=self.request_url, 
                                                                      offset=result_offset) 
                feature_type_ids.extend(jsn_resp['result'])
            
            feature_types = {}

            for feat_type_id in feature_type_ids:
                jsn_resp = self.retrieveJson("{}/{}".format(self.request_url,feat_type_id))

                if jsn_resp:

                    feature_name = jsn_resp['result']['name']
                    feature_types[feature_name] = feat_type_id

        return feature_types
    
    
    def create(self, featureTypeName, featureTypeDescription=''):
        ''' Creates a FeatureType
        Returns:
            'result' from the json response
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        post_data = {"featureType": {
            "description": featureTypeDescription, "name": featureTypeName}}

        resp = super(FeatureTypeApi, self).create(
                "{}{}".format(self.base_url, self.request_url), 
                json.dumps(post_data))
        jsn_resp = resp.json()

        if ('result' not in jsn_resp and 
            ('statusCode' and 'message' not in jsn_resp['metadata'])):
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result']

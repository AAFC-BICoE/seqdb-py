'''
Created on April 06, 2016

@author: Oksana Korol
'''

import json

from api.BaseSeqdbApi import UnexpectedContent
from api.BaseSequenceEntity import BaseSequenceEntity


class ConsensusSequenceApi(BaseSequenceEntity):
    '''
    classdocs
    '''

    def __init__(self, api_key, base_url):
        super(ConsensusSequenceApi, self).__init__(api_key=api_key, base_url=base_url, request_url="consensus")
        
        
    # TODO verify which payload values are required by the SeqDB WS API
    def create(self, name, sequence, qualities=None,
            seqType="N", readNum=0, additional=None):
        
        post_data = {
            'consensus': {
                'name': name,
                'seq': sequence,
                'seqType': seqType,
                'readNum': readNum,
                'qualities': qualities,
            },
        }

        if additional is not None:
            post_data['consensus'].update(additional)

        resp = super(ConsensusSequenceApi, self).create(self.base_url + "/consensus", json.dumps(post_data))
        jsn_resp = resp.json()

        if 'result' and 'metadata' not in jsn_resp:
            raise UnexpectedContent(response=jsn_resp)
        
        if 'statusCode' and 'message' not in jsn_resp['metadata']:
            raise UnexpectedContent(response=jsn_resp)


        return jsn_resp['result'], jsn_resp['metadata']['statusCode'], jsn_resp['metadata']['message']

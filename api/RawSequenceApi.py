'''
Created on Apr 1, 2016

@author: Oksana Korol
'''
from api.BaseSeqdbApi import UnexpectedContent
from api.BaseSequenceEntity import BaseSequenceEntity


class RawSequenceApi(BaseSequenceEntity):
    '''
    classdocs
    '''


    def __init__(self, api_key, base_url):
        super(RawSequenceApi, self).__init__(api_key=api_key, base_url=base_url, request_url="sequence")
    
    # This method is only applicable to raw sequences, so it is not in BaseSequenceEntity
    def getFastqSequencesWithOffset(self, offset, limit):
        fastq_resp, result_offset = self._getSequencesWithOffset(offset, limit, "fastq")
        if fastq_resp and fastq_resp[0]!="@":
            raise UnexpectedContent("Response is not in fastq format.")
        
        return fastq_resp, result_offset
    

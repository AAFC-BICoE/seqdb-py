"""
Created on Apr 1, 2016

@author: Oksana Korol
"""
import base64
import gzip
import json
import logging
import mimetypes
import os

from BaseSeqdbApi import UnexpectedContent
from BaseSequenceEntity import BaseSequenceEntity


class RawSequenceApi(BaseSequenceEntity):
    """
    classdocs
    """

    def __init__(self, api_key, base_url):
        super(RawSequenceApi, self).\
            __init__(api_key=api_key,
                     base_url=base_url,
                     request_url='sequence'
                     )
    
    # This method is only applicable to raw sequences, so it is not in
    # BaseSequenceEntity
    def get_fastq_sequences_with_offset(self, offset, limit):
        fastq_resp, result_offset = self.get_sequences_with_offset(
            sequence_format='fastq', offset=offset, limit=limit)
        if fastq_resp and fastq_resp[0] != '@':
            raise UnexpectedContent('Response is not in fastq format.')
        
        return fastq_resp, result_offset
    
    def get_fastq_sequence(self, seq_id):
        """
        Returns a sequence in a fasta format.
            Note that to get multiple fasta sequences it is bests to use 
            getFastaSequencesWithOffset() method, to avoid multiple call
            to the SeqDB API.
        """
        return self._get_formatted_seq(seq_id, 'fastq')

    def import_chromat_sequences(
            self, blob, dest_file_name,
            notes='', trace_file_path=''):
        """ Imports a binary blob (i.e. chromatogram) to seqdb
        Kwargs:
            dest_file_name name with which the chromatogram will be saved on
                           sedDB side; i.e. expected to have .ab1 extension
        Returns:
            a list of seqdb sequence ids for the created sequences OR empty
            list if creation failed
        Raises:
            UnexpectedContent
        """

        # Base 64-encode chromat as required by SeqDB api
        chromat_b64 = base64.b64encode(blob)
        chromat_b64_ascii = chromat_b64.decode('ascii')

        post_data = {
            'sequenceImportPayload': {
                'base64File': chromat_b64_ascii,
                'fileName': dest_file_name,
                'plateType': 1,
                'createLocation': False,
                'traceFilePath': trace_file_path,
                'notes': notes
            }
        }
        resp = self.create(
            self.base_url + 'sequenceImport', json.dumps(post_data))
        jsn_resp = resp.json()
        
        if 'statusCode' and 'message' not in jsn_resp['metadata']:
            raise UnexpectedContent(response=jsn_resp)

        result = ''
        if 'result' in jsn_resp:
            if type(jsn_resp['result']) == list and \
                    len(jsn_resp['result']) != 1:
                raise UnexpectedContent('When creating a chromatogram '
                                        'sequence, response should contain '
                                        'only one result.')
            result = jsn_resp['result'][0]
        else:
            logging.warning(
                'Importing chromatogram failed with status code: {} '
                'and message: ()'.format(
                    jsn_resp['metadata']['statusCode'],
                    jsn_resp['metadata']['message']))

        return result

    def import_chromat_sequences_from_file(
            self, chromat_file, notes='',
            trace_file_path='', dest_file_name=''):
        """ Imports a chromatogram from a file
        Args:
            chromat_file name of the chromatogram file
        Returns:
            list of sequence ids of the sequences that were imported from the
            chromatogram
        Raises:
            IOError
            UnexpectedContent
        """

        if not os.path.isfile(chromat_file):
            raise IOError('Expecting a file, but got a directory.')

        chromat_file_name = os.path.basename(chromat_file)

        if not dest_file_name:
            # get just the file name, no extension
            dest_file_name = os.path.splitext(
                os.path.basename(chromat_file_name))[0]

        if mimetypes.guess_type(chromat_file_name)[1] == 'gzip':

            file_stream = gzip.open(chromat_file, 'rb')

        else:
            file_stream = open(chromat_file, 'rb')

        blob = file_stream.read()

        file_stream.close()

        return self.import_chromat_sequences(
                blob=blob, dest_file_name=dest_file_name,
                notes=notes, trace_file_path=trace_file_path)

    # Note that this method should be in BaseSequenceEntity, but because
    # ../region/123/consensus is not implemented, it has been moved here

    def get_sequence_ids_by_region_with_offset(self, region_id, offset=0):
        """ Given a region id, return sequence ids, belonging to this region
        Returns:
            a list of seqdb sequence ids for the created sequences OR empty
            list id created failed
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        """
        
        region_request_url = '/region/' + str(region_id) + '/' + self.request_url
        jsn_resp, result_offset = self.retrieve_json_with_offset(
            request_url=region_request_url, offset=offset)
        
        sequence_ids = ''
        
        if jsn_resp:            
            sequence_ids = jsn_resp['result']            
        
        return sequence_ids, result_offset

    def get_accepted_specimen_determination(self, sequence_id):
        """
        Returns:
             an accepted determination of a SPECIMEN, which is associated
             with this sequence
        """

        jsn_resp = self.retrieve_json('{}/{}/acceptedSpecimenDetermination'
                                      .format(self.request_url, sequence_id))
        
        if jsn_resp:
            return jsn_resp['result']
        else:
            return ''

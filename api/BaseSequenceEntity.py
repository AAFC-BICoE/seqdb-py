'''
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
'''

import base64
import gzip
import json
import logging
import mimetypes
import os

from api.BaseApiEntity import BaseApiEntity
from api.BaseSeqdbApi import UnexpectedContent
from api.ProjectTagApi import ProjectTagApi


class BaseSequenceEntity(BaseApiEntity):

    def __init__(self, api_key, base_url, request_url):
        self.clearAllFilters()
        super(BaseSequenceEntity, self).__init__(api_key=api_key, base_url=base_url, request_url=request_url)
        
    @property
    def specimenNumFilter(self):
        return self.__specimenNumFilter
    
    @specimenNumFilter.setter
    def specimenNumFilter(self, specimenNum):
        self.__specimenNumFilter = specimenNum

    @property
    def sequenceNameFilter(self):
        return self.__sequenceNameFilter
    
    @sequenceNameFilter.setter
    def sequenceNameFilter(self, sequenceName):
        self.__sequenceNameFilter = sequenceName
    
    @property
    def sampleNameFilter(self):
        return self.__sampleNameFilter
    
    @sampleNameFilter.setter
    def sampleNameFilter(self, sampleName):
        self.__sampleNameFilter = sampleName
    
    @property
    def pubRefSeqFilter(self):
        return self.__pubRefSeqFilter
    
    @pubRefSeqFilter.setter
    def pubRefSeqFilter(self, isSubmittedToInsdc):
        if not(isSubmittedToInsdc) or isinstance(isSubmittedToInsdc, bool):
            self.__pubRefSeqFilter = isSubmittedToInsdc
        else:
            pass
            #TODO handle exceptions
    
    @property
    def genBankGIFilter(self):
        return self.__genBankGIFilter
    
    @genBankGIFilter.setter
    def genBankGIFilter(self, genBankGI):
        self.__genBankGIFilter = genBankGI
    
    @property
    def regionNameFilter(self):
        return self.__regionNameFilter
    
    @regionNameFilter.setter
    def regionNameFilter(self, regionName):
        self.__regionNameFilter = regionName
    
    @property
    def projectNameFilter(self):
        return self.__projectNameFilter
    
    @projectNameFilter.setter
    def projectNameFilter(self, projectName):
        self.__projectNameFilter = projectName
    
    @property
    def collectionCodeFilter(self):
        return self.__collectionCodeFilter
    
    @collectionCodeFilter.setter
    def collectionCodeFilter(self, collectionCode):
        self.__collectionCodeFilter = collectionCode
    
    @property
    def taxonomyRankFilter(self):
        return self.__taxonomyRankFilter
    
    @taxonomyRankFilter.setter
    def taxonomyRankFilter(self, taxonomyRank):
        self.__taxonomyRankFilter = taxonomyRank
    
    @property
    def taxonomyValueFilter(self):
        return self.__taxonomyValueFilter
    
    @taxonomyValueFilter.setter
    def taxonomyValueFilter(self, taxonomyValue):
        self.__taxonomyValueFilter = taxonomyValue
    
        
    def getParamsStr(self):
        params = ''

        if self.specimenNumFilter:
            params = params + "filterName=specimen.number&filterValue={}&filterWildcard=false&".format(self.specimenNumFilter)

        if self.sequenceNameFilter:
            params = params + "filterName=sequence.name&filterValue={}&filterWildcard=true&".format(self.sequenceNameFilter)
            
        if self.sampleNameFilter:
            params = params + "filterName=sample.name&filterValue={}&filterWildcard=true&".format(self.sampleNameFilter)
            
        if self.pubRefSeqFilter:
            params = params + "filterName=sequence.submittedToInsdc&filterValue=true&filterWildcard=false&"

        if self.genBankGIFilter:
            params = params + "filterName=sequence.genBankGI&filterValue={}&filterWildcard=false&".format(self.genBankGIFilter)
        
        if self.regionNameFilter:
            params = params + "filterName=region.name&filterValue={}&filterWildcard=true&".format(self.regionNameFilter)
            
        if self.projectNameFilter:
            projectTag = ProjectTagApi(api_key=self.api_key, base_url=self.base_url)
            projectTag.nameFilter = self.projectNameFilter
            project_ids = projectTag.getIds()
            
            if len(project_ids) == 1:
                params = params + "tagId=%s&" %project_ids[0]
                
            elif len(project_ids) > 1:
                raise "More than one project is associated with the name '%s'. Currently sequences can only be filtered on one project only. Please refine your search."
            
        if self.collectionCodeFilter:
            params = params + "filterName=biologicalCollection.name&filterValue={}&filterWildcard=false&".format(self.collectionCodeFilter)
            
        if self.taxonomyRankFilter and self.taxonomyValueFilter:
            filter_name = self.convertNcbiToSeqdbTaxRank(self.taxonomyRankFilter)
            params = params + "filterName=specimen.identification.taxonomy.{}&filterValue={}&filterOperator=and&filterWildcard=true&".format(filter_name, self.taxonomyValueFilter)
        
        return params
    
    def clearAllFilters(self):
        self.specimenNumFilter = None
        self.sequenceNameFilter = None
        self.sampleNameFilter = None
        self.pubRefSeqFilter = None
        self.genBankGIFilter = None
        self.regionNameFilter = None
        self.projectNameFilter = None
        self.collectionCodeFilter = None
        self.taxonomyRankFilter = None
        self.taxonomyValueFilter = None
    
    '''    
    def getNumber(self):
        return super(BaseSequenceEntity, self).getNumber(self.getParamsStr())
    
    def getIds(self):
        return super(BaseSequenceEntity, self).getIds(self.getParamsStr())
    '''
    
    def getFastaSequencesWithOffset(self, offset, limit):
        fasta_resp, result_offset = self._getSequencesWithOffset(offset, limit, "fasta")
        if fasta_resp and fasta_resp[0]!=">":
            raise UnexpectedContent("Response is not in fasta format.")
        
        return fasta_resp, result_offset
        
    
    def getFastqSequencesWithOffset(self, offset, limit):
        #TODO: Only raw sequences? Then move this method to SequenceApi
        fastq_resp, result_offset = self._getSequencesWithOffset(offset, limit, "fastq")
        if fastq_resp and fastq_resp[0]!="@":
            raise UnexpectedContent("Response is not in fastq format.")
        
        return fastq_resp, result_offset
        
    
    def _getSequencesWithOffset(self, offset, limit, sequence_format):
        
        if offset < 0:
            raise "Negative offset: either you've retrieved all sequences or the method usage is incorrect."
        
        if sequence_format not in {"fasta", "fastq"}:
            raise SystemError("Sequences can only be in fasta or fastq format.")
        
        #TODO: refactor this method so that the api call to get number of sequences
        #    is not called eache time this method is called
        # Define a global var for total number of sequences to avoid calling the 
        # api method every time __getSequencesWithOffset is called
        '''
        global globalSequenceNumber
        if not globalSequenceNumber:
            globalSequenceNumber = self.getNumber()
        '''
        
        params = self.getParamsStr() + "limit={}&offset={}&".format(limit,offset)
        
        resp = self.retrieve("{}/sequence.{}".format(self.base_url,sequence_format), params=params)
        
        """
        jsn_response, result_offset = self.retrieveJsonWithOffset(request_url=self.request_url + ".fasta", 
                                                   params=self.getParamsStr(), 
                                                   offset=offset, 
                                                   limit=limit)
        """
        
        fasta_resp = resp.content
        
        # Calculate the new offset
        result_offset = -1
        #TODO: refactor for efficiency, see TODO above.
        result_total_num = self.getNumber()
        result_returned_so_far =  limit + offset
        if (result_total_num > result_returned_so_far):    
            result_offset = result_returned_so_far 
        
        return fasta_resp, result_offset
         
    def getFastaSequence(self, seqId):
        ''' Returns a sequence in a fasta format.
            Note that to get multiple fasta sequences it is bests to use 
            getFastaSequencesWithOffset() method, to avoid multiple call to the SeqDB API.
        '''
        return self._getFormattedSeq(seqId, "fasta")

    def getFastqSequence(self, seqId):
        ''' Returns a sequence in a fasta format.
            Note that to get multiple fasta sequences it is bests to use 
            getFastaSequencesWithOffset() method, to avoid multiple call to the SeqDB API.
        '''
        return self._getFormattedSeq(seqId, "fastq")
    
    def _getFormattedSeq(self, seq_id, sequence_format):
        ''' Gets sequence in fasta or fastq format 
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
        '''
        if sequence_format not in {"fasta", "fastq"}:
            raise SystemError("Sequences can only be in fasta or fastq format.")
        
        url = "{}/{}/{}.{}".format(self.base_url, self.request_url, str(seq_id), sequence_format)
        response = self.retrieve(url)
        return response.content

    
    def importChromatSequences(
            self, blob, dest_file_name,
            notes="", trace_file_path=""):
        ''' Imports a binary blob (i.e. chromatogram) to seqdb
        Kwargs:
            dest_file_name name with which the chromatogram will be saved on
                           sedDB side; i.e. expected to have .ab1 extension
        Returns:
            a list of seqdb sequence ids for the created sequences OR empty
            list if creation failed
        Raises:
            UnexpectedContent
        '''

        chromat_b64 = base64.b64encode(blob)

        post_data = {
            "sequenceImportPayload": {
                "base64File": chromat_b64,
                "fileName": dest_file_name,
                "plateType": 1,
                "createLocation": False,
                "traceFilePath": trace_file_path,
                "notes": notes

            }
        }

        resp = self.create(
            self.base_url + '/sequenceImport', json.dumps(post_data))
        jsn_resp = resp.json()

        if 'statusCode' and 'message' not in jsn_resp['metadata']:
            raise UnexpectedContent(response=jsn_resp)

        result = ""
        if 'result' in jsn_resp:
            if type(jsn_resp['result']) == list and len(jsn_resp['result']) != 1:
                raise UnexpectedContent("When creating a chromatogram sequence, response should contain only one result.")
            result = jsn_resp['result'][0]
        else:
            logging.warn(
                "Importing chromatogram failed with status code: '%s' "
                "and message: '%s'" % (
                    jsn_resp['metadata']['statusCode'], jsn_resp['metadata']['message']))

        return result


    def importChromatSequencesFromFile(
            self, chromat_file, notes="",
            trace_file_path="", dest_file_name=""):
        ''' Imports a chromatogram from a file
        Args:
            chromat_file name of the chromatogram file
        Returns:
            list of sequence ids of the sequences that were imported from the
            chromatogram
        Raises:
            IOError
            UnexpectedContent
        '''

        if not os.path.isfile(chromat_file):
            raise IOError("Expecting a file, but got a directory.")

        chromat_file_name = os.path.basename(chromat_file)

        if not dest_file_name:
            # get just the file name, no extension
            dest_file_name = os.path.splitext(
                os.path.basename(chromat_file_name))[0]

        if mimetypes.guess_type(chromat_file_name)[1] == "gzip":

            file_strem = gzip.open(chromat_file, "rb")

            # Remove .gz extension for file name
            chromat_file_name = chromat_file_name[:-3]

        else:
            file_strem = open(chromat_file, "r")

        blob = file_strem.read()

        file_strem.close()

        return self.importChromatSequences(
                blob=blob, dest_file_name=dest_file_name,
                notes=notes, trace_file_path=trace_file_path)
        
        
    ########################################################################
    #                     Taxonomy Helpers
    ########################################################################
    
    # This is a conversion dictionary for NCBI taxonomy to be imported to SeqDB. 
    # Keys are seqdb taxonomy ranks, and values are NCBI
    _seqdb_to_ncbi_taxonomy = {'kingdom':'kingdom', 
                              'phylum':'phylum', 
                              'taxanomicClass':'class',
                              'taxanomicOrder':'order',
                              'family':'family',
                              'genus':'genus', 
                              'subgenera':'subgenus', 
                              'species':'species', 
                              'variety':'varietas',
                            }
    # NOTE that following SeqDB taxonomic ranks currently do not have equivalent in NCBI taxonomy:
    # 'strain','authors','division','synonym','typeSpecimen','taxanomicGroup','commonName','state'
    
    
    def vetTaxonomy(self, taxonomy):
        ''' Checks taxonomy dictionary and returns it in the format that can be accepted by 
            seqDB (i.e. appropriate ranks)
        '''
        vettedTaxonomy = dict()
        if taxonomy:
            for seqdb_rank in self._seqdb_to_ncbi_taxonomy:
                ncbi_rank = self._seqdb_to_ncbi_taxonomy[seqdb_rank]
                if ncbi_rank in taxonomy:
                    vettedTaxonomy[seqdb_rank] = taxonomy[ncbi_rank]
        
        return vettedTaxonomy
    
    
    def convertNcbiToSeqdbTaxRank(self, tax_rank):
        for seqdb_tax_rank in self._seqdb_to_ncbi_taxonomy:
            if self._seqdb_to_ncbi_taxonomy[seqdb_tax_rank] == tax_rank:
                return seqdb_tax_rank
        return "No matches found."

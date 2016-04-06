'''
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
'''


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
    
    
    def getFastaSequencesWithOffset(self, offset, limit=20):
        fasta_resp, result_offset = self._getSequencesWithOffset(offset, limit, "fasta")
        if fasta_resp and fasta_resp[0]!=">":
            raise UnexpectedContent("Response is not in fasta format.")
        
        return fasta_resp, result_offset
            
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
        
        resp = self.retrieve("{}/{}.{}".format(self.base_url, self.request_url, sequence_format), params=params)
        
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

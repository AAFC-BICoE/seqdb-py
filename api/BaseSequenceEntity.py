'''
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
'''

from api.BaseApiEntity import BaseApiEntity
from api.BaseSeqdbApi import BaseSeqdbApi
from api.BaseSeqdbApi import UnexpectedContent


class BaseSequenceEntity(BaseApiEntity):

    def __init__(self, api_key, base_url):
        super(BaseSequenceEntity, self).__init__(api_key=api_key, base_url=base_url, request_url="sequence")
        
        
    @property
    def specimenNumFilter(self):
        return self.specimenNumFilter
    
    @specimenNumFilter.setter
    def specimenNumFilter(self, specimenNum):
        self.specimenNumFilter = specimenNum
    
    @property
    def sequenceNameFilter(self):
        return self.sequenceNameFilter
    
    @sequenceNameFilter.setter
    def sequenceNameFilter(self, sequenceName):
        self.sequenceNameFilter = sequenceName
    
    @property
    def sampleNameFilter(self):
        return self.specimenNumFilter
    
    @sampleNameFilter.setter
    def sampleNameFilter(self, sampleName):
        self.sampleNameFilter = sampleName
    
    @property
    def pubRefSeqFilter(self):
        return self.pubRefSeqFilter
    
    @pubRefSeqFilter.setter
    def pubRefSeqFilter(self, isSubmittedToInsdc):
        if isinstance(isSubmittedToInsdc, bool):
            self.pubRefSeqFilter = isSubmittedToInsdc
        else:
            pass
            #TODO handle exceptions
    
    @property
    def genBankGIFilter(self):
        return self.genBankGIFilter
    
    @genBankGIFilter.setter
    def genBankGIFilter(self, genBankGI):
        self.genBankGIFilter = genBankGI
    
    @property
    def regionNameFilter(self):
        return self.regionNameFilter
    
    @regionNameFilter.setter
    def regionNameFilter(self, regionName):
        self.regionNameFilter = regionName
    
        
    def getParamsStr(self):
        params = ''
            
        if self.specimenNumFilter:
            params = params + "filterName=specimen.number&filterValue={}&filterWildcard=false&".format(self.specimenNumFilter)
    
        if self.sequenceNameFilter:
            params = params + "filterName=sequence.name&filterValue={}&filterWildcard=true&".format(self.sequenceName)
            
        if self.sampleNameFilter:
            params = params + "filterName=sample.name&filterValue=%s&filterWildcard=true&".format(self.sampleName)
            
        if self.pubRefSeqFilter:
            params = params + "filterName=sequence.submittedToInsdc&filterValue=true&filterWildcard=false"

        if self.genBankGIFilter:
            params = params + "filterName=sequence.genBankGI&filterValue=%s&filterWildcard=false&".format(self.genBankGI)
        
        if self.regionNameFilter:
            params = params + "filterName=region.name&filterValue=%s&filterWildcard=true&".format(self.regionName)
            
        '''    
        if self.projectNameFilter:
            project_ids = self.getProjectTagIds(projectName)
            
            if len(project_ids) == 1:
                params = params + "tagId=%s&" %project_ids[0]
                
            elif len(project_ids) > 1:
                raise "More than one project is associated with the name '%s'. Currently sequences can only be filtered on one project only. Please refine your search."
            
        if collectionCodeFilter:
            params = params + "filterName=biologicalCollection.name&filterValue=%s&filterWildcard=false&".format(self.collectionCode
            
        if taxonomyRankFilter:
            filter_name = self.convertNcbiToSeqdbTaxRank(taxonomyRank)
            params = params + "filterName=specimen.identification.taxonomy.%s&filterValue=%s&filterOperator=and&filterWildcard=true&" %(filter_name, taxonomyValue)
        '''
        
        return params
    
    def clearAllFilters(self):
        self.specimenNumFilter = None
        self.sequenceNameFilter = None
        self.sampleNameFilter = None
        self.pubRefSeqFilter = None
        self.genBankGIFilter = None
        self.regionNameFilter = None
        
'''
Created on Apr 1, 2016

@author: korolo
'''
import unittest

import yaml

from config import config_root
from api.ConsensusSequenceApi import ConsensusSequenceApi


class TestConsensusSequenceApi(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = ConsensusSequenceApi(api_key=config['seqdb']['api_key'],
                                                   base_url=config['seqdb']['api_url'])
    

    def setUp(self):
        pass


    def tearDown(self):
        pass
    
    def testFilters(self):
        # specimenNumber - tested below, skipping
        # sequenceName
        self.fixture.sequenceNameFilter = "Pyt_arrhenomanes_BR0"
        self.assertEqual([358381,358485], self.fixture.getIds(), "Expecting 2 consensus sequences filtered by sequenceName = Pyt_arrhenomanes_BR0")
        self.fixture.clearAllFilters()
        
        # sampleNameFilter
        self.fixture.sampleNameFilter = "LEV5508"
        self.assertEqual(1, self.fixture.getNumber(), 
                         "Expecting 1 consensus sequences filtered by sampleName = LEV5508, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()
        
        # pubRefSeqFilter
        #TODO
        '''
        self.fixture.pubRefSeqFilter = True
        #TODO: this fails, i.e. curl -H "apikey: ***REMOVED***" "***REMOVED***/sequence?filterName=sequence.submittedToInsdc&filterValue=true&filterWildcard=false"
        # Investigate why
        va = self.fixture.getNumber()
        self.assertEqual(15, self.fixture.getNumber(), 
                         "Expecting 15 raw sequences filtered by pubRefSeq = , but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()
        '''
        
        # genBankGIFilter
        #TODO: fix this test
        '''
        self.fixture.genBankGIFilter = "gi_"
        self.assertEqual(15, self.fixture.getNumber(), 
                         "Expecting 15 raw sequences filtered by genBankGI = LEV4277, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()
        '''
        
        # regionNameFilter
        self.fixture.regionNameFilter = "ITS2-28S"
        self.assertEqual(4, self.fixture.getNumber(), 
                         "Expecting 4 consensus sequences filtered by regionName = ITS2-28S, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()
        
        # collectionCodeFilter
        self.fixture.collectionCodeFilter = "lev"
        self.assertEqual(211, self.fixture.getNumber(), 
                         "Expecting 211 consensus sequences filtered by collectionCode = lev, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()
        
        # taxonomyRankFilter
        self.fixture.taxonomyRankFilter = "species"
        self.fixture.taxonomyValueFilter = "megasperma"
        self.assertEqual(3, self.fixture.getNumber(), 
                         "Expecting 3 consensus sequences filtered by taxonomy, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()
        
        # TODO: test combinations of filters
        

    def testGetConsensusSequenceIds(self):
        self.fixture.specimenNumFilter = 4405
        actual = self.fixture.getIds()
        self.assertTrue(actual, "No Sequence ids returned.")
        self.assertEqual(1, len(actual),"Expecting 1 consensus sequence associated with this specimen.")
        self.assertIn(358301, actual, "Sequence id 358301 is expected to be in results, since it is consensus.")  
        self.assertNotIn(27755, actual, "Sequence id 27755 is not expected to be in results, since it is consensus.")
        
    def testCreateGetDeleteSequence(self):
        # create
        seqId, errCod, msg = self.fixture.create(name="Test", sequence="ACGTCTGATCGATC")
        self.assertTrue(seqId, "Creating consensus sequence did not return an id.")
        self.assertEqual(errCod, 201, "Did not get successful exit code for create consensus sequence.")
        
        # get
        self.fixture.sequenceNameFilter="Test"
        seqIds = self.fixture.getIds()
        self.assertTrue(seqIds, "Creating consensus sequence did not return an id.")
        self.assertIn(seqId, seqIds, "Expected sequence id was not in the result.")
        
        # delete
        delete_jsn_resp = self.fixture.delete(seqId)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete feature type.")
    
        
    def testGetFastaSequencesWithOffset(self):
        self.fixture.specimenNumFilter = 4405
        actual_fasta, result_offset = self.fixture.getFastaSequencesWithOffset(offset=0)
        self.assertTrue(actual_fasta, "No Sequences returned.")
        self.assertIn(">seqdb|358301", actual_fasta, "Fasta return is expected to have sequence 358301, since it is consensus.")
        self.assertNotIn(">seqdb|27755", actual_fasta,"Fasta return is not expected to have sequence id 27755, since it is raw.")

            
    def testGetFastaSeq(self):
        actual = self.fixture.getFastaSequence("358301")
        self.assertTrue(actual, "Fasta sequence is empty.")
        self.assertIn(">seqdb|358301", actual, "Fasta does not contain >seqdb|358301.")
              
            

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
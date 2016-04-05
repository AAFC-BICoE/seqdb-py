'''
Created on Apr 1, 2016

@author: korolo
'''
import unittest

import yaml

from api.SequenceApi import SequenceApi
from config import config_root


class TestSequenceApi(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = SequenceApi(api_key=config['seqdb']['api_key'],
                                                   base_url=config['seqdb']['api_url'])
    

    def setUp(self):
        pass


    def tearDown(self):
        pass
    
    def testGetFastaSequencesWithOffset(self):
        self.fixture.specimenNumFilter = 4405
        # This filter should have 22 raw sequences associated with it. Therefore there should
        # be 2 calls with limit 15
        actual_fasta, result_offset = self.fixture.getFastaSequencesWithOffset(offset=0, limit=15)
        self.assertTrue(actual_fasta, "No Sequences returned.")
        self.assertIn(">seqdb|27755", actual_fasta,"Expecting that fasta return will contain id 27755.")
        self.assertNotIn(">seqdb|358301", actual_fasta, "Fasta return is not expected to have sequence 358301, since it is consensus.")
        self.assertEquals(result_offset, 15, "First offset should be 15 (there are 22 raw sequences for specimenNum=4405).")
        
        actual_fasta, result_offset = self.fixture.getFastaSequencesWithOffset(offset=result_offset, limit=15)
        self.assertTrue(actual_fasta, "No Sequences returned.")
        self.assertIn(">seqdb|160946", actual_fasta,"Expecting that fasta return will contain id 160946.")
        self.assertNotIn(">seqdb|358301", actual_fasta, "Fasta return is not expected to have sequence 358301, since it is consensus.")
        self.assertEquals(result_offset, -1, "Second offset should be -1, since there are no more sequences to retrieve.")
        
    def testGetSequenceIds(self):
        self.fixture.specimenNumFilter = 4405
        actual = self.fixture.getIds()
        self.assertTrue(actual, "No Sequence ids returned.")
        self.assertEqual(22, len(actual),"Expecting 22 sequences associated with this specimen.")
        self.assertIn(27755, actual, "Sequence id 27755 is expected to be associated with specimen 4405.")
        self.assertNotIn(358301, actual, "Sequence id 358301 is not expected to be in results, since it is consensus.")  

    
    def testCreateChromatSequence_wrong_path(self):
        self.assertRaises(IOError, self.fixture.importChromatSequencesFromFile, "data/")
        self.assertRaises(IOError, self.fixture.importChromatSequencesFromFile, "zzz/non-existent.ab1")

    def testCreateDeleteChromatSequence(self):
        """ Test creating a sequence with binary .abi or .ab1 file (chromatogram) """
        seq_id = self.fixture.importChromatSequencesFromFile(chromat_file = "data/GRDI_test_seq.ab1", notes="This is a test upload.", 
                    trace_file_path="test_path",dest_file_name="test.ab1",)
        
        self.assertTrue(seq_id, "Persisting chromatogram did not return an id.")
    
        # Delete
        delete_jsn_resp = self.fixture.delete(seq_id)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete feature type.")

    

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
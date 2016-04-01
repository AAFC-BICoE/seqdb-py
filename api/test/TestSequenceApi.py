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
    
    '''
    def testGetFastaSequencesWithOffset(self):
        actual = self.fixture.getRawSequencesWithOffset(offset=0, limit=30, sequence_format="fasta", specimenNum=4405)
        self.assertTrue(actual, "No Sequences returned.")
        self.assertIn(">seqdb|27755", actual,"Expecting that fasta return will contain id 27755.")
        self.assertNotIn(">seqdb|358301", actual, "Fasta return is not expected to have sequence 358301, since it is consensus.")
    '''
        
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
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
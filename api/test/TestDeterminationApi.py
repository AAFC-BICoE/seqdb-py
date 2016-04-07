'''
Created on Apr 1, 2016

@author: korolo
'''
import unittest

import yaml

from config import config_root
from api.DeterminationApi import DeterminationApi


class TestDeterminationApi(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = DeterminationApi(api_key=config['seqdb']['api_key'],
                                                   base_url=config['seqdb']['api_url'])
    

    def setUp(self):
        pass


    def tearDown(self):
        pass
    

    def testCreateGetDeleteDetermination(self):
        #test_taxonomy = {'superkingdom': 'Eukaryota', 'genus': 'Phytophthora', 'species': 'Phytophthora ramorum', 'order': 'Peronosporales'}
        test_taxonomy = {'genus': 'Phytophthora', 'species': 'Phytophthora ramorum', 'taxanomicOrder': 'Peronosporales'}
        det_id = self.fixture.createSequenceDetermination(28954, test_taxonomy)
        self.assertTrue(det_id, "Creating determination on a sequence did not return a determination id.")
        
        get_determination = self.fixture.getEntity(det_id)
        self.assertTrue(get_determination, "Did not get determination back.")
        self.assertEquals(False, get_determination['accepted'], "")
        self.assertEquals('Phytophthora', get_determination['taxonomy']['genus'], "")
        
        # delete
        delete_jsn_resp = self.fixture.delete(det_id)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete determination.")
    
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
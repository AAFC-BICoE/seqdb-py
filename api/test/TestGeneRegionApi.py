'''
Created on Apr 1, 2016

@author: korolo
'''
import unittest

import yaml

from config import config_root
from api.GeneRegionApi import GeneRegionApi


class TestGeneRegionApi(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = GeneRegionApi(api_key=config['seqdb']['api_key'],
                                                   base_url=config['seqdb']['api_url'])
    

    def setUp(self):
        pass


    def tearDown(self):
        pass
    
    ''' Test not working" "Missing or malformed request object". Following up with Nazir.
    '''
    def testCreateGetDeleteITSRegion(self):
        # create
        create_regId = self.fixture.create("ITS_test", "Test", 168)
        self.assertTrue(create_regId, "Creating region did not return an id.")
        
        # get
        get_regIds = self.fixture.getItsRegionIds()
        self.assertTrue(get_regIds, "No ITS region ids returned.")
        self.assertIn(create_regId, get_regIds, "Expected region is not in the results.")

        # delete
        delete_jsn_resp = self.fixture.delete(create_regId)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete region.")     
        
        # get
        get_regIds = self.fixture.getItsRegionIds()
        self.assertNotIn(create_regId, get_regIds, "Not expecting this region to be in the results. Region was deleted.")
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
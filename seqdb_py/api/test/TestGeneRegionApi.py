"""
Created on Apr 1, 2016

@author: korolo
"""
import unittest

import yaml

from config import config_root
from api.GeneRegionApi import GeneRegionApi


class TestGeneRegionApi(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        with open(config_root.path() +
                  '/config4tests.yaml', 'r') as config_file:
            config = yaml.safe_load(config_file)
            cls.fixture = GeneRegionApi(
                api_key=config['seqdb']['api_key'],
                base_url=config['seqdb']['base_url'])

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    ''' Test not working' 'Missing or malformed request object'. Following up with Nazir.
    '''
    def test_create_get_delete_its_region(self):
        # create
        create_reg_id = self.fixture.create('ITS_test', 'Test', 168)
        self.assertTrue(create_reg_id, 'Creating region did not return an id.')
        
        # get
        get_reg_ids = self.fixture.get_its_region_ids()
        self.assertTrue(get_reg_ids, 'No ITS region ids returned.')
        self.assertIn(create_reg_id, get_reg_ids,
                      'Expected region is not in the results.')

        # delete
        delete_jsn_resp = self.fixture.delete(create_reg_id)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'],
                         'Could not delete region.')
        
        # get
        get_reg_ids = self.fixture.get_its_region_ids()
        self.assertNotIn(create_reg_id, get_reg_ids,
                         'Not expecting this region to be in '
                         'the results. Region was deleted.')


if __name__ == '__main__':
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

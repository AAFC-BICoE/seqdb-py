"""
Created on Mar 12, 2015

Note: these tests should be run against a local SeqDB instance only. They do not mock the API. 
    Since the work on SeqDB API is ongoing, they serve as both tests for the request methods 
    as well as tests for the SeqDB API itself.
    
    Therefore, these tests are brittle. They depend not only on the current content of your 
    local SeqDB, but also on the information that you have access to (i.e. groups you're 
    part of). Here, by 'you' I mean the user whose API Key is used for testing (seqdb_api_key)

Before running the tests:
    - Make suse that config/config4tests.yaml is created from sample and the modified to point
     to local instance of SeqDB
    - Make sure that local instance of SeqDB is running

@author: Oksana Korol
"""

import requests
import unittest 
import yaml

from api.ProjectTagApi import ProjectTagApi
from config import config_root


class TestProjectTagApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(config_root.path() +
                  '/config4tests.yaml', 'r') as config_file:
            config = yaml.safe_load(config_file)
            cls.fixture = ProjectTagApi(
                api_key=config['seqdb']['api_key'],
                base_url=config['seqdb']['base_url'])
    
    def test_retrieve(self):
        # Test faulty connection
        self.assertRaises(requests.exceptions.ConnectionError,
                          self.fixture.retrieve, 'http://jibberish')
        # TODO: test wrong api key
        
    def test_retrieve_json(self):
        actual = self.fixture.retrieve_json('/projectTag/12')
        self.assertTrue(actual, 'retrieve_json: no result was returned.')
        # self.assertIn('Cranberry', actual,
        # 'Expecting a Cranberry project tag.')
        # Cranberry
        
    def test_get_entity(self):
        actual = self.fixture.get_entity(12)
        self.assertTrue(actual, 'No Project Tag returned.')
        self.assertIn('Cranberry', actual['name'],
                      'Expecting Cranberry project tag.')
        
    def test_get_ids_with_offset(self):
        self.fixture.name_filter = 'a'
        actual_entity_ids = self.fixture.get_ids()
        self.assertEquals(3, len(actual_entity_ids),
                          'Expecting 3 ids, but got {}.'
                          .format(len(actual_entity_ids)))


if __name__ == '__main__':
    unittest.main()

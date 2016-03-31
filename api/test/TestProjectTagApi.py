'''
Created on Mar 12, 2015

Note: these tests should be run against a local SeqDB instance only. They do not mock the API. 
    Since the work on SeqDB API is ongoing, they serve as both tests for the request methods 
    as well as tests for the SeqDB API itself.
    
    Therefore, these tests are brittle. They depend not only on the current content of your 
    local SeqDB, but also on the information that you have access to (i.e. groups you're 
    part of). Here, by "you" I mean the user whose API Key is used for testing (seqdb_api_key)

Before running the tests:
    - Make suse that config/config4tests.yaml is created from sample and the modified to point
     to local instance of SeqDB
    - Make sure that local instance of SeqDB is running

@author: Oksana Korol
'''
import requests
import unittest 
import yaml

from api import ProjectTagApi
from config import config_root


class TestProjectTagApi(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = ProjectTagApi.ProjectTagApi(api_key=config['seqdb']['api_key'],
                                                   base_url=config['seqdb']['api_url'])
    
    def testRetrieve(self):
        # Test faulty connection
        self.assertRaises(requests.exceptions.ConnectionError, self.fixture.retrieve, "http://jibberish")
        # TODO: test wrong api key
        
    def testRetrieveJson(self):
        actual = self.fixture.retrieveJson("/projectTag/12")
        self.assertTrue(actual, "retrieveJson: no result was returned.")
        #self.assertIn("Cranberry", actual, "Expecting a Cranberry project tag.")
        #Cranberry 
        
    def testGetEntity(self):
        actual = self.fixture.getEntity(12)
        self.assertTrue(actual, "No Project Tag returned.")
        self.assertIn("Cranberry", actual["name"],"Expecting Cranberry project tag.")
        
    def testGetIdsWithOffset(self):
        self.fixture.nameFilter="a"
        actualEntityIds, offset = self.fixture.getIdsWithOffset(0);
        self.assertEquals(2, len(actualEntityIds), "Expecting 10 ids, but got {}.".format(len(actualEntityIds)))
                
if __name__ == "__main__":
    unittest.main()
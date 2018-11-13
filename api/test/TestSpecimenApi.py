"""
Created on Mar 12, 2015

Note:
    These tests should be run against a local SeqDB instance. They do not
    mock the API. Since the work on SeqDB API is ongoing, they serve as both
    tests for the request methods as well as tests for the SeqDB API itself.
    
    These tests are therefore brittle. They not only depend on the current
    content of the local SeqDB instance, but also on the information to which
    the user has access (i.e. groups they are part of). What is meant by "user"
    is the SeqDB user whose API Key is used for testing (seqdb_api_key).

Before running the tests:
    - Ensure that config/config4tests.yaml is created from .sample and that
      this modified configuration file points to the local instance of SeqDB
    - Make sure that local instance of SeqDB is running

@author: Oksana Korol
"""
import requests
import unittest 
import yaml

from api.SpecimenApi import SpecimenApi
from config import config_root


class TestSpecimenApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(config_root.path() +
                  '/config4tests.yaml', 'r') as config_file:
            config = yaml.load(config_file)
            cls.fixture = SpecimenApi(
                api_key=config['seqdb']['api_key'],
                base_url=config['seqdb']['base_url'])
    
    def test_retrieve(self):
        # Test faulty connection
        self.assertRaises(requests.exceptions.ConnectionError,
                          self.fixture.retrieve, ')http://jibberish')

    def testGetEntity(self):
        actual = self.fixture.get_entity(6601)
        self.assertTrue(actual, 'No Specimen returned.')
        self.assertEqual(6601, actual['id'], 'Expecting specimen 6601.')
          
    def testGetIdsWithOffset(self):
        # TODO: fix this test
        # self.fixture.otherIds='|CV-F:CV547|'
        self.fixture.other_ids_filter = 'M.'
        actual_entity_ids = self.fixture.get_ids()
        self.assertEquals(2, len(actual_entity_ids),
                          'Expecting 10 ids, but got {}.'
                          .format(len(actual_entity_ids)))

                
if __name__ == '__main__':
    unittest.main()

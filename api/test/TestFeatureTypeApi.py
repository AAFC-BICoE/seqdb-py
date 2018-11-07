'''
Created on Apr 1, 2016

@author: korolo
'''
import requests
import unittest
import yaml

from config import config_root
from api.FeatureTypeApi import FeatureTypeApi


class TestFeatureTypeApi(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))
        self.fixture = FeatureTypeApi(api_key=config['seqdb']['api_key'],
                                                   base_url=config['seqdb']['base_url'])


    def setUp(self):
        pass


    def tearDown(self):
        pass

    def testGetFeatureTypes(self):
        # "http://localhost:2002/seqdb\/api/v1/featureType"
        actual = self.fixture.getFeatureTypesWithIds()
        self.assertTrue(actual, "No feature types returned.")
        self.assertIn("Quality Trim", actual, "No feature type 'Quality Trim' found.")


    def testCreateDeleteFeatureType(self):
        #curl -X POST -H "apikey:***REMOVED***" -H "Content-Type: application/json" -d '{"featureType":{"featureDescription":"test description 1231q","featureName":"test type 123123"}}' "***REMOVED***/featureType"

        #self.fixture.deleteFeatureType("19")

        featureTypeID = self.fixture.create("test type", "test type description")
        self.assertTrue(featureTypeID, "Feature type ID was not returned after feature type creation.")

        #curl -X DELETE -H "apikey: ***REMOVED***" "***REMOVED***/featureType/6"
        actual = self.fixture.delete(featureTypeID)
        self.assertEqual(200, actual['metadata']['statusCode'], "Could not delete feature type.")


    def testCreateDeleteFeatureType_multiple_duplicate(self):
        # select * from FeatureTypes where Name like "Test%";
        #self.fixture.delete("18")

        # Create
        featureTypeId1 = self.fixture.create("Test1", "Test")
        featureTypeId2 = self.fixture.create("Test2", "Test")
        self.assertTrue(featureTypeId2, "Second feature type ID was not returned after feature type creation.")
        """Test creation of a feature - expected to fail because of duplicate feature type id"""
        self.assertRaises(requests.exceptions.HTTPError, self.fixture.create, "Test1", "Duplicate of Test1")

        # Delete
        actual = self.fixture.delete(featureTypeId1)
        self.assertEqual(200, actual['metadata']['statusCode'], "Could not delete feature type.")
        actual = self.fixture.delete(featureTypeId2)
        self.assertEqual(200, actual['metadata']['statusCode'], "Could not delete feature type.")




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

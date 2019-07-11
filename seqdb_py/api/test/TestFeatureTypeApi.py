"""
Created on Apr 1, 2016

@author: korolo
"""
import requests
import unittest
import yaml

from config import config_root
from api.FeatureTypeApi import FeatureTypeApi


class TestFeatureTypeApi(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(config_root.path() +
                  '/config4tests.yaml', 'r') as config_file:
            config = yaml.safe_load(config_file)
            cls.fixture = FeatureTypeApi(
                api_key=config['seqdb']['api_key'],
                base_url=config['seqdb']['base_url'])

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_get_feature_types(self):
        # 'http:/***REMOVED***\/api/v1/featureType'
        actual = self.fixture.getFeatureTypesWithIds()
        self.assertTrue(actual, 'No feature types returned.')
        self.assertIn('Quality Trim', actual,
                      "No feature type 'Quality Trim' found.")

    def test_create_delete_feature_type(self):
        # curl -X POST -H 'apikey:***REMOVED***' -H
        # 'Content-Type: application/json' -d '{'featureType':
        # {'featureDescription':'test description 1231q',
        # 'featureName':'test type 123123'}}'
        # '***REMOVED***/api/v1/featureType'

        # self.fixture.deleteFeatureType('19')

        feature_type_id = self.fixture.create('test type',
                                              'test type description')
        self.assertTrue(feature_type_id,
                        'Feature type ID was not returned after '
                        'feature type creation.')

        # curl -X DELETE -H 'apikey: ***REMOVED***'
        # '***REMOVED***/api/v1/featureType/6'
        actual = self.fixture.delete(feature_type_id)
        self.assertEqual(200, actual['metadata']['statusCode'],
                         'Could not delete feature type.')

    def test_create_delete_feature_type_multiple_duplicate(self):
        # select * from FeatureTypes where Name like 'Test%';
        # self.fixture.delete('18')

        # Create
        feature_type_id1 = self.fixture.create('Test1', 'Test')
        feature_type_id2 = self.fixture.create('Test2', 'Test')
        self.assertTrue(feature_type_id2,
                        'Second feature type ID was not returned '
                        'after feature type creation.')
        '''Test creation of a feature - expected to fail because of duplicate feature type id'''
        self.assertRaises(requests.exceptions.HTTPError,
                          self.fixture.create, 'Test1', 'Duplicate of Test1')

        # Delete
        actual = self.fixture.delete(feature_type_id1)
        self.assertEqual(200, actual['metadata']['statusCode'],
                         'Could not delete feature type.')
        actual = self.fixture.delete(feature_type_id2)
        self.assertEqual(200, actual['metadata']['statusCode'],
                         'Could not delete feature type.')


if __name__ == '__main__':
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

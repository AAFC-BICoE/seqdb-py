'''
Created on Apr 1, 2016

@author: korolo
'''
import unittest

import yaml

from config import config_root
from api.FeatureApi import FeatureApi


class TestFeatureApi(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = FeatureApi(api_key=config['seqdb']['api_key'],
                                                   base_url=config['seqdb']['api_url'])
    

    def setUp(self):
        pass


    def tearDown(self):
        pass
    
    def testCreateGetDeleteFeature(self):
        
        # Create feature
        sampleFeatureLocations1 = [{"start":3,"end":5,"frame":1,"strand":1}]
        sampleFeatureLocations2 = [{"start":3,"end":5,"frame":1,"strand":1}, {"start":334,"end":454,"frame":2,"strand":1}]
        
        featureId = self.fixture.create("testName", 1, sampleFeatureLocations1, 1, "sample description", True)
        self.assertTrue(featureId, "Feature id was not returned after feature creation.")
        
        # Get feature
        retrieved_feature = self.fixture.getEntity(featureId)
        self.assertTrue(retrieved_feature, "No feature was retrieved.") 
        self.assertEqual("testName", retrieved_feature['name'], "Name of the retrieved feature does not match.")
        self.assertEqual("sample description", retrieved_feature['description'], "Feature description does not match")
        self.assertEqual(sampleFeatureLocations1, retrieved_feature['featureLocations'], "")
        self.assertEqual(1, retrieved_feature['featureType']['id'], "Feature type id does not match.")
        self.assertTrue(retrieved_feature['featureDefault'], "Feature default does not match")
        
        # Delete feature
        delete_jsn_resp = self.fixture.delete(featureId)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete feature type.")
        
        # Get feature again
        retrieved_feature = self.fixture.getEntity(featureId)
        self.assertFalse(retrieved_feature, "Unexpected: Feature was found after being deleted.") 
    
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
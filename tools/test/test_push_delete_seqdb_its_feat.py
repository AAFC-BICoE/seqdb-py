'''
Created on Feb 12, 2015

@author: korolo
'''
import unittest, os
from tools import push_to_seqdb, delete_seqdb_features
from api.seqdbWebService import seqdbWebService

test_url = '***REMOVED***:2002/seqdb/api/v1'
test_api_key = '***REMOVED***'

class Test(unittest.TestCase):


    def setUp(self):
        self.fixture = seqdbWebService(test_api_key, test_url)
        
    
    def tearDown(self):
        pass

        

    def testMain(self):    
        """ This test needs to be re-implemented
        ok_feature_ids = push_to_seqdb.push_its_features(test_api_key,"data/test.positions.txt", test_url )
        
        self.assertEqual(6, len(ok_feature_ids), "Expected 6 feature ids, but created %i feature ids." % len(ok_feature_ids))
        
        
        delete_seqdb_features.delete_features(test_api_key, push_to_seqdb.output_file_name, test_url)
         
        for feature_id in ok_feature_ids:
            self.assertFalse(self.fixture.getFeature(feature_id), "Feature was found after being delete.")   
        
        """ 
                
if __name__ == "__main__":
    unittest.main()
'''
Created on Feb 12, 2015

@author: korolo
'''
import unittest, os
from tools import push_seqdb_its_feat, delete_seqdb_features
from api.seqdbWebService import seqdbWebService

test_url = 'http://localhost:8181/seqdb.web-2.5/api/v1'
test_api_key = '***REMOVED***'

class Test(unittest.TestCase):


    def setUp(self):
        pass
        
    
    def tearDown(self):
        try:
            os.remove(push_seqdb_its_feat.output_file_name)
        except:
            pass

        

    def testMain(self):    
      
        ok_feature_ids, failed_feature_ids = push_seqdb_its_feat.main(test_api_key,"data/test.positions.txt", test_url )
        
        self.assertEqual(6, len(ok_feature_ids), "Expected 6 feature ids, but created %i feature ids." % len(ok_feature_ids))
        
        
        delete_seqdb_features.main(test_api_key, push_seqdb_its_feat.output_file_name, test_url)
         
        for feature_id in ok_feature_ids:
            self.assertFalse(self.seqdbWS.getFeature(feature_id), "Feature was found after being delete.")   
        
            
                
if __name__ == "__main__":
    unittest.main()
    
 
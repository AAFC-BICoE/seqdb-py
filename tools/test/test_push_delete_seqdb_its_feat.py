'''
Created on Feb 12, 2015

@author: korolo
'''
import unittest, os
from tools import push_seqdb_its_feat, delete_seqdb_features
from api.seqdbWebService import seqdbWebService

#test_url = 'http://10.117.203.138:8181/seqdb.web-2.5/api/v1'
test_url = '***REMOVED***:8181/seqdb.web-2.5/api/v1'
test_api_key = '***REMOVED***'

class Test(unittest.TestCase):


    def setUp(self):
        pass
        #local: "***REMOVED***:8181/seqdb.web-2.5/api/v1"
        #prod: "***REMOVED***/api/v1"
        
        self.seqdbWS = seqdbWebService(test_api_key, test_url)
    
    def tearDown(self):
        try:
            os.remove('blah.txt')
        except:
            pass
        try:
            os.remove(push_seqdb_its_feat.output_file_name)
        except:
            pass

    def testDelete(self):
        pass
        #self.seqdbWS.deleteFeature(85888)
        

    def testMain(self):    
      
        ok_feature_ids, failed_feature_ids = push_seqdb_its_feat.main(test_api_key,"data/test.positions.txt", test_url )
        
        self.assertEqual(6, len(ok_feature_ids), "Expected 6 feature ids, but created %i feature ids." % len(ok_feature_ids))
        
        
        delete_seqdb_features.main(test_api_key, push_seqdb_its_feat.output_file_name, test_url)
         
        for feature_id in ok_feature_ids:
            self.assertFalse(self.seqdbWS.getFeature(feature_id), "Feature was found after being delete.")   
        
            
                
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
    
 
'''
Created on Feb 12, 2015

@author: korolo
'''
import unittest, os
from tools import push_seqdb_its_feat, delete_seqdb_features
from api.seqdbWebService import seqdbWebService

#from seqdb_ws import json_seqdb_request, seqdb_ws_request

class Test(unittest.TestCase):


    def setUp(self):
        pass
        #local: "***REMOVED***:8181/seqdb.web-2.5/api/v1"
        #prod: "***REMOVED***/api/v1"
        
        self.seqdbWS = seqdbWebService('***REMOVED***', '***REMOVED***:8181/seqdb.web-2.5/api/v1')
    
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
        api_key = '***REMOVED***'
        base_url = '***REMOVED***:8181/seqdb.web-2.5/api/v1'  

        feature_ids = push_seqdb_its_feat.main(api_key,"data/test.positions.txt", base_url )
        
        self.assertEqual(6, len(feature_ids), "Expected 6 feature ids, but created %i feature ids." % len(feature_ids))
        
        
        delete_seqdb_features.main(api_key, push_seqdb_its_feat.output_file_name, base_url)
         
        for feature_id in feature_ids:
            self.assertFalse(self.seqdbWS.getFeature(feature_id), "Feature was found after being delete.")   
        

    
                
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
    
 
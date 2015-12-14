'''
Created on Feb 12, 2015

@author: korolo
'''
import unittest
import yaml

from tools import push_to_seqdb, delete_seqdb_features
from api.seqdbWebService import seqdbWebService
from config import config_root


class Test(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = seqdbWebService(api_key=config['seqdb_api_key'],
                                                   base_url=config['seqdb_api_url'])

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
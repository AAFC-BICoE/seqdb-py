'''
Created on Feb 12, 2015

@author: korolo
'''
import unittest
import yaml

from tools import push_to_seqdb, delete_seqdb_features
from api.seqdbWebService import seqdbWebService
from config import config_root
import os.path


class Test(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = seqdbWebService(api_key=config['seqdb']['api_key'],
                                                   base_url=config['seqdb']['api_url'])
        self.itsx_positions_file_name = "data/test.positions.txt"
        self.output_file_name = "seqdb_feature_ids.txt"
        self.failed_ids_output_file_name = "delete_failed_feature_ids.txt"
        
    
    @classmethod
    def tearDownClass(self):
        if os.path.isfile(self.output_file_name):
            os.remove(self.output_file_name)
        
        if os.path.isfile(self.failed_ids_output_file_name):
            os.remove(self.failed_ids_output_file_name) 
     
      
    def test_open_file(self):
        push_to_seqdb.open_file(self.itsx_positions_file_name)
        self.assertTrue(os.path.isfile(self.itsx_positions_file_name), "ITSx positions file was not opened")    


    def test_delete_from_seqdb(self):
        # time: 1.477s
        
        created_feature_ids = push_to_seqdb.push_its_features(self.fixture, self.itsx_positions_file_name)
        self.assertTrue(os.path.isfile(self.output_file_name), "SeqDB feature IDs file was not created")
        self.assertEqual(6, len(created_feature_ids), "Expected 6 feature IDs to be successfully written to SeqDB, but wrote %i feature IDs." % len(created_feature_ids))
        success_ids = delete_seqdb_features.delete_from_seqdb(self.fixture, self.output_file_name, "features")
        self.assertEqual(6, len(success_ids[0]), "Expected 6 feature IDs to be successfully deleted from SeqDB, but deleted %i feature IDs." % len(success_ids[0]))
        for created_feature_id in created_feature_ids:
            self.assertFalse(self.fixture.getFeature(created_feature_id), "Feature was found after being deleted.")

                
if __name__ == "__main__":
    unittest.main()
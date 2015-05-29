'''
Created on May 28, 2015

@author: korolo
'''
import unittest, os
from tools.seqdb_config_maker import SeqdbConfigMaker


class Test(unittest.TestCase):


    def setUp(self):
        self.config_file_name = "test_config.yaml"
        
    
    def tearDown(self):
        if os.path.isfile(self.config_file_name):
            os.remove(self.config_file_name)
        

    def testCreateConfigFile(self):
        config_file_abs = SeqdbConfigMaker(api_url="http://***REMOVED***.ca/seqdb", config_file_name=self.config_file_name).createConfigFile(api_key="blah")
        
        self.assertTrue(config_file_abs, "Config file was not created.")
        self.assertTrue(self.config_file_name in config_file_abs, "Created config file is not named as expected.")
        
        config_file_abs = open(self.config_file_name, 'r') 
        actual = config_file_abs.read()
        config_file_abs.close()
        
        actual_lines = actual.splitlines()
        expected_lines = ['seqdb:','    api_key: "blah"','    api_url: "http://***REMOVED***.ca/seqdb"']

        self.assertEqual(expected_lines, actual_lines, "User config file does not match expected content.")
                
if __name__ == "__main__":
    unittest.main()
'''
Created on Oct 16, 2015

@author: korolo
'''
import unittest
import yaml

from tools import pull_seqdb_seqs
from api.seqdbWebService import seqdbWebService
from config import config_root


class Test(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = seqdbWebService(api_key=config['seqdb_api_key'],
                                                   base_url=config['seqdb_api_url'])
    
    '''    
    def test_get_ITS_seq_ids(self): 
        # time: 1117.424s    
        seq_ids = pull_seqdb_seqs.get_ITS_seq_ids(self.fixture)
        self.assertEqual(17787, len(seq_ids), "Expected 17787 ITS sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    '''
    """        
    def test_get_seq_ids_all(self):    
        # Note: this test takes REALLY long time, like 40 hours long time
        seq_ids = pull_seqdb_seqs.get_seq_ids(fixture=self.fixture, pull_type="all")
        self.assertEqual(4782554 , len(seq_ids), "Expected 4782554 sequences (total), but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    """
        
    def test_get_seq_ids_consensus(self):
        # time: 8.324s 
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus")
        self.assertEqual(5555 , len(seq_ids), "Expected 5555 consensus sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
        
    def test_get_seq_ids_specimen(self):  
        # time: 9.248s
        
        # consensus  
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, "consensus", [4405])  
        self.assertEqual(1 , len(seq_ids), "Expected 1 consensus sequences, but got {}. ".format(len(seq_ids)))

        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, "consensus", [4405,4264])  
        self.assertTrue(seq_ids, "No Sequence ids returned.")
        self.assertEqual(3 , len(seq_ids), "Expected 3 consensus sequences, but got {}. ".format(len(seq_ids)))
        self.assertIn(358301, seq_ids, "Sequence id 358301 is expected to be associated with specimen 4405.")
        self.assertIn(358302, seq_ids, "Sequence id 358302 is expected to be associated with specimen 4264.")
   
        # raw        
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, "raw", [4405])  
        self.assertEqual(22 , len(seq_ids), "Expected 22 sequences, but got {}. ".format(len(seq_ids)))
        
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, "raw", [4405,4264])  
        self.assertEqual(33 , len(seq_ids), "Expected 33 sequences, but got {}. ".format(len(seq_ids)))
        self.assertIn(27755, seq_ids, "Sequence id 27755 is expected to be associated with specimen 4405.")
        self.assertIn(28262, seq_ids, "Sequence id 358302 is expected to be associated with specimen 4264.")
      

    def test_get_seq_ids_sequenceName(self):  
        # time: 0.182s  
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", sequence_name="Pyt_arrhenomanes_")  
        self.assertEqual(5, len(seq_ids), "Expected 5 sequences, but got %i. " % len(seq_ids))
           
    def test_get_ids_sampleName(self):
        #time:
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", sample_name="LEV4183")
        self.assertEquals(1, len(seq_ids), "Expected 1 sequence, but got %i. " % len(seq_ids))
        
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", sample_name="LEV6103")
        self.assertEquals(60, len(seq_ids), "Expected 60 sequences, but got %i. " % len(seq_ids))   

    def test_get_seq_ids_regionName(self): 
        # time: 23.918s
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", region_name="28s")  
        self.assertEqual(69, len(seq_ids), "Expected 69 sequences, but got %i. " % len(seq_ids))
 
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", region_name="ef-1a")  
        self.assertEqual(492, len(seq_ids), "Expected 492 sequences, but got %i. " % len(seq_ids))

    def test_get_seq_ids_projectName(self):  
        # time: 36.331s  
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", project_name="grdi")  
        self.assertEqual(5555, len(seq_ids), "Expected 5555 sequences, but got %i. " % len(seq_ids))

        # TODO: takes too long
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", project_name="Pythium Type Specimens")  
        self.assertEqual(4331, len(seq_ids), "Expected 4331 sequences, but got %i. " % len(seq_ids))

    def test_get_seq_ids_colletionCode(self): 
        # time: 1.35s
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", collection_code="lev")  
        self.assertEqual(211, len(seq_ids), "Expected 211 sequences, but got %i. " % len(seq_ids))
  
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", collection_code="pm")  
        self.assertEqual(148, len(seq_ids), "Expected 32936 sequences, but got %i. " % len(seq_ids))


    def test_get_seq_ids_taxonomy(self):    
        # time: 15.135s
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", taxonomy_rank="species", taxonomy_value="megasperma")  
        self.assertEqual(3, len(seq_ids), "Expected 3 sequences, but got %i. " % len(seq_ids))

        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", taxonomy_rank="species", taxonomy_value="megasperma")  
        self.assertEqual(215, len(seq_ids), "Expected 215 sequences, but got %i. " % len(seq_ids))

                
if __name__ == "__main__":
    unittest.main()
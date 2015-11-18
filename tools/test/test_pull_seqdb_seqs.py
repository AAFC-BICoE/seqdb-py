'''
Created on Oct 16, 2015

@author: korolo
'''
import unittest
from tools import pull_seqdb_seqs
from api.seqdbWebService import seqdbWebService

test_url = '***REMOVED***:2002/seqdb/api/v1/'
test_api_key = '***REMOVED***'

class Test(unittest.TestCase):


    def setUp(self):
        self.fixture = seqdbWebService(test_api_key, test_url)
        
    
    def tearDown(self):
        pass
        
    """
    def test_get_ITS_seq_ids(self): 
        # Takes a few minutes to run   
        seq_ids = pull_seqdb_seqs.get_ITS_seq_ids(self.fixture)
        self.assertEqual(18070, len(seq_ids), "Expected 18070 ITS sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    """    
    """        
    def test_get_seq_ids_all(self):    
        # Note: this test takes REALLY long time, like 40 hours long time
        seq_ids = pull_seqdb_seqs.get_seq_ids(fixture=self.fixture, pull_type="all")
        self.assertEqual(4782554 , len(seq_ids), "Expected 4782554 sequences (total), but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    """
        
    def test_get_seq_ids_consensus(self):
        # time: 4.493s 
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus")
        self.assertEqual(3047 , len(seq_ids), "Expected 3047 consensus sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
        
    def test_get_seq_ids_specimen(self):  
        # time: 3.131s  
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, "consensus", 4405)  
        self.assertEqual(1 , len(seq_ids), "Expected 1 consensus sequences, but got %i. " % len(seq_ids))

        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, "all", 4405)  
        self.assertEqual(23 , len(seq_ids), "Expected 23 sequences, but got %i. " % len(seq_ids))

    def test_get_seq_ids_sequenceName(self):  
        # time: 1.106s  
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", sequence_name="Pyt_arrhenomanes_")  
        self.assertEqual(5, len(seq_ids), "Expected 5 sequences, but got %i. " % len(seq_ids))

    def test_get_seq_ids_regionName(self): 
        # time: 0.091s   
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", region_name="28s")  
        self.assertEqual(25, len(seq_ids), "Expected 25 sequences, but got %i. " % len(seq_ids))

        # time: 0.091s   
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="all", region_name="ef-1a")  
        self.assertEqual(492, len(seq_ids), "Expected 492 sequences, but got %i. " % len(seq_ids))

    def test_get_seq_ids_projectName(self):  
        # time: 4.696s  
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", project_name="grdi")  
        self.assertEqual(3047, len(seq_ids), "Expected 3047 sequences, but got %i. " % len(seq_ids))

        # TODO: takes too long
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="all", project_name="Pythium Type Specimens")  
        self.assertEqual(4373, len(seq_ids), "Expected 4373 sequences, but got %i. " % len(seq_ids))

    def test_get_seq_ids_colletionCode(self): 
        # time: 0.439s
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", collection_code="lev")  
        self.assertEqual(211, len(seq_ids), "Expected 211 sequences, but got %i. " % len(seq_ids))

        # time: 1.810s   
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="all", collection_code="pm")  
        self.assertEqual(148, len(seq_ids), "Expected 32936 sequences, but got %i. " % len(seq_ids))


    def test_get_seq_ids_taxonomy(self):    
        # time: 0.664s
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", taxonomy_rank="species", taxonomy_value="megasperma")  
        self.assertEqual(3, len(seq_ids), "Expected 3 sequences, but got %i. " % len(seq_ids))

        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="all", taxonomy_rank="species", taxonomy_value="megasperma")  
        self.assertEqual(218, len(seq_ids), "Expected 218 sequences, but got %i. " % len(seq_ids))

                
if __name__ == "__main__":
    unittest.main()
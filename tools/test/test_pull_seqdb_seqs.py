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
        self.seqdbWS = seqdbWebService(test_api_key, test_url)
        
    
    def tearDown(self):
        pass
        

    def test_get_ITS_seq_ids(self):    
        seq_ids = pull_seqdb_seqs.get_ITS_seq_ids(self.seqdbWS)
        self.assertEqual(18070, len(seq_ids), "Expected 18070 ITS sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
        
    def test_get_consensus_seq_ids(self):    
        # Note that this function is not used for some reason. Evaluate if this is the way to go.
        seq_ids = pull_seqdb_seqs.get_consensus_seq_ids(self.seqdbWS)        
        self.assertEqual(566, len(seq_ids), "Expected 566 consensus sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))

    """        
    def test_get_seq_ids_all(self):    
        # Note: this test takes REALLY long time, like 40 hours long time
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.seqdbWS, pull_type="all")
        self.assertEqual(4782554 , len(seq_ids), "Expected 4782554 sequences (total), but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    """
        
    def test_get_seq_ids_consensus(self):    
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.seqdbWS, pull_type="consensus")
        self.assertEqual(566 , len(seq_ids), "Expected 566 consensus sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
        
            
                
if __name__ == "__main__":
    unittest.main()
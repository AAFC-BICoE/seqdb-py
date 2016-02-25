'''
Created on Oct 16, 2015

@author: korolo
'''
import unittest, os
import yaml

from tools import pull_seqdb_seqs
from api.seqdbWebService import seqdbWebService
from config import config_root
from tools.pull_seqdb_seqs import return_types


class Test(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = seqdbWebService(api_key=config['seqdb_api_key'],
                                                   base_url=config['seqdb_api_url'])
        # Create a temporary file
        self.file_name = "test_seqdb_sequences."
        
        if return_types:
            self.file_type = "fasta"
        else:
            self.file_type = "fastq"
        
    def tearDown(self):
        # Remove the file after the test
        if os.path.isfile(self.file_name + self.file_type):
            os.remove(self.file_name + self.file_type)
            
    '''    
    def test_get_ITS_seq_ids(self): 
        # time: 1117.424s    
        seq_ids = pull_seqdb_seqs.get_ITS_seq_ids(self.fixture)
        self.assertEqual(17787, len(seq_ids), "Expected 17787 ITS sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    '''
    
    """       
    def test_get_seq_ids_all(self):
        # time: 9666.604s    
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="all")
        self.assertEqual(485643 , len(seq_ids), "Expected 485,643 sequences (total), but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    """
        
    def test_get_seq_ids_consensus(self):
        # time: 15.324s 
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus")
        self.assertEqual(5555 , len(seq_ids), "Expected 5,555 consensus sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    
    def test_get_seq_ids_raw(self):
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw")
        self.assertEquals(480088, len(seq_ids), "Expected 480,088 raw sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
             
    def test_get_seq_ids_specimen(self):  
        # time: 14.322s
        
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
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, "all", [4405])
        self.assertEqual(23, len(seq_ids), "Expected 23 sequences, but got {}.".format(len(seq_ids)))

    def test_get_seq_ids_sequenceName(self):  
        # time: 11.840s
        
        # consensus  
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", sequence_name="Pyt_arrhenomanes_")  
        self.assertEqual(5, len(seq_ids), "Expected 5 sequences, but got %i. " % len(seq_ids))
        
        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, pull_type="raw" , sequence_name="Y10_16_F167")
        self.assertEqual(2, len(seq_ids), "Expected 2 sequences, but got %i" % len(seq_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, pull_type="all", sequence_name="Y10_16_F167")
        self.assertEqual(2, len(seq_ids), "Expected 2 sequences, but got %i" % len(seq_ids))
        
    def test_get_ids_sampleName(self):
        #time: 12.737s
        
        # consensus
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", sample_name="LEV4183")
        self.assertEquals(1, len(seq_ids), "Expected 1 sequence, but got %i. " % len(seq_ids))
        
        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", sample_name="LEV6103")
        self.assertEquals(60, len(seq_ids), "Expected 60 sequences, but got %i. " % len(seq_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="all", sample_name="LEV6103")
        self.assertEquals(61, len(seq_ids), "Expected 61 sequences, but got %i. " % len(seq_ids))

    def test_get_seq_ids_regionName(self): 
        # time: 78.108s
        
        # consensus
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", region_name="28s")  
        self.assertEqual(69, len(seq_ids), "Expected 69 sequences, but got %i. " % len(seq_ids))
 
        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", region_name="ef-1a")  
        self.assertEqual(492, len(seq_ids), "Expected 492 sequences, but got %i. " % len(seq_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="all", region_name="ACA")  
        self.assertEqual(1041, len(seq_ids), "Expected 1,041 sequences, but got %i. " % len(seq_ids))

    def test_get_seq_ids_projectName(self):  
        # time: 94.259s 
        
        # consensus 
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", project_name="grdi")  
        self.assertEqual(5555, len(seq_ids), "Expected 5,555 sequences, but got %i. " % len(seq_ids))

        # raw
        # TODO: takes too long
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", project_name="Pythium Type Specimens")  
        self.assertEqual(4331, len(seq_ids), "Expected 4,331 sequences, but got %i. " % len(seq_ids))
        
        # all 
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, pull_type="all", project_name="Pythium Type Specimens")
        self.assertEqual(4373, len(seq_ids), "Expected 4,373 sequences, but got %i." % len(seq_ids))

    def test_get_seq_ids_colletionCode(self): 
        # time: 2.164s
        
        # consensus
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", collection_code="lev")  
        self.assertEqual(211, len(seq_ids), "Expected 211 sequences, but got %i. " % len(seq_ids))
  
        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", collection_code="pm")  
        self.assertEqual(148, len(seq_ids), "Expected 148 sequences, but got %i. " % len(seq_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, pull_type="all", collection_code="BISS")
        self.assertEqual(25, len(seq_ids), "Expected 25 sequences, but got %i." % len(seq_ids))

    def test_get_seq_ids_taxonomy(self):    
        # time: 45.438s
        
        # consensus
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", taxonomy_rank="species", taxonomy_value="megasperma")  
        self.assertEqual(3, len(seq_ids), "Expected 3 sequences, but got %i. " % len(seq_ids))

        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", taxonomy_rank="species", taxonomy_value="megasperma")  
        self.assertEqual(215, len(seq_ids), "Expected 215 sequences, but got %i. " % len(seq_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, pull_type="all", taxonomy_rank="species", taxonomy_value="megasperma")
        self.assertEqual(218, len(seq_ids), "Expected 218 sequences, but got %i." % len(seq_ids))

           
    def test_write_sequence_file(self):
        # time 14.067s
        # consensus
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", sample_name="LEV4183")
        success_ids = pull_seqdb_seqs.write_sequence_file(self.fixture, seq_ids, file_name=self.file_name, file_type="fasta")
        self.assertEqual(1, len(success_ids), "The output file is created. The file is expected to contain 1 sequence, but contains %i." % len(success_ids))
        
        """
        output_file = open(self.file_name + "fasta", 'r')
        actual = output_file.readline()
        output_file.close()
        print actual
        expected_first_line = ">seqdb|358455 Pythium scleroteichum Phy_operculata_CBS24183_ACA \n"
        self.assertEqual(expected_first_line, actual, "File does not match expected content.")
        """
        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", sample_name="LEV6103")
        success_ids = pull_seqdb_seqs.write_sequence_file(self.fixture, seq_ids, file_name=self.file_name, file_type="fastq")
        self.assertEqual(60, len(success_ids), "The output file is created. The file is expected to contain 60 sequence, but contains %i." % len(success_ids))
    
        """
        output_file = open(self.file_name + "fastq", 'r')
        actual = output_file.readline()
        output_file.close()
        print actual
        expected_first_line = "@seqdb|266400 Myzocytiopsis ? sp. affin. intermedia AL_HM_H047_09_SEQ_LEV6103_18S_NS1\n"
        self.assertEqual(expected_first_line, actual, "File does not match expected content.")        
        """
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="all", sample_name="LEV6103")
        success_ids = pull_seqdb_seqs.write_sequence_file(self.fixture, seq_ids, file_name=self.file_name, file_type="fasta")
        self.assertEqual(61, len(success_ids), "The output file was created. It is expected to contain 61 sequence, but contains %i." % len(success_ids))
        
        
        output_file = open(self.file_name + "fasta", 'r')
        actual = output_file.readline()
        output_file.close()
        print actual
        expected_first_line = ">seqdb|266400 Myzocytiopsis ? sp. affin. intermedia AL_HM_H047_09_SEQ_LEV6103_18S_NS1 \n"
        self.assertEqual(expected_first_line, actual, "File does not match expected content.")        
        
 
if __name__ == "__main__":
    unittest.main()
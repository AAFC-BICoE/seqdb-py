'''
Created on Oct 16, 2015

@author: korolo
'''
import os.path
import unittest

import yaml

from api.seqdbWebService import seqdbWebService
from config import config_root
from tools import pull_seqdb_seqs


class TestPullSeqdbSeqs(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = seqdbWebService(api_key=config['seqdb']['api_key'],
                                                   base_url=config['seqdb']['api_url'])
        self.output_file_name = "test_output_file."
        self.output_fasta_file_name = self.output_file_name + "fasta"
        self.output_fastq_file_name = self.output_file_name + "fastq"
        self.output_taxon_file_name = "test_output_taxonomy.txt"    

    @classmethod
    def tearDownClass(self):
        if os.path.isfile(self.output_taxon_file_name):
            os.remove(self.output_taxon_file_name) 
        
        if os.path.isfile(self.output_fasta_file_name):
            os.remove(self.output_fasta_file_name) 
        
        if os.path.isfile(self.output_fastq_file_name):
            os.remove(self.output_fastq_file_name) 
        
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
        # time: 8.395s 
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus")
        self.assertEqual(5555 , len(seq_ids), "Expected 5,555 consensus sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    
    """
    def test_get_seq_ids_raw(self):
        # time: stopped
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw")
        self.assertEquals(480088, len(seq_ids), "Expected 480,088 raw sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    """    
    def test_execute_script_consensus(self):
        # time: 127.712s
        
        # Testing fasta file creation for consensus sequences 
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "consensus"], 
                             self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        self.assertEqual(13983492, os.stat(self.output_fasta_file_name).st_size, "File size expected to be 13983492 bytes, but is %i" %os.stat(self.output_fasta_file_name).st_size)
        output_file = open(self.output_fasta_file_name)
        actual_lines = output_file.readlines()
        expected_line_0 = '>seqdb|358301 Pythium arrhenomanes Pyt_arrhenomanes_BR1028_ACA \n'
        self.assertEqual(expected_line_0, actual_lines[0], "Expected line is not the same as the actual line.")
        expected_line_89597 = '>seqdb|4823203 Couesius plumbeus NXG2013358A_RHO \n'
        self.assertEqual(expected_line_89597, actual_lines[89597], "Expected line is not the same as the actual line.")
        expected_line_161610 = '>seqdb|4829279   12G4130_GRSPaV_partial \n'
        self.assertEqual(expected_line_161610, actual_lines[161610], "Expected line is not the same as the actual line.")
        
    def test_execute_script_raw(self):
        # time: 28.616
        
        # Testing fasta file creation for raw sequences 
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "raw", "--seqName", "S-SH-"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        self.assertEqual(138996, os.stat(self.output_fasta_file_name).st_size, "File size expected to be 138996 bytes, but is %i" %os.stat(self.output_fasta_file_name).st_size)
        output_file = open(self.output_fasta_file_name)
        actual_lines = output_file.readlines()
        expected_line_0 = '>seqdb|1 G07_S-SH-27_R25_its4_3.ab1 \n'
        self.assertEqual(expected_line_0, actual_lines[0], "Expected line is not the same as the actual line.")
        expected_line_1026 = '>seqdb|79390 D08_S-SH-98_RS13_LROR_4.ab1 \n'
        self.assertEqual(expected_line_1026, actual_lines[1026], "Expected line is not the same as the actual line.")
        expected_line_1838 = '>seqdb|126059 B01_S-SH-36_M2_its4a_2.ab1 \n'
        self.assertEqual(expected_line_1838, actual_lines[1838], "Expected line is not the same as the actual line.")
        
        # Testing fastq file creation for raw sequences 
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fastq", "raw", "--sampleName", "LEV6103"], 
                             self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fastq_file_name), "Fastq file was not created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        self.assertEqual(106623, os.stat(self.output_fastq_file_name).st_size, "File size expected to be 106623 bytes, but is %i" %os.stat(self.output_fastq_file_name).st_size)
        output_file = open(self.output_fastq_file_name)
        actual_lines = output_file.readlines()
        expected_line_0 = '@seqdb|266400 Myzocytiopsis ? sp. affin. intermedia AL_HM_H047_09_SEQ_LEV6103_18S_NS1\n'
        self.assertEqual(expected_line_0, actual_lines[0], "Expected line is not the same as the actual line.")
        expected_line_140 = '@seqdb|301609 Myzocytiopsis ? sp. affin. intermedia AL_TR_T126_41_SEQ_LEV6103A_29_TF_M13_LO\n'
        self.assertEqual(expected_line_140, actual_lines[140], "Expected line is not the same as the actual line.")
        expected_line_228 = '@seqdb|331086 Myzocytiopsis ? sp. affin. intermedia AL_HM_H080_25_SEQ_LEV6103_BTUB_OOM-BTUB-LO-1401\n'
        self.assertEqual(expected_line_228, actual_lines[228], "Expected line is not the same as the actual line.")
        
    def test_execute_script_all(self):
        # time: 72.028s
        
        # Testing fasta file creation for all sequences
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "all", "--geneRegion", "ACA"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        output_file = open(self.output_fasta_file_name)
        actual_lines = output_file.readlines()
        expected_line_0 = '>seqdb|6818 Pythium irregulare G079_08_BR1000_ACA_ACAUP16 \n'
        self.assertEqual(expected_line_0, actual_lines[0], "Expected line is not the same as the actual line.")
        expected_line_7794 = '>seqdb|112673 Pythium dimorphum T041_014_M0236_ACA_ACAUP829PP \n'
        self.assertEqual(expected_line_7794, actual_lines[7794], "Expected line is not the same as the actual line.")
        expected_line_14604 = '>seqdb|358505 Phytophthora porri Phy_porri_Lev1964_ACA \n'
        self.assertEqual(expected_line_14604, actual_lines[14604], "Expected line is not the same as the actual line.")
        
             
    def test_get_seq_ids_specimen(self):  
        # time: 13.833s
        
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
        # time: 8.732s
        
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
        #time: 9.368s
        
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
        # time: 61.073s
        
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
        # time: 63.894s 
        
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
        # time: 1.508s
        
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
        # time: 31.726s
        
        # consensus
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", taxonomy_rank="species", taxonomy_value="megasperma")  
        self.assertEqual(3, len(seq_ids), "Expected 3 sequences, but got %i. " % len(seq_ids))

        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", taxonomy_rank="species", taxonomy_value="megasperma")  
        self.assertEqual(215, len(seq_ids), "Expected 215 sequences, but got %i. " % len(seq_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.fixture, pull_type="all", taxonomy_rank="species", taxonomy_value="megasperma")
        self.assertEqual(218, len(seq_ids), "Expected 218 sequences, but got %i." % len(seq_ids))
     
    """
    def test_write_sequence_file(self):
        # time 19.587ss
        
        # consensus
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus", sample_name="LEV4183")
        success_ids = pull_seqdb_seqs.write_sequence_file(self.fixture, seq_ids, file_name="test_1_seqdb_sequences.", file_type="fasta")
        self.assertEqual(1, len(success_ids), "The output file is created. The file is expected to contain 1 sequence, but contains %i." % len(success_ids))
        
        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", sample_name="LEV6103")
        success_ids = pull_seqdb_seqs.write_sequence_file(self.fixture, seq_ids, file_name="test_2_seqdb_sequences.", file_type="fastq")
        self.assertEqual(60, len(success_ids), "The output file is created. The file is expected to contain 60 sequence, but contains %i." % len(success_ids))
        
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw", sample_name="LEV6103")
        success_ids = pull_seqdb_seqs.write_sequence_file(self.fixture, seq_ids, "test_3_seqdb_sequences.", file_type="fasta")
        self.assertEqual(60, len(success_ids), "The output file is created. The file is expected to contain 60 sequence, but contains %i." % len(success_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="all", sample_name="LEV6103")
        success_ids = pull_seqdb_seqs.write_sequence_file(self.fixture, seq_ids, file_name="test_4_seqdb_sequences.", file_type="fasta")
        self.assertEqual(61, len(success_ids), "The output file was created. It is expected to contain 61 sequence, but contains %i." % len(success_ids))
    """
 
if __name__ == "__main__":
    unittest.main()
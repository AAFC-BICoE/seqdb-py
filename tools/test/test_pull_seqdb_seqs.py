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
                                                   base_url=config['seqdb']['base_url'])
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
    #Method no longer used
    
    def test_get_ITS_seq_ids(self): 
        #OK. Time: 237.375s     
        seq_ids = pull_seqdb_seqs.get_ITS_seq_ids(self.rawSeqFixture)
        self.assertEqual(17787, len(seq_ids), "Expected 17787 ITS sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    
    def test_get_seq_ids_consensus(self):
        #OK. Time 6.462s 
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, pull_type="consensus")
        self.assertEqual(5555 , len(seq_ids), "Expected 5,555 consensus sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
           
    def test_get_seq_ids_all(self):
        #OK. Time: 1926.634s    
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, pull_type="all")
        self.assertEqual(485643 , len(seq_ids), "Expected 485,643 sequences (total), but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))

    def test_get_seq_ids_raw(self):
        #OK. Time: 1885.130s
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, pull_type="raw")    
        self.assertEqual(480088 , len(seq_ids), "Expected 480,088 sequences (total), but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    '''
            
    '''
    #Method no longer used
    
    ###TESTING FILTERS
    
    ### pull_seqdb_seqs.get_seq_ids is removed. TODO: make sure that the equivalent tests exist in api tests for filters

    def test_get_seq_ids_specimen(self):  
        # OK. Time: 4.761s
        
        # consensus  
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "consensus", [4405])  
        self.assertEqual(1 , len(seq_ids), "Expected 1 consensus sequences, but got {}.".format(len(seq_ids)))

        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "consensus", [4405,4264])  
        self.assertTrue(seq_ids, "No Sequence ids returned.")
        self.assertEqual(3 , len(seq_ids), "Expected 3 consensus sequences, but got {}.".format(len(seq_ids)))
        self.assertIn(358301, seq_ids, "Sequence id 358301 is expected to be associated with specimen 4405.")
        self.assertIn(358302, seq_ids, "Sequence id 358302 is expected to be associated with specimen 4264.")
   
        # raw        
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "raw", [4405])  
        self.assertEqual(22 , len(seq_ids), "Expected 22 sequences, but got {}.".format(len(seq_ids)))
        
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "raw", [4405,4264])  
        self.assertEqual(33 , len(seq_ids), "Expected 33 sequences, but got {}.".format(len(seq_ids)))
        self.assertIn(27755, seq_ids, "Sequence id 27755 is expected to be associated with specimen 4405.")
        self.assertIn(28262, seq_ids, "Sequence id 358302 is expected to be associated with specimen 4264.")
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", [4405])
        self.assertEqual(23, len(seq_ids), "Expected 23 sequences, but got {}.".format(len(seq_ids)))

    def test_get_seq_ids_sequenceName(self):  
        # OK. Time: 6.107s
        
        # consensus  
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "consensus", sequence_name="Pyt_arrhenomanes_")  
        self.assertEqual(5, len(seq_ids), "Expected 5 sequences, but got %i." % len(seq_ids))
        
        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "raw", sequence_name="Y10_16_F167")
        self.assertEqual(2, len(seq_ids), "Expected 2 sequences, but got %i" % len(seq_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", sequence_name="Y10_16_F167")
        self.assertEqual(2, len(seq_ids), "Expected 2 sequences, but got %i" % len(seq_ids))
        
    def test_get_ids_sampleName(self):
        #OK. Time 3.663s
        
        # consensus
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "consensus", sample_name="LEV4183")
        self.assertEquals(1, len(seq_ids), "Expected 1 sequence, but got %i." % len(seq_ids))
        
        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "raw", sample_name="LEV6103")
        self.assertEquals(60, len(seq_ids), "Expected 60 sequences, but got %i." % len(seq_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", sample_name="LEV6103")
        self.assertEquals(61, len(seq_ids), "Expected 61 sequences, but got %i." % len(seq_ids))

    def test_get_seq_ids_regionName(self): 
        #OK. Time: 13.739s
        
        # consensus
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "consensus", region_name="28s")  
        self.assertEqual(69, len(seq_ids), "Expected 69 sequences, but got %i." % len(seq_ids))
 
        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "raw", region_name="ef-1a")  
        self.assertEqual(492, len(seq_ids), "Expected 492 sequences, but got %i." % len(seq_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", region_name="ACA")  
        self.assertEqual(1041, len(seq_ids), "Expected 1,041 sequences, but got %i." % len(seq_ids))

    
    #Note: Failing. Error message: exceptions must be old-style classes or derived from BaseException, not str
    def test_get_seq_ids_projectName(self):  
        # time: 63.894s 
        
        # consensus 
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "consensus", project_name="grdi")  
        self.assertEqual(5555, len(seq_ids), "Expected 5,555 sequences, but got %i." % len(seq_ids))

        # raw
        # TODO: takes too long
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "raw", project_name="Pythium Type Specimens")  
        self.assertEqual(4331, len(seq_ids), "Expected 4,331 sequences, but got %i." % len(seq_ids))
        
        # all 
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", project_name="Pythium Type Specimens")
        self.assertEqual(4373, len(seq_ids), "Expected 4,373 sequences, but got %i." % len(seq_ids))
    
        
    def test_get_seq_ids_colletionCode(self): 
        #OK. Time: 1.042s
        
        # consensus
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "consensus", collection_code="lev")  
        self.assertEqual(211, len(seq_ids), "Expected 211 sequences, but got %i." % len(seq_ids))
  
        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "raw", collection_code="pm")  
        self.assertEqual(148, len(seq_ids), "Expected 148 sequences, but got %i." % len(seq_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", collection_code="BISS")
        self.assertEqual(25, len(seq_ids), "Expected 25 sequences, but got %i." % len(seq_ids))

    def test_get_seq_ids_taxonomy(self):    
        #OK. Time: 11.161s
        
        # consensus
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "consensus", taxonomy_rank="species", taxonomy_value="megasperma")  
        self.assertEqual(3, len(seq_ids), "Expected 3 sequences, but got %i." % len(seq_ids))

        # raw
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "raw", taxonomy_rank="species", taxonomy_value="megasperma")  
        self.assertEqual(215, len(seq_ids), "Expected 215 sequences, but got %i." % len(seq_ids))
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", taxonomy_rank="species", taxonomy_value="megasperma")
        self.assertEqual(218, len(seq_ids), "Expected 218 sequences, but got %i." % len(seq_ids))
    '''


    ### TESTING FASTA FILE CREATION 

    def test_execute_script_consensus_fasta(self):
        #OK. Time: 50.61s
        
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "consensus"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        #Better test for the number of sequences in the file
        #self.assertEqual(13983492, os.stat(self.output_fasta_file_name).st_size, "File size expected to be 13983492 bytes, but is %i." %os.stat(self.output_fasta_file_name).st_size)
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    idList.append(line.split()[0])
        self.assertIn('>seqdb|358301', idList, "Expected sequence ID 358301 is not found in the file")
        self.assertIn('>seqdb|4823203', idList, "Expected sequence ID 4823203 is not found in the file")
        self.assertIn('>seqdb|4829279', idList, "Expected sequence ID 4829279 is not found in the file")
        
        
    def test_execute_script_raw_fasta(self):
        #OK. Time: 26.66s
        
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "raw", "--seqName", "S-SH-"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        #Better test for the number of sequences in the file
        #self.assertEqual(138996, os.stat(self.output_fasta_file_name).st_size, "File size expected to be 138996 bytes, but is %i." %os.stat(self.output_fasta_file_name).st_size)
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    idList.append(line.split()[0])        
        self.assertIn('>seqdb|1', idList, "Expected sequence ID 1 is not found in the file")
        self.assertIn('>seqdb|79390', idList, "Expected sequence ID 79390 is not found in the file")
        self.assertIn('>seqdb|126059', idList, "Expected sequence ID 126059 is not found in the file")
    
 
    def test_execute_script_all_fasta(self):
        #OK. Time: 28.818s
        
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "all", "--geneRegion", "EF-1a"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        #Better test for the number of sequences in the file
        #self.assertEqual(557033, os.stat(self.output_fasta_file_name).st_size, "File size expected to 557033 bytes, but is %i." %os.stat(self.output_fasta_file_name).st_size)
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    idList.append(line.split()[0])        
        self.assertIn('>seqdb|1689', idList, "Expected sequence ID 1689 is not found in the file")
        self.assertIn('>seqdb|103372', idList, "Expected sequence ID 103372 is not found in the file")
        self.assertIn('>seqdb|149807', idList, "Expected sequence ID 149807 is not found in the file")
        
        
    ### TESTING FASTQ FILE CREATION
        
    def test_execute_script_raw_fastq(self):
        #OK. Time: 5.16s
        
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fastq", "raw", "--sampleName", "LEV6103"], 
                             self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fastq_file_name), "Fastq file was not created.")
        #Better test for the number of sequences in the file
        #self.assertEqual(106623, os.stat(self.output_fastq_file_name).st_size, "File size expected to be 106623 bytes, but is %i." %os.stat(self.output_fastq_file_name).st_size)
        idList = []
        with open(self.output_fastq_file_name) as f:
            for line in f:
                if line.startswith('@'):
                    idList.append(line.split()[0])
        self.assertIn('@seqdb|266400', idList, "Expected sequence ID 266400 is not found in the file")
        self.assertIn('@seqdb|301609', idList, "Expected sequence ID 301609 is not found in the file")
        self.assertIn('@seqdb|331086', idList, "Expected sequence ID 331086 is not found in the file")   
 
     
    ### TESTING TAXONOMY FILE CREATION
    
    def test_execute_script_consensus_taxonomy(self):    
        #OK. Time: 1.06s

        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "-t", "consensus", "--seqName", "Pyt_arrhenomanes_"], 
                                self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertTrue(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        #Better test for the number of sequences in the file
        #self.assertEqual(515, os.stat(self.output_taxon_file_name).st_size, "File size expected to be 515 bytes, but is %i." %os.stat(self.output_taxon_file_name).st_size)
        idList = []
        with open(self.output_taxon_file_name) as f:
            for line in f:
                idList.append(line.split()[0])
        self.assertIn('358301', idList, "Expected taxonomy ID 358301 is not found in the file")
        self.assertIn('358327', idList, "Expected taxonomy ID 358327 is not found in the file")
        self.assertIn('358485', idList, "Expected taxonomy ID 358485 is not found in the file")
      
      
    def test_execute_script_raw_taxonomy(self):
        #OK. Time: 3.40s

        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "raw", "-t", "--sampleName", "INVITRO221"], 
                             self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertTrue(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was not created.")
        #Better test for the number of sequences in the file        
        #self.assertEqual(608, os.stat(self.output_taxon_file_name).st_size, "File size expected to be 608 bytes, but is %i." %os.stat(self.output_taxon_file_name).st_size)
        idList = []
        with open(self.output_taxon_file_name) as f:
            for line in f:
                idList.append(line.split()[0])
        self.assertIn('961', idList, "Expected taxonomy ID 961 is not the found in the file")
        self.assertIn('97830', idList, "Expected taxonomy ID 97830 is not found in the file")
        self.assertIn('97847', idList, "Expected taxonomy ID 97847 is not found in the file")
        
        
    def test_execute_script_all_taxonomy(self):
        #OK. Time: 124.98s
        
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "-t", "all", "--geneRegion", "ACA"], 
                             self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertTrue(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        #Better test for the number of sequences in the file
        #self.assertEqual(111174, os.stat(self.output_taxon_file_name).st_size, "File size expected to be 111174 bytes, but is %i." %os.stat(self.output_taxon_file_name).st_size)
        idList = []
        with open(self.output_taxon_file_name) as f:
            for line in f:
                idList.append(line.split()[0])
        self.assertIn('358301', idList, "Expected taxonomy ID 358301 is not found in the file")
        self.assertIn('37674', idList, "Expected taxonomy ID 37674 is not found in the file")
        self.assertIn('148710', idList, "Expected taxonomy ID 148710 is not found in the file")

    
    ### TESTING ITS SEQUENCES
        
    def test_execute_script_its(self):
        #OK. Time: 963.44s
        
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "its"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        #Better test for the number of sequences in the file
        #self.assertEqual(19377145, os.stat(self.output_fasta_file_name).st_size, "File size expected to be 19377145 bytes, but is %i." %os.stat(self.output_fasta_file_name).st_size)
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    idList.append(line.split()[0])
        self.assertIn('>seqdb|131072', idList, "Expected sequence ID 131072 is not found in the file")
        self.assertIn('>seqdb|111872', idList, "Expected sequence ID 11187 is not found in the file")
        self.assertIn('>seqdb|131058', idList, "Expected sequence ID 131071 is not found in the file")             


if __name__ == "__main__":
    unittest.main()
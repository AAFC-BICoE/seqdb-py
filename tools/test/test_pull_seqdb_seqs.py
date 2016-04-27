'''
Created on Oct 16, 2015

@author: korolo
'''
import os.path
import unittest

from config import config_root
from tools import pull_seqdb_seqs


class TestPullSeqdbSeqs(unittest.TestCase):

    @classmethod
    def setUpClass(self):
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
    #Method no longer used'
    
    ** NOTE: pull_seqdb_seqs.get_seq_ids is removed. TODO: make sure that the equivalent tests exist in API tests for filters

    # The remaining tests for "All Sequences" do not exist yet in API tests for filters. This is why they have been left commented out in here.

    def test_get_seq_ids_specimen(self):  
     
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", [4405])
        self.assertEqual(23, len(seq_ids), "Expected 23 sequences, but got {}.".format(len(seq_ids)))


    def test_get_seq_ids_sequenceName(self):  

        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", sequence_name="Y10_16_F167")
        self.assertEqual(2, len(seq_ids), "Expected 2 sequences, but got %i" % len(seq_ids))
        
        
    def test_get_ids_sampleName(self):

        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", sample_name="LEV6103")
        self.assertEquals(61, len(seq_ids), "Expected 61 sequences, but got %i." % len(seq_ids))


    def test_get_seq_ids_regionName(self): 
 
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", region_name="ACA")  
        self.assertEqual(1041, len(seq_ids), "Expected 1,041 sequences, but got %i." % len(seq_ids))

    
    #Note: Failing. Error message: exceptions must be old-style classes or derived from BaseException, not str
    def test_get_seq_ids_projectName(self):  
    
       # all 
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", project_name="Pythium Type Specimens")
        self.assertEqual(4373, len(seq_ids), "Expected 4,373 sequences, but got %i." % len(seq_ids))
    
    
    def test_get_seq_ids_colletionCode(self): 
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", collection_code="BISS")
        self.assertEqual(25, len(seq_ids), "Expected 25 sequences, but got %i." % len(seq_ids))

    def test_get_seq_ids_taxonomy(self):    
        
        # all
        seq_ids = pull_seqdb_seqs.get_seq_ids(self.rawSeqFixture, self.consensusSeqFixture, "all", taxonomy_rank="species", taxonomy_value="megasperma")
        self.assertEqual(218, len(seq_ids), "Expected 218 sequences, but got %i." % len(seq_ids))
    '''


    ### TESTING FASTA FILE CREATION 

    def test_execute_script_consensus_fasta(self):
        
        # Getting all Consensus Sequences. Time: 56.56s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "consensus"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):     
                    count = count + 1
                    idList.append(line.split()[0])
        self.assertEqual(5593, count, "Expected 5,593 sequences but got {}".format(count))
        self.assertIn('>seqdb|358301', idList, "Expected sequence ID 358301 is not found in the file")
        self.assertIn('>seqdb|4823203', idList, "Expected sequence ID 4823203 is not found in the file")
        self.assertIn('>seqdb|4829279', idList, "Expected sequence ID 4829279 is not found in the file")
        
        
        # Filtering on Gene Region Name. Time: 2.71s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "consensus", "--geneRegion", "28s"],
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):     
                    count = count + 1
                    idList.append(line.split()[0])
        self.assertEqual(69, count, "Expected 69 sequences but got {}".format(count))
        self.assertIn('>seqdb|1582548', idList, "Expected sequence ID 358301 is not found in the file")
        self.assertIn('>seqdb|4825579', idList, "Expected sequence ID 4823203 is not found in the file")
        self.assertIn('>seqdb|4827758', idList, "Expected sequence ID 4829279 is not found in the file")


        # Filtering on Specimen Identifier. Time: 1.32s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "consensus", "--specNums", "4405,4264"],
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):     
                    count = count + 1
                    idList.append(line.split()[0])
        self.assertEqual(3, count, "Expected 3 sequences but got {}".format(count))
        self.assertIn('>seqdb|358301', idList, "Expected sequence ID 358301 is not found in the file")
        self.assertIn('>seqdb|358302', idList, "Expected sequence ID 358301 is not found in the file")
        self.assertIn('>seqdb|4825628', idList, "Expected sequence ID 358301 is not found in the file")
              
        
    def test_execute_script_raw_fasta(self):
        
        '''
        #Getting all Raw Sequences. Time: TOO LONG
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "raw"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    count = count + 1
        self.assertEqual(480088, count, "Expected 480,088 sequences but got {}".format(count))                
        '''
        
        # Filtering on Sequence Name. Time: 26.66s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "raw", "--seqName", "S-SH-"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    count = count + 1
                    idList.append(line.split()[0])
        self.assertEqual(134, count, "Expected 134 sequences but got {}".format(count))        
        self.assertIn('>seqdb|1', idList, "Expected sequence ID 1 is not found in the file")
        self.assertIn('>seqdb|79390', idList, "Expected sequence ID 79390 is not found in the file")
        self.assertIn('>seqdb|126059', idList, "Expected sequence ID 126059 is not found in the file")
    
    
        # Filtering on Collection Code. Time: 1.91s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "raw", "--collectionCode", "pm"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    count = count + 1
                    idList.append(line.split()[0])
        self.assertEqual(148, count, "Expected 148 sequences but got {}".format(count))        
        self.assertIn('>seqdb|268749', idList, "Expected sequence ID 1 is not found in the file")
        self.assertIn('>seqdb|308734', idList, "Expected sequence ID 79390 is not found in the file")
        self.assertIn('>seqdb|356572', idList, "Expected sequence ID 126059 is not found in the file")
        
        
        # Filtering on Specimen Identifier. Time: 6.03s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "raw", "--specNums", "4405,4264"],
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):     
                    count = count + 1
                    idList.append(line.split()[0])
        self.assertEqual(33, count, "Expected 33 sequences but got {}".format(count))
        self.assertIn('>seqdb|27755', idList, "Expected sequence ID 358301 is not found in the file")
        self.assertIn('>seqdb|155033', idList, "Expected sequence ID 358301 is not found in the file")
        self.assertIn('>seqdb|239733', idList, "Expected sequence ID 358301 is not found in the file")
        
 
    def test_execute_script_all_fasta(self):
        
        '''
        #Getting All Sequences. Time: TOO LONG
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "all"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    count = count + 1
        self.assertEqual(485643, count, "Expected 485,643 sequences but got {}".format(count))     
        '''
                
        # Filtering on Gene Region Name. Time: 28.82s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "all", "--geneRegion", "EF-1a"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    count = count + 1
                    idList.append(line.split()[0]) 
        self.assertEqual(492, count, "Expected 492 sequences but got {}".format(count))     
        self.assertIn('>seqdb|1689', idList, "Expected sequence ID 1689 is not found in the file")
        self.assertIn('>seqdb|103372', idList, "Expected sequence ID 103372 is not found in the file")
        self.assertIn('>seqdb|149807', idList, "Expected sequence ID 149807 is not found in the file")
        
        
        # Filtering on Project Name. Time: 83.99s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "all", "--projectName", "Pythium Type Specimens"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    count = count + 1
                    idList.append(line.split()[0]) 
        self.assertEqual(4373, count, "Expected 4,373 sequences but got {}".format(count))     
        self.assertIn('>seqdb|358305', idList, "Expected sequence ID 1689 is not found in the file")
        self.assertIn('>seqdb|196715', idList, "Expected sequence ID 103372 is not found in the file")
        self.assertIn('>seqdb|356858', idList, "Expected sequence ID 149807 is not found in the file")
        
        
        
    ### TESTING FASTQ FILE CREATION
        
    def test_execute_script_raw_fastq(self):

        # Filtering on Sample Name. Time: 5.16s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fastq", "raw", "--sampleName", "LEV6103"], 
                             self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fastq_file_name), "Fastq file was not created.")
        self.assertFalse(os.path.isfile(self.output_fasta_file_name), "Fasta file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        idList = []
        with open(self.output_fastq_file_name) as f:
            for line in f:
                if line.startswith('@'):
                    count = count + 1
                    idList.append(line.split()[0])
        self.assertEqual(60, count, "Expected 60 sequences but got {}".format(count))
        self.assertIn('@seqdb|266400', idList, "Expected sequence ID 266400 is not found in the file")
        self.assertIn('@seqdb|301609', idList, "Expected sequence ID 301609 is not found in the file")
        self.assertIn('@seqdb|331086', idList, "Expected sequence ID 331086 is not found in the file")   

 
     
    ### TESTING TAXONOMY FILE CREATION
    
    def test_execute_script_consensus_taxonomy(self):    

        # Filtering on Sequence Name. Time: 1.06s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "-t", "consensus", "--seqName", "Pyt_arrhenomanes_"], 
                                self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertTrue(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        count = 0
        idList = []
        with open(self.output_taxon_file_name) as f:
            for line in f:
                count = count + 1
                idList.append(line.split()[0])
        self.assertEqual(5, count, "Expected 5 sequence but got {}".format(count))
        self.assertIn('358301', idList, "Expected taxonomy ID 358301 is not found in the file")
        self.assertIn('358327', idList, "Expected taxonomy ID 358327 is not found in the file")
        self.assertIn('358485', idList, "Expected taxonomy ID 358485 is not found in the file")

        
        # Filtering on Taxonomy Rank and Taxonomy Value. Time: 2.49s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "-t", "consensus", "--taxRank", "species", "--taxValue", "megasperma"], 
                                self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertTrue(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        count = 0
        idList = []
        with open(self.output_taxon_file_name) as f:
            for line in f:
                count = count + 1
                idList.append(line.split()[0])
        self.assertEqual(3, count, "Expected 3 sequence but got {}".format(count))
        self.assertIn('358368', idList, "Expected taxonomy ID 358301 is not found in the file")
        self.assertIn('358385', idList, "Expected taxonomy ID 358327 is not found in the file")
        self.assertIn('358394', idList, "Expected taxonomy ID 358485 is not found in the file")
      
      
    def test_execute_script_raw_taxonomy(self):

        # Filtering on Sample Name. Time: 3.40s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "raw", "-t", "--sampleName", "INVITRO221"], 
                             self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertTrue(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was not created.")
        count = 0
        idList = []
        with open(self.output_taxon_file_name) as f:
            for line in f:
                count = count + 1
                idList.append(line.split()[0])
        self.assertEqual(6, count, "Expected 6 sequences but got {}".format(count))
        self.assertIn('961', idList, "Expected taxonomy ID 961 is not the found in the file")
        self.assertIn('97830', idList, "Expected taxonomy ID 97830 is not found in the file")
        self.assertIn('97847', idList, "Expected taxonomy ID 97847 is not found in the file")

          
    def test_execute_script_all_taxonomy(self):
        
        # Filtering on Gene Region Name. Time: 124.98s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "-t", "all", "--geneRegion", "ACA"], 
                             self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertTrue(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        count = 0
        idList = []
        with open(self.output_taxon_file_name) as f:
            for line in f:
                count = count + 1
                idList.append(line.split()[0])
        self.assertEqual(1041, count, "Expected 1,041 sequences but got {}".format(count))
        self.assertIn('358301', idList, "Expected taxonomy ID 358301 is not found in the file")
        self.assertIn('37674', idList, "Expected taxonomy ID 37674 is not found in the file")
        self.assertIn('148710', idList, "Expected taxonomy ID 148710 is not found in the file")


    
    ### TESTING ITS SEQUENCES
        
    def test_execute_script_its(self):
        
        # Getting all ITS Sequences. Time: 963.44s
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "its"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertFalse(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was created.")
        count = 0
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    count = count + 1
                    idList.append(line.split()[0])
        self.assertEqual(17787, count, "Expected 17787 sequences but got {}".format(count))
        self.assertIn('>seqdb|131072', idList, "Expected sequence ID 131072 is not found in the file")
        self.assertIn('>seqdb|111872', idList, "Expected sequence ID 11187 is not found in the file")
        self.assertIn('>seqdb|131058', idList, "Expected sequence ID 131071 is not found in the file")             


if __name__ == "__main__":
    unittest.main()
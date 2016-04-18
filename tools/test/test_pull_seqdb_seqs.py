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
    
    ''' Method removed    
    def test_get_seq_ids_consensus(self):
        # time: 8.395s 
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="consensus")
        self.assertEqual(5555 , len(seq_ids), "Expected 5,555 consensus sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    '''
    
    """
    def test_get_seq_ids_raw(self):
        # time: stopped
        seq_ids = pull_seqdb_seqs.get_seq_ids(seqdbWS=self.fixture, pull_type="raw")
        self.assertEquals(480088, len(seq_ids), "Expected 480,088 raw sequences, but got %i. Doublecheck test db to make sure the numbers haven't changed there." % len(seq_ids))
    """    
    def test_execute_script_consensus(self):
        # time: 132.482s
        
        # Testing fasta file creation for consensus sequences 
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "consensus"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertEqual(13983492, os.stat(self.output_fasta_file_name).st_size, "File size expected to be 13983492 bytes, but is %i." %os.stat(self.output_fasta_file_name).st_size)
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                idList.append(line.split()[0])
        self.assertEqual('>seqdb|358301', idList[0], "Expected sequence ID 358301 is not the same as the ID %s in the file." %idList[0])
        self.assertEqual('>seqdb|4823203', idList[89597], "Expected sequence ID 4823203 is not the same as the ID %s in the file." %idList[89597])
        self.assertEqual('>seqdb|4829279', idList[161610], "Expected sequence ID 4829279 is not the same as the ID %s in the file." %idList[161610])
        
        # Testing taxonomy file creation for consensus sequences
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "-t", "consensus", "--seqName", "Pyt_arrhenomanes_"], 
                                self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertTrue(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertEqual(515, os.stat(self.output_taxon_file_name).st_size, "File size expected to be 515 bytes, but is %i." %os.stat(self.output_taxon_file_name).st_size)
        idList = []
        with open(self.output_taxon_file_name) as f:
            for line in f:
                idList.append(line.split()[0])
        self.assertEqual('358301' , idList[0], "Expected taxonomy ID 358301 is not the same as the ID %s in the file." %idList[0] )
        self.assertEqual('358327' , idList[2], "Expected taxonomy ID 358327 is not the same as the ID %s in the file." %idList[2])
        self.assertEqual('358485' , idList[4], "Expected taxonomy ID 358485 is not the same as the ID %s in the file." %idList[4])
        
        
    def test_execute_script_raw(self):
        # time: 44.538s
        
        # Testing fasta file creation for raw sequences 
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "raw", "--seqName", "S-SH-"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        
        # Very arbitrary test. Better test for the number of sequences in the file
        #self.assertEqual(138996, os.stat(self.output_fasta_file_name).st_size, "File size expected to be 138996 bytes, but is %i." %os.stat(self.output_fasta_file_name).st_size)
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                idList.append(line.split()[0])        
        self.assertEqual('>seqdb|1' , idList[0], "Expected sequence ID 1 is not the same as the ID %s in the file." %idList[0])
        self.assertEqual('>seqdb|79390', idList[1026], "Expected sequence ID 79390 is not the same as the ID %s in the file." %idList[1026])
        self.assertEqual('>seqdb|126059', idList[1838], "Expected sequence ID 126059 is not the same as the ID %s in the file." %idList[1838])
        
        # Testing fastq file creation for raw sequences 
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fastq", "raw", "--sampleName", "LEV6103"], 
                             self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fastq_file_name), "Fastq file was not created.")
        #self.assertEqual(106623, os.stat(self.output_fastq_file_name).st_size, "File size expected to be 106623 bytes, but is %i." %os.stat(self.output_fastq_file_name).st_size)
        idList = []
        with open(self.output_fastq_file_name) as f:
            for line in f:
                idList.append(line.split()[0])
        self.assertEqual('@seqdb|266400', idList[0], "Expected sequence ID 266400 is not the same as the ID %s in the file." %idList[0])
        self.assertEqual('@seqdb|301609', idList[140], "Expected sequence ID 301609 is not the same as the ID %s in the file." %idList[140])
        self.assertEqual('@seqdb|331086', idList[228], "Expected sequence ID 331086 is not the same as the ID %s in the file." %idList[228])
        
        # Testing taxonomy file creation for raw sequences
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fastq", "raw", "-t", "--sampleName", "LEV6103"], 
                             self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fastq_file_name), "Fastq file was not created.")
        self.assertTrue(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was not created.")
        #self.assertEqual(8822, os.stat(self.output_taxon_file_name).st_size, "File size expected to be 8822 bytes, but is %i." %os.stat(self.output_taxon_file_name).st_size)
        idList = []
        with open(self.output_taxon_file_name) as f:
            for line in f:
                idList.append(line.split()[0])
        self.assertEqual('266400' , idList[0], "Expected taxonomy ID 266400 is not the same as the ID %s in the file." %idList[0])
        self.assertEqual('299528' , idList[32], "Expected taxonomy ID 299528 is not the same as the ID %s in the file." %idList[32])
        self.assertEqual('4777769' , idList[59], "Expected taxonomy ID 4777769 is not the same as the ID %s in the file." %idList[59])
        
        
    def test_execute_script_all(self):
        # time: 212.517s
        
        # Testing fasta file creation for all sequences
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "all", "--geneRegion", "ACA"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertEqual(1129727, os.stat(self.output_fasta_file_name).st_size, "File size expected to 1129727 bytes, but is %i." %os.stat(self.output_fasta_file_name).st_size)
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    idList.append(line.split()[0])
        self.assertIn('>seqdb|6818', idList, "Expected sequence ID 6818 is not found the file." )        
        self.assertIn('>seqdb|112673', idList, "Expected sequence ID 112673 is not found the file." )        
        self.assertIn('>seqdb|358505', idList, "Expected sequence ID 358505 is not found the file." )        
        
        # Testing taxonomy file creation for all sequences
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "-t", "all", "--geneRegion", "ACA"], 
                             self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertTrue(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertEqual(111174, os.stat(self.output_taxon_file_name).st_size, "File size expected to be 111174 bytes, but is %i." %os.stat(self.output_taxon_file_name).st_size)
        idList = []
        with open(self.output_taxon_file_name) as f:
            for line in f:
                idList.append(line.split()[0])
        self.assertIn('6818', idList, "Expected sequence ID 6818 is not found the file." )        
        self.assertIn('112600', idList, "Expected sequence ID 112600 is not found the file." )        
        self.assertIn('358502', idList, "Expected sequence ID 358502 is not found the file." )        
        
        
    def test_execute_script_taxonomy(self):
        # time: 13.551s
        
        # Testing fasta file creation for consensus sequences
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "consensus", "--taxRank", "species", "--taxValue", "arrhenomanes"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertEqual(6063, os.stat(self.output_fasta_file_name).st_size, "File size expected to be 6063 bytes, but is %i." %os.stat(self.output_fasta_file_name).st_size)
        idList = []
        with open(self.output_fasta_file_name) as f:
            for line in f:
                if line.startswith('>'):
                    idList.append(line.split()[0]) 
        self.assertIn('>seqdb|358301', idList, "Expected sequence ID 358301 is not found the file." )        
        self.assertIn('>seqdb|358381', idList, "Expected sequence ID 358381 is not found the file." )        
        self.assertIn('>seqdb|358485', idList, "Expected sequence ID 358485 is not found the file." )        
        
        # Testing taxonomy file creation for raw sequences
        pull_seqdb_seqs.execute_script(["-c", config_root.path() + '/config4tests.yaml', "-r", "fasta", "-t", "consensus", "--taxRank", "genus", "--taxValue", "Pythium"], 
                            self.output_file_name, self.output_taxon_file_name)
        self.assertTrue(os.path.isfile(self.output_fasta_file_name), "Fasta file was not created.")
        self.assertTrue(os.path.isfile(self.output_taxon_file_name), "Taxonomy file was not created.")
        self.assertFalse(os.path.isfile(self.output_fastq_file_name), "Fastq file was created.")
        self.assertEqual(11293, os.stat(self.output_taxon_file_name).st_size, "File size expected to be 11293 bytes, but is %i." %os.stat(self.output_taxon_file_name).st_size)
        idList = []
        with open(self.output_taxon_file_name) as f:
            for line in f:
                idList.append(line.split()[0])
        self.assertIn('358301', idList, "Expected sequence ID 358301 is not found the file." )        
        self.assertEqual('358301' , idList[0], "Expected taxonomy ID 358301 is not the same as the ID %s in the file." %idList[0])
        self.assertEqual('358416' , idList[62], "Expected taxonomy ID 358416 is not the same as the ID %s in the file." %idList[62])
        self.assertEqual('1582548' , idList[108], "Expected taxonomy ID 1582548 is not the same as the ID %s in the file." %idList[108])
        
    
    """ pull_seqdb_seqs.get_seq_ids is removed. TODO: make sure that the equivalet test exists in api tests for filters
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
 
if __name__ == "__main__":
    unittest.main()
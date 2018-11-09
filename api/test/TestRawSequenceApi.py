'''
Created on Apr 1, 2016

@author: korolo
'''
import unittest

import yaml

from api.RawSequenceApi import RawSequenceApi
from config import config_root


class TestRawSequenceApi(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))
        self.fixture = RawSequenceApi(api_key=config['seqdb']['api_key'],
                                                   base_url=config['seqdb']['base_url'])


    def setUp(self):
        pass


    def tearDown(self):
        pass

    def testFilters(self):
        # specimenNumber - tested below, skipping

        # sequenceName
        self.fixture.sequenceNameFilter = "SH-254"
        self.assertEqual([13808,65268], self.fixture.getIds(), "Expecting 2 raw sequences filtered by sequenceName = SH-254")
        self.fixture.clearAllFilters()

        self.fixture.sequenceNameFilter = "Y10_16_F167"
        self.assertEqual(2, self.fixture.getNumber(), "Expecting 2 raw sequences filtered by sequenceName = Y10_16_F167, but got {}.".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()

        # sampleNameFilter
        self.fixture.sampleNameFilter = "LEV4277"
        self.assertEqual(15, self.fixture.getNumber(), "Expecting 15 raw sequences filtered by sampleName = LEV4277, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()

        self.fixture.sampleNameFilter = "LEV6103"
        self.assertEqual(60, self.fixture.getNumber(), "Expecting 60 raw sequences filtered by sampleName = LEV6103, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()

        # pubRefSeqFilter
        #TODO
        '''
        self.fixture.pubRefSeqFilter = True
        #TODO: this fails, i.e. curl -H "apikey: ***REMOVED***" "***REMOVED***?filterName=sequence.submittedToInsdc&filterValue=true&filterWildcard=false"
        # Investigate why
        va = self.fixture.getNumber()
        self.assertEqual(15, self.fixture.getNumber(),
                         "Expecting 15 raw sequences filtered by pubRefSeq = , but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()
        '''

        # genBankGIFilter
        #TODO: fix this test
        '''
        self.fixture.genBankGIFilter = "gi_"
        self.assertEqual(15, self.fixture.getNumber(),
                         "Expecting 15 raw sequences filtered by genBankGI = LEV4277, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()
        '''

        # regionNameFilter
        self.fixture.regionNameFilter = "ITS/LR6"
        self.assertEqual(58, self.fixture.getNumber(), "Expecting 58 raw sequences filtered by regionName = ITS/LR6, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()

        self.fixture.regionNameFilter = "ef-1a"
        self.assertEqual(492, self.fixture.getNumber(), "Expecting 492 raw sequences filtered by regionName = ef-1a, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()

        # collectionCodeFilter
        self.fixture.collectionCodeFilter = "lev"
        self.assertEqual(33156, self.fixture.getNumber(), "Expecting 33156 raw sequences filtered by collectionCode = lev, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()

        self.fixture.collectionCodeFilter = "pm"
        self.assertEqual(148, self.fixture.getNumber(), "Expecting 148 raw sequences filtered by collectionCode = pm, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()

        # taxonomyRankFilter
        self.fixture.taxonomyRankFilter = "species"
        self.fixture.taxonomyValueFilter = "megasperma"
        self.assertEqual(215, self.fixture.getNumber(), "Expecting 215 raw sequences filtered by taxonomy, but got {}".format(self.fixture.getNumber()))
        self.fixture.clearAllFilters()

        # TODO: test combinations of filters

    def testGetFastaSequencesWithOffset(self):
        self.fixture.specimenNumFilter = 4405
        # This filter should have 22 raw sequences associated with it. Therefore there should
        # be 2 calls with limit 15
        actual_fasta, result_offset = self.fixture.getFastaSequencesWithOffset(offset=0, limit=15)
        self.assertTrue(actual_fasta, "No Sequences returned.")
        self.assertIn(">seqdb|27755", actual_fasta,"Expecting that fasta return will contain id 27755.")
        self.assertNotIn(">seqdb|358301", actual_fasta, "Fasta return is not expected to have sequence 358301, since it is consensus.")
        self.assertEquals(result_offset, 15, "First offset should be 15 (there are 22 raw sequences for specimenNum=4405).")

        actual_fasta, result_offset = self.fixture.getFastaSequencesWithOffset(offset=result_offset, limit=15)
        self.assertTrue(actual_fasta, "No Sequences returned.")
        self.assertIn(">seqdb|160946", actual_fasta,"Expecting that fasta return will contain id 160946.")
        self.assertNotIn(">seqdb|358301", actual_fasta, "Fasta return is not expected to have sequence 358301, since it is consensus.")
        self.assertEquals(result_offset, -1, "Second offset should be -1, since there are no more sequences to retrieve.")

    def testGetFastqSequencesWithOffset(self):
        self.fixture.specimenNumFilter = 4405
        # This filter should have 22 raw sequences associated with it. Therefore there should
        # be 2 calls with limit 15
        actual_fastq, result_offset = self.fixture.getFastqSequencesWithOffset(offset=0, limit=15)
        self.assertTrue(actual_fastq, "No Sequences returned.")
        self.assertIn("@seqdb|27755", actual_fastq,"Expecting that fasta return will contain id 27755.")
        self.assertNotIn("@seqdb|358301", actual_fastq, "Fasta return is not expected to have sequence 358301, since it is consensus.")
        self.assertEquals(result_offset, 15, "First offset should be 15 (there are 22 raw sequences for specimenNum=4405).")

        actual_fastq, result_offset = self.fixture.getFastqSequencesWithOffset(offset=result_offset, limit=15)
        self.assertTrue(actual_fastq, "No Sequences returned.")
        self.assertIn("@seqdb|160946", actual_fastq,"Expecting that fasta return will contain id 160946.")
        self.assertNotIn("@seqdb|358301", actual_fastq, "Fasta return is not expected to have sequence 358301, since it is consensus.")
        self.assertEquals(result_offset, -1, "Second offset should be -1, since there are no more sequences to retrieve.")

    def testGetSequenceIds(self):
        self.fixture.specimenNumFilter = 4405
        actual = self.fixture.getIds()
        self.assertTrue(actual, "No Sequence ids returned.")
        self.assertEqual(22, len(actual),"Expecting 22 sequences associated with this specimen.")
        self.assertIn(27755, actual, "Sequence id 27755 is expected to be associated with specimen 4405.")
        self.assertNotIn(358301, actual, "Sequence id 358301 is not expected to be in results, since it is consensus.")


    def testCreateChromatSequence_wrong_path(self):
        self.assertRaises(IOError, self.fixture.importChromatSequencesFromFile, "data/")
        self.assertRaises(IOError, self.fixture.importChromatSequencesFromFile, "zzz/non-existent.ab1")

    def testCreateDeleteChromatSequence(self):
        """ Test creating a sequence with binary .abi or .ab1 file (chromatogram) """
        seq_id = self.fixture.importChromatSequencesFromFile(chromat_file = "data/GRDI_test_seq.ab1", notes="This is a test upload.",
                    trace_file_path="test_path",dest_file_name="test.ab1",)

        self.assertTrue(seq_id, "Persisting chromatogram did not return an id.")

        # Delete
        delete_jsn_resp = self.fixture.delete(seq_id)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete feature type.")

    def testCreateDeleteChromatSequenceGz_valid(self):
        """ Test creating a sequence with a zipped binary file (chromatogram) i.e. ab1.gz """
        seq_id = self.fixture.importChromatSequencesFromFile(chromat_file = "data/blob_db.ab1.gz")
        self.assertTrue(seq_id, "Persisting chromatogram did not return an id.")

        # Delete
        delete_jsn_resp = self.fixture.delete(seq_id)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete feature type.")

    def testGetFastaSeq(self):
        actual = self.fixture.getFastaSequence("1")
        self.assertTrue(actual, "Fasta sequence is empty.")
        self.assertIn(">seqdb|1", actual, "Fasta does not contain >.")

    def testGetFastqSeq(self):
        actual = self.fixture.getFastqSequence("1")
        self.assertTrue(actual, "Fastq sequence is empty.")
        self.assertIn("@seqdb|1", actual, "Fastq does not contain @seqdb.")

    def testSequenceIdsByRegion(self):
        actual_seq_ids, result_offset = self.fixture.getSequenceIdsByRegionWithOffset(97)
        self.assertTrue(actual_seq_ids, "No sequence ids returned for region id=97.")
        self.assertEquals(20, len(actual_seq_ids), "Expecting 20 sequences, returned in the first offset query and associated with region id=97 , but got {}.".format(len(actual_seq_ids)))

    def testGetAcceptedSpecimenDetermination(self):
        actual = self.fixture.getAcceptedSpecimenDetermination(27755)
        self.assertTrue(actual, "Expecting accepted determination, but got none.")
        self.assertEquals("Fungi", actual['taxonomy']['kingdom'], "Expecting kingdom Fungi, but got {}".format(actual['taxonomy']['kingdom']))
        self.assertEquals("arrhenomanes", actual['taxonomy']['species'], "Expecting species arrhenomanes, but got {}".format(actual['taxonomy']['species']))

        actual = self.fixture.getAcceptedSpecimenDetermination(358301)
        self.assertTrue(actual, "Expecting accepted determination, but got none.")
        self.assertEquals("Oomycota", actual['taxonomy']['phylum'], "Expecting phylum Oomycota, but got {}".format(actual['taxonomy']['phylum']))
        self.assertEquals("Pythiales", actual['taxonomy']['taxanomicOrder'], "Expecting order Pythiales, but got {}".format(actual['taxonomy']['taxanomicOrder']))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

"""
Created on Apr 1, 2016

@author: korolo
"""
import unittest
import os
import yaml
from context import RawSequenceApi
from config import config_root


class TestRawSequenceApi(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(config_root.path() +
                  '/config4tests.yaml', 'r') as config_file:
            config = yaml.safe_load(config_file)
            cls.fixture = RawSequenceApi(
                api_key=config['seqdb']['api_key'],
                base_url=config['seqdb']['base_url'])
            cls.dir_path=os.path.dirname(os.path.realpath(__file__))
            #cls.base_path = config['basepath']

    def setUp(self):
        pass

    def tearDown(self):
        pass

    """
    def testFilters(self):
        # specimenNumber - tested below, skipping

        # sequenceName
        self.fixture.sequence_name_filter = 'SH-254'
        self.assertEqual([13808, 65268], self.fixture.get_ids(),
                         'Expecting 2 raw sequences filtered by '
                         'sequenceName = SH-254')
        self.fixture.clear_all_filters()

        self.fixture.sequence_name_filter = 'Y10_16_F167'
        self.assertEqual(2, self.fixture.get_number(),
                         'Expecting 2 raw sequences filtered by '
                         'sequenceName = Y10_16_F167, but got {}.'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        # sample_name_filter
        self.fixture.sample_name_filter = 'LEV4277'
        self.assertEqual(15, self.fixture.get_number(),
                         'Expecting 15 raw sequences filtered by '
                         'sampleName = LEV4277, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        self.fixture.sample_name_filter = 'LEV6103'
        self.assertEqual(60, self.fixture.get_number(),
                         'Expecting 60 raw sequences filtered by '
                         'sampleName = LEV6103, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        # pub_ref_seq_filter
        # TODO
        '''
        self.fixture.pub_ref_seq_filter = True
        #TODO: this fails, i.e. curl -H 'apikey: ***REMOVED***' '***REMOVED***?
        filterName=sequence
        .submittedToInsdc&filterValue=true&filterWildcard=false'
        # Investigate why
        va = self.fixture.get_number()
        self.assertEqual(15, self.fixture.get_number(),
                         'Expecting 15 raw sequences filtered by pubRefSeq = , 
                         but got {}'.format(self.fixture.get_number()))
        self.fixture.clear_all_filters()
        '''

        # gen_bank_GI_filter
        # TODO: fix this test
        '''
        self.fixture.gen_bank_GI_filter = 'gi_'
        self.assertEqual(15, self.fixture.get_number(),
                         'Expecting 15 raw sequences filtered by 
                         genBankGI = LEV4277, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()
        '''

        # region_name_filter
        self.fixture.region_name_filter = 'ITS/LR6'
        self.assertEqual(58, self.fixture.get_number(),
                         'Expecting 58 raw sequences filtered by '
                         'regionName = ITS/LR6, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        self.fixture.region_name_filter = 'ef-1a'
        self.assertEqual(492, self.fixture.get_number(),
                         'Expecting 492 raw sequences filtered by '
                         'regionName = ef-1a, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        # collection_code_filter
        self.fixture.collection_code_filter = 'lev'
        self.assertEqual(33156, self.fixture.get_number(),
                         'Expecting 33156 raw sequences filtered by '
                         'collectionCode = lev, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        self.fixture.collection_code_filter = 'pm'
        self.assertEqual(148, self.fixture.get_number(),
                         'Expecting 148 raw sequences filtered by '
                         'collectionCode = pm, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        # taxonomy_rank_filter
        self.fixture.taxonomy_rank_filter = 'species'
        self.fixture.taxonomy_value_filter = 'megasperma'
        self.assertEqual(215, self.fixture.get_number(),
                         'Expecting 215 raw sequences filtered by '
                         'taxonomy, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        # TODO: test combinations of filters

    def test_get_fasta_sequences_with_offset(self):
        self.fixture.specimen_num_filter = 4405
        # This filter should have 22 raw sequences associated with it.
        # Therefore there should be 2 calls with limit 15
        actual_fasta, result_offset = \
            self.fixture.get_fasta_sequences_with_offset(offset=0, limit=15)
        self.assertTrue(actual_fasta, 'No Sequences returned.')
        self.assertIn('>seqdb|27755', actual_fasta,
                      'Expecting that fasta return will contain id 27755.')
        self.assertNotIn('>seqdb|358301', actual_fasta,
                         'Fasta return is not expected to have sequence '
                         '358301, since it is consensus.')
        self.assertEquals(result_offset, 15,
                          'First offset should be 15 (there are 22 raw '
                          'sequences for specimenNum=4405).')

        actual_fasta, result_offset = \
            self.fixture.get_fasta_sequences_with_offset(
                offset=result_offset, limit=15)
        self.assertTrue(actual_fasta, 'No Sequences returned.')
        self.assertIn('>seqdb|160946', actual_fasta,
                      'Expecting that fasta return will contain id 160946.')
        self.assertNotIn('>seqdb|358301', actual_fasta,
                         'Fasta return is not expected to have sequence '
                         '358301, since it is consensus.')
        self.assertEquals(result_offset, -1,
                          'Second offset should be -1, since there '
                          'are no more sequences to retrieve.')

    def test_get_fastq_sequences_with_offset(self):
        self.fixture.specimen_num_filter = 4405
        # This filter should have 22 raw sequences associated with it.
        # Therefore there should be 2 calls with limit 15
        actual_fastq, result_offset = \
            self.fixture.get_fastq_sequences_with_offset(offset=0, limit=15)
        self.assertTrue(actual_fastq, 'No Sequences returned.')
        self.assertIn('@seqdb|27755', actual_fastq,
                      'Expecting that fasta return will contain id 27755.')
        self.assertNotIn('@seqdb|358301', actual_fastq,
                         'Fasta return is not expected to have sequence '
                         '358301, since it is consensus.')
        self.assertEquals(result_offset, 15,
                          'First offset should be 15 (there are 22 '
                          'raw sequences for specimenNum=4405).')

        actual_fastq, result_offset = \
            self.fixture.get_fastq_sequences_with_offset(offset=result_offset,
                                                         limit=15)
        self.assertTrue(actual_fastq, 'No Sequences returned.')
        self.assertIn('@seqdb|160946', actual_fastq,
                      'Expecting that fasta return will contain id 160946.')
        self.assertNotIn('@seqdb|358301', actual_fastq,
                         'Fasta return is not expected to have sequence '
                         '358301, since it is consensus.')
        self.assertEquals(result_offset, -1,
                          'Second offset should be -1, since there are '
                          'no more sequences to retrieve.')

    def test_get_sequence_ids(self):
        self.fixture.specimen_num_filter = 4405
        actual = self.fixture.get_ids()
        self.assertTrue(actual, 'No Sequence ids returned.')
        self.assertEqual(22, len(actual),
                         'Expecting 22 sequences associated with '
                         'this specimen.')
        self.assertIn(27755, actual,
                      'Sequence id 27755 is expected to be associated '
                      'with specimen 4405.')
        self.assertNotIn(358301, actual,
                         'Sequence id 358301 is not expected to be in '
                         'results, since it is consensus.')

    def test_create_chromat_sequence_wrong_path(self):
        self.assertRaises(IOError,
                          self.fixture.import_chromat_sequences_from_file,
                          'data/')
        self.assertRaises(IOError,
                          self.fixture.import_chromat_sequences_from_file,
                          'zzz/non-existent.ab1')
    """

    def test_create_delete_chromat_sequence(self):
        """
        Test creating a sequence with binary .abi or .ab1 file (chromatogram)
        """
        seq_id = self.fixture.import_chromat_sequences_from_file(
            chromat_file='{}/data/asdfjka.abi'.format(self.dir_path),
            #chromat_file='{}/seqdb_py/api/test/data/asdfjka.abi'.format(self.base_path),
            notes='This is a test upload.',
            trace_file_path='test_path',
            dest_file_name='test2.abi'
        )

        self.assertTrue(seq_id, 'Persisting chromat did not return an id.')
        """
        # Delete

        delete_jsn_resp = self.fixture.delete(seq_id)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'],
                         'Could not delete feature type.')
        """

    def test_create_delete_chromat_sequence_gzip_valid(self):
        """
        Test creating a sequence with a zipped binary file
        (chromatogram) i.e. ab1.gz
        """

        seq_id = self.fixture.import_chromat_sequences_from_file(
            chromat_file='{}/data/blob_db.ab1.gz'.format(self.dir_path)
            #chromat_file='{}/seqdb_py/api/test/data/blob_db.ab1.gz'.format(self.base_path)
        )
        self.assertTrue(seq_id,
                        'Persisting chromatogram did not return an id.')

        # Delete
        # curl -X DELETE -H [header] [api] currently returns 404 resource not found
        
        #delete_jsn_resp = self.fixture.delete(seq_id)
        #self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'],
        #                 'Could not delete feature type.')
    """    
    
    def test_get_fasta_seq(self):
        actual = self.fixture.get_fasta_sequence('24')
        self.assertTrue(actual, 'Fasta sequence is empty.')
        self.assertIn('>seqdb|24', actual, 'Fasta does not contain >.')

    def test_get_fastq_seq(self):
        actual = self.fixture.get_fastq_sequence('24')
        self.assertTrue(actual, 'Fastq sequence is empty.')
        self.assertIn('@seqdb|24', actual, 'Fastq does not contain @seqdb.')

    def test_sequence_ids_by_region(self):
        actual_seq_ids, result_offset = \
            self.fixture.get_sequence_ids_by_region_with_offset(97)
        self.assertTrue(actual_seq_ids,
                        'No sequence ids returned for region id=97.')
        self.assertEquals(20, len(actual_seq_ids),
                          'Expecting 20 sequences, returned in the first '
                          'offset query and associated with region id=97 , '
                          'but got {}.'
                          .format(len(actual_seq_ids)))

    def test_get_accepted_specimen_determination(self):
        actual = self.fixture.get_accepted_specimen_determination(27755)
        self.assertTrue(actual,
                        'Expecting accepted determination, but got none.')
        self.assertEquals('Fungi', actual['taxonomy']['kingdom'],
                          'Expecting kingdom Fungi, but got {}'
                          .format(actual['taxonomy']['kingdom']))
        self.assertEquals('arrhenomanes', actual['taxonomy']['species'],
                          'Expecting species arrhenomanes, but got {}'
                          .format(actual['taxonomy']['species']))

        actual = self.fixture.get_accepted_specimen_determination(358301)
        self.assertTrue(actual,
                        'Expecting accepted determination, but got none.')
        self.assertEquals('Oomycota', actual['taxonomy']['phylum'],
                          'Expecting phylum Oomycota, but got {}'
                          .format(actual['taxonomy']['phylum']))
        self.assertEquals('Pythiales', actual['taxonomy']['taxanomicOrder'],
                          'Expecting order Pythiales, but got {}'
                          .format(actual['taxonomy']['taxanomicOrder']))
    """


if __name__ == '__main__':
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

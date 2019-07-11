"""
Created on Apr 1, 2016

@author: korolo
"""
import unittest

import yaml

from config import config_root
from api.ConsensusSequenceApi import ConsensusSequenceApi


class TestConsensusSequenceApi(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(config_root.path() +
                  '/config4tests.yaml', 'r') as config_file:
            config = yaml.safe_load(config_file)
            cls.fixture = ConsensusSequenceApi(
                api_key=config['seqdb']['api_key'],
                base_url=config['seqdb']['base_url'])

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_filters(self):
        # specimenNumber - tested below, skipping

        # sequenceName
        self.fixture.sequence_name_filter = 'Pyt_arrhenomanes_BR0'
        self.assertEqual([358381, 358485], self.fixture.get_ids(),
                         'Expecting 2 consensus sequences filtered by '
                         'sequenceName = Pyt_arrhenomanes_BR0')
        self.fixture.clear_all_filters()

        self.fixture.sequence_name_filter = 'Pyt_arrhenomanes_'
        self.assertEqual(5, self.fixture.get_number(),
                         'Expecting 5 consensus sequences filtered by '
                         'sequenceName = Pyt_arrhenomanes_, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        # sample_name_filter
        self.fixture.sample_name_filter = 'LEV5508'
        self.assertEqual(1, self.fixture.get_number(),
                         'Expecting 1 consensus sequence filtered by '
                         'sampleName = LEV5508, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        self.fixture.sample_name_filter = 'LEV4183'
        self.assertEqual(1, self.fixture.get_number(),
                         'Expecting 1 consensus sequence filtered by '
                         'sampleName = LEV4183, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        # pub_ref_seq_filter
        # TODO
        '''
        self.fixture.pub_ref_seq_filter = True
        #TODO: this fails, i.e. curl -H 'apikey: ***REMOVED***' '***REMOVED***?
        filterName=sequence.submittedToInsdc&filterValue=true&filterWildcard=false'
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
        self.fixture.region_name_filter = 'ITS2-28S'
        self.assertEqual(4, self.fixture.get_number(),
                         'Expecting 4 consensus sequences filtered by '
                         'regionName = ITS2-28S, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        self.fixture.region_name_filter = '18s-28s'
        self.assertEqual(1, self.fixture.get_number(),
                         'Expecting 1 consensus sequences filtered by '
                         'regionName = 18s-28s, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        # projectName
        # Note: This test is failing
        '''
        self.fixture.project_name_filter = 'grdi'
        self.assertEqual(5555, self.fixture.get_number(), 
        'Expecting 5,555 consensus sequences filtered by 
        projectName = grdi, but got {}'.format(self.fixture.get_number()))
        self.fixture.clear_all_filters()
        '''

        # collection_code_filter
        self.fixture.collection_code_filter = 'lev'
        self.assertEqual(229, self.fixture.get_number(),
                         'Expecting 229 consensus sequences filtered by '
                         'collectionCode = lev, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        # taxonomy_rank_filter
        self.fixture.taxonomy_rank_filter = 'species'
        self.fixture.taxonomy_value_filter = 'megasperma'
        self.assertEqual(3, self.fixture.get_number(),
                         'Expecting 3 consensus sequences filtered by '
                         'taxonomy, but got {}'
                         .format(self.fixture.get_number()))
        self.fixture.clear_all_filters()

        # TODO: test combinations of filters

    def test_get_consensus_sequence_ids(self):
        self.fixture.specimen_num_filter = 4405
        actual = self.fixture.get_ids()
        self.assertTrue(actual, 'No Sequence IDs returned.')
        self.assertEqual(1, len(actual),
                         'Expecting 1 consensus sequence associated '
                         'with this specimen.')
        self.assertIn(358301, actual,
                      'Sequence ID 358301 is expected to be in results, '
                      'since it is consensus.')
        self.assertNotIn(27755, actual,
                         'Sequence ID 27755 is not expected to '
                         'be in results, since it is consensus.')

        self.fixture.specimen_num_filter = 4264
        actual = self.fixture.get_ids()
        self.assertTrue(actual, 'No Sequence IDs returned.')
        self.assertEqual(2, len(actual),
                         'Expecting 2 consensus sequences associated with '
                         'this specimen, but got {}'.format(len(actual)))
        self.assertIn(358302, actual,
                      'Sequence ID 358302 is expected to be in results, '
                      'since it is consensus')
        self.assertIn(4825628, actual,
                      'Sequence ID 4825628 is expected to be in results, '
                      'since it is consensus')
        self.assertNotIn(27755, actual,
                         'Sequence ID 27755 is not expected to be in '
                         'results, since it is consensus.')

    def test_create_get_delete_sequence(self):
        # create
        seq_id, err_cod, msg = self.fixture.create(name='Test',
                                                   sequence='ACGTCTGATCGATC')
        self.assertTrue(seq_id,
                        'Creating consensus sequence did not return an id.')
        self.assertEqual(err_cod, 201,
                         'Did not get successful exit code for '
                         'create consensus sequence.')

        # get
        self.fixture.sequence_name_filter = 'Test'
        seq_ids = self.fixture.get_ids()
        self.assertTrue(seq_ids,
                        'Creating consensus sequence did not return an id.')
        self.assertIn(seq_id, seq_ids,
                      'Expected sequence id was not in the result.')

        # delete
        delete_jsn_resp = self.fixture.delete(seq_id)
        self.assertEqual(200,
                         delete_jsn_resp['metadata']['statusCode'],
                         'Could not delete feature type.')

    def test_get_fasta_sequences_with_offset(self):
        self.fixture.specimen_num_filter = 4405
        actual_fasta, result_offset = self.fixture.get_fasta_sequences_with_offset(offset=0)
        self.assertTrue(actual_fasta, 'No Sequences returned.')
        self.assertIn('>seqdb|358301',
                      actual_fasta,
                      'Fasta return is expected to have sequence 3'
                      '58301, since it is consensus.')
        self.assertNotIn('>seqdb|27755',
                         actual_fasta,
                         'Fasta return is not expected to have '
                         'sequence id 27755, since it is raw.')

    def test_get_fasta_seq(self):
        actual = self.fixture.get_fasta_sequence('358301')
        self.assertTrue(actual, 'Fasta sequence is empty.')
        self.assertIn('>seqdb|358301', actual, 'Fasta does not contain >seqdb|358301.')


if __name__ == '__main__':
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

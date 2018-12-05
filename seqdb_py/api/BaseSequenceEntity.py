"""
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
"""

from .BaseApiEntity import BaseApiEntity
from .BaseSeqdbApi import UnexpectedContent
from .ProjectTagApi import ProjectTagApi


class BaseSequenceEntity(BaseApiEntity):

    def __init__(self, api_key, base_url, request_url):
        self.clear_all_filters()
        super(BaseSequenceEntity, self).__init__(
            api_key=api_key,
            base_url=base_url,
            request_url=request_url
        )

    @property
    def specimen_num_filter(self):
        return self.__specimen_num_filter

    @specimen_num_filter.setter
    def specimen_num_filter(self, specimen_num):
        self.__specimen_num_filter = specimen_num

    @property
    def sequence_name_filter(self):
        return self.__sequence_name_filter

    @sequence_name_filter.setter
    def sequence_name_filter(self, sequence_name):
        self.__sequence_name_filter = sequence_name

    @property
    def sample_name_filter(self):
        return self.__sample_name_filter

    @sample_name_filter.setter
    def sample_name_filter(self, sample_name):
        self.__sample_name_filter = sample_name

    @property
    def pub_ref_seq_filter(self):
        return self.__pub_ref_seq_filter

    @pub_ref_seq_filter.setter
    def pub_ref_seq_filter(self, is_submitted_to_insdc):
        if not (is_submitted_to_insdc) or isinstance(is_submitted_to_insdc, bool):
            self.__pub_ref_seq_filter = is_submitted_to_insdc
        else:
            pass
            # TODO handle exceptions

    @property
    def gen_bank_GI_filter(self):
        return self.__gen_bank_GI_filter

    @gen_bank_GI_filter.setter
    def gen_bank_GI_filter(self, gen_bank_GI):
        self.__gen_bank_GI_filter = gen_bank_GI

    @property
    def region_name_filter(self):
        return self.__region_name_filter

    @region_name_filter.setter
    def region_name_filter(self, region_name):
        self.__region_name_filter = region_name

    @property
    def project_name_filter(self):
        return self.__project_name_filter

    @project_name_filter.setter
    def project_name_filter(self, project_name):
        self.__project_name_filter = project_name

    @property
    def collection_code_filter(self):
        return self.__collection_code_filter

    @collection_code_filter.setter
    def collection_code_filter(self, collection_code):
        self.__collection_code_filter = collection_code

    @property
    def taxonomy_rank_filter(self):
        return self.__taxonomy_rank_filter

    @taxonomy_rank_filter.setter
    def taxonomy_rank_filter(self, taxonomy_rank):
        self.__taxonomy_rank_filter = taxonomy_rank

    @property
    def taxonomy_value_filter(self):
        return self.__taxonomy_value_filter

    @taxonomy_value_filter.setter
    def taxonomy_value_filter(self, taxonomy_value):
        self.__taxonomy_value_filter = taxonomy_value

    def get_param_str(self):
        params = ''

        if self.specimen_num_filter:
            params = params + 'filterName=specimen.number&filterValue={}' \
                              '&filterWildcard=false&' \
                .format(self.specimen_num_filter)

        if self.sequence_name_filter:
            params = params + 'filterName=sequence.name&filterValue={}' \
                              '&filterWildcard=true&' \
                .format(self.sequence_name_filter)

        if self.sample_name_filter:
            params = params + 'filterName=sample.name&filterValue={}' \
                              '&filterWildcard=true&' \
                .format(self.sample_name_filter)

        if self.pub_ref_seq_filter:
            params = params + 'filterName=sequence.submittedToInsdc' \
                              '&filterValue=true&filterWildcard=false&'

        if self.gen_bank_GI_filter:
            params = params + 'filterName=sequence.genBankGI&filterValue={}' \
                              '&filterWildcard=false&' \
                .format(self.gen_bank_GI_filter)

        if self.region_name_filter:
            params = params + 'filterName=region.name&filterValue={}' \
                              '&filterWildcard=true&' \
                .format(self.region_name_filter)

        if self.project_name_filter:
            project_tag = ProjectTagApi(
                api_key=self.api_key, base_url=self.base_url)
            project_tag.name_filter = self.project_name_filter
            project_ids = project_tag.get_ids()

            if len(project_ids) == 1:
                params = params + 'tagId=%s&' % project_ids[0]

            elif len(project_ids) > 1:
                raise Exception('More than one project is associated with the '
                                'name "%s". Currently sequences can only be '
                                'filtered on one project only. '
                                'Please refine your search.')

        if self.collection_code_filter:
            params = params + 'filterName=biologicalCollection.name&' \
                              'filterValue={}&filterWildcard=false&' \
                .format(self.collection_code_filter)

        if self.taxonomy_rank_filter and self.taxonomy_value_filter:
            filter_name = self.convert_ncbi_to_seqdb_tax_rank(
                self.taxonomy_rank_filter)
            params = params + 'filterName=specimen.identification.taxonomy.' \
                              '{}&filterValue={}&filterOperator=and' \
                              '&filterWildcard=true&' \
                .format(filter_name, self.taxonomy_value_filter)

        return params

    def clear_all_filters(self):
        self.specimen_num_filter = None
        self.sequence_name_filter = None
        self.sample_name_filter = None
        self.pub_ref_seq_filter = None
        self.gen_bank_GI_filter = None
        self.region_name_filter = None
        self.project_name_filter = None
        self.collection_code_filter = None
        self.taxonomy_rank_filter = None
        self.taxonomy_value_filter = None

    def get_fasta_sequences_with_offset(self, offset, limit=20):
        fasta_resp, result_offset = self.get_sequences_with_offset(
            sequence_format='fasta', offset=offset, limit=limit)
        if fasta_resp and fasta_resp[0] != '>':
            raise UnexpectedContent('Response is not in fasta format.')

        return fasta_resp, result_offset

    def get_sequences_with_offset(self, sequence_format, offset, limit=20):

        if offset < 0:
            raise Exception('Negative offset: either you have retrieved all '
                            'sequences or the method usage is incorrect.')

        if sequence_format not in {'fasta', 'fastq'}:
            raise SystemError('Sequences can only be in fasta '
                              'or fastq format.')

        # TODO: refactor this method so that the api call to get number of sequences
        #    is not called each time this method is called
        # Define a global var for total number of sequences to avoid calling the 
        # api method every time __getSequencesWithOffset is called
        '''
        global globalSequenceNumber
        if not globalSequenceNumber:
            globalSequenceNumber = self.get_number()
        '''

        params = self.get_param_str() + 'limit={}&offset={}&' \
            .format(limit, offset)

        resp = self.retrieve('{}{}.{}'.format(self.base_url, self.request_url,
                                              sequence_format), params=params)

        sequence_formatted_resp = resp.content

        # This whole weirdness with getting the number of results from response, instead of just moving
        # one frame forward (i.e. move by a number of sequences that is specified in 'limit') is because
        # SeqDB has a weird bug where the rest of the frame response is not returned if there is a problem 
        # with a data record, i.e. it was not deleted properly and quering for its fasta sequence returns 
        # no results.
        curr_result_num = limit
        if sequence_formatted_resp:
            if sequence_format == 'fasta':
                curr_result_num = sequence_formatted_resp.count('>')
            elif sequence_format == 'fastq':
                curr_result_num = sequence_formatted_resp.count('@seqdb')
        else:
            curr_result_num = 1

        # Calculate the new offset
        result_offset = -1
        # TODO: refactor for efficiency, see TODO above.
        result_total_num = self.get_number()
        result_returned_so_far = curr_result_num + offset
        if result_total_num > result_returned_so_far:
            result_offset = result_returned_so_far

        # print 'Offset: {}\tReturned results: {}\t'
        # .format(offset, sequence_formatted_resp.count('>'))

        return sequence_formatted_resp, result_offset

    def get_fasta_sequence(self, seq_id):
        """ Returns a sequence in a fasta format.
            Note that to get multiple fasta sequences it is bests to use 
            getFastaSequencesWithOffset() method, to avoid multiple call
            to the SeqDB API.
        """
        return self._get_formatted_seq(seq_id, 'fasta')

    def _get_formatted_seq(self, seq_id, sequence_format):
        """ Gets sequence in fasta or fastq format 
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
        """
        if sequence_format not in {'fasta', 'fastq'}:
            raise SystemError('Sequences can only be in fasta or fastq format.')

        url = '{}{}/{}.{}'.format(self.base_url, self.request_url, str(seq_id), sequence_format)
        response = self.retrieve(url)
        return response.content

    ########################################################################
    #                     Taxonomy Helpers
    ########################################################################
    # Used in this entity, but should they be in DeterminationApi?

    # This is a conversion dictionary for NCBI taxonomy to be imported to SeqDB. 
    # Keys are seqdb taxonomy ranks, and values are NCBI
    _seqdb_to_ncbi_taxonomy = {'kingdom': 'kingdom',
                               'phylum': 'phylum',
                               'taxanomicClass': 'class',
                               'taxanomicOrder': 'order',
                               'family': 'family',
                               'genus': 'genus',
                               'subgenera': 'subgenus',
                               'species': 'species',
                               'variety': 'varietas',
                               }

    # NOTE that following SeqDB taxonomic ranks currently do not have equivalent in NCBI taxonomy:
    # 'strain','authors','division','synonym','typeSpecimen','taxanomicGroup','commonName','state'

    def vet_taxonomy(self, taxonomy):
        """ Checks taxonomy dictionary and returns it in the format that can be accepted by 
            seqDB (i.e. appropriate ranks)
        """
        vetted_taxonomy = dict()
        if taxonomy:
            for seqdb_rank in self._seqdb_to_ncbi_taxonomy:
                ncbi_rank = self._seqdb_to_ncbi_taxonomy[seqdb_rank]
                if ncbi_rank in taxonomy:
                    vetted_taxonomy[seqdb_rank] = taxonomy[ncbi_rank]

        return vetted_taxonomy

    def convert_ncbi_to_seqdb_tax_rank(self, tax_rank):
        for seqdb_tax_rank in self._seqdb_to_ncbi_taxonomy:
            if self._seqdb_to_ncbi_taxonomy[seqdb_tax_rank] == tax_rank:
                return seqdb_tax_rank
        return 'No matches found.'

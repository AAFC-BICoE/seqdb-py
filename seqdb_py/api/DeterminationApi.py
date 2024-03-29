"""
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
"""

import json

from api.BaseApiEntity import BaseApiEntity
from api.BaseSeqdbApi import UnexpectedContent


class DeterminationApi(BaseApiEntity):

    def __init__(self, api_key, base_url):
        super(DeterminationApi, self).__init__(api_key=api_key,
                                               base_url=base_url,
                                               request_url='determination')

    def get_param_str(self):
        return ''
    
    # Curerntly accepted taxonomy fields in SeqDB:
    # 'notes', 'typeSpecimen', 'phylum', 'subgenera', 'variety', 'kingdom', 'division', 'genus', 
    # 'taxanomicGroup', 'state', 'strain', 'commonName', 'id', 'species', 'taxanomicOrder', 
    # 'superfamily', 'family', 'authors', 'synonym', 'lastModified', 'taxanomicClass'
    ###
    
    def create_sequence_determination(self, sequence_id, taxonomy,
                                      is_accepted=False, ncbi_taxon_id=None,
                                      notes=None):

        """ Creates a determination for a sequence
        Args:
            sequence_id: id of a sequence for which determination is writtemn
            is_accepted: boolean, whether this determination is an accepted determination
            ncbi_taxon_id: integer of a ncbi taxon id
            taxonomy:  list of tuples (taxonomy rank, taxonomy rank name). Usually from NCBI.
                 i.e. [['species', 'Phytophthora lateralis'], ['genus', 'Phytophthora']]
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        """
        # taxonomy = self.vet_taxonomy(taxonomy)
        post_data = {
            'identification': {
                'accepted': is_accepted,
                'sequence': {'id': sequence_id},
                'taxonomy': taxonomy,
                'ncbiTaxonId': ncbi_taxon_id,
                }
            }

        resp = super(DeterminationApi, self)\
            .create('{}{}'
                    .format(self.base_url, self.request_url),
                    json.dumps(post_data))
        jsn_resp = resp.json()
        
        if 'result' and 'metadata' not in jsn_resp:
            raise UnexpectedContent(response=jsn_resp)
        
        if 'statusCode' and 'message' not in jsn_resp['metadata']:
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result']

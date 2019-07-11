"""
Created on Feb 16, 2016

@author: Oksana Korol

Class that extracts common functionality for all SeqDB API entities
"""

# import pdb
from abc import abstractmethod
# from wsgiref.util import request_uri

from BaseSeqdbApi import BaseSeqdbApi, UnexpectedContent


class BaseApiEntity(BaseSeqdbApi):

    def __init__(self, api_key, base_url, request_url):
        """ Initializes the object with API access attributes and
            specific API call (request url)
        Args:
            api_key: user api key for accessing SeqDB
            base_url: base SeqDB API url (ex.***REMOVED***/api/v1/)
            request_url: specific entity request param that will be added at
            the end of base_url
                    (ex. 'sequence')
        """
        super(BaseApiEntity, self).__init__(api_key, base_url)
        self.request_url = request_url

    @abstractmethod
    def get_param_str(self):
        """ Based on the object specified filter parameters, create a parameter
            string, which will be added to the request
        """
        raise NotImplementedError('Must override getParamStr')

    def get_entity(self, entity_id):
        """ Retrieves an entity
        Args:
            entity_id: entity id
        Returns:
            'result' from the json response OR nothing if entity was not found
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        """
        new_request_url = self.request_url + '/' + str(entity_id)
        jsn_resp = self.retrieve_json(request_url=new_request_url)

        if jsn_resp:
            return jsn_resp['result']
        else:
            return ''

    def get_number(self):
        """ Returns a number of entities
        """
        param_str = self.get_param_str()
        # jsn_resp = self.retrieve_json(request_url=request_uri, params=param_str)
        jsn_resp = self.retrieve_json(request_url=self.request_url, params=param_str)
        result_num = int(jsn_resp['metadata']['resultCount'])

        return result_num

    def delete(self, entity_id):
        """ Deletes an entity from SeqDB
        Args:
            entity_id: id of the entity to be deleted
        Returns:
            json response
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.HTTPError
            UnexpectedContent
        """
        request_url = self.request_url + '/' + str(entity_id)
        jsn_resp = super(BaseApiEntity, self)\
            .delete(self.base_url + request_url).json()

        if 'statusCode' and 'message' not in jsn_resp['metadata']:
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp

    def bulk_delete(self, entity_ids):
        """ Deletes a list of entities from SeqDB
        Args:
            entity_ids: list of seqdb entity ids to be deleted
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        """
        for entityId in entity_ids:
            self.delete(entityId)

    def param_str(self, offset=0, limit=0):
        """ Get entity IDs with offset, filtered by specified filter parameters
        Args:
            offset: nothing if it is a first query, then number of records
            from which to load the next set of ids
        Returns:
            a list of entity ids
            offset of results. If 0 then all/last set of results have been
            retrieved, if > 0, then the function has to be called again with
            this offset to retrieve more results
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        """
        param_str = self.get_param_str()
        jsn_resp, result_offset = self.retrieve_json_with_offset(
            request_url=self.request_url,
            params=param_str,
            offset=offset,
            limit=limit
        )

        entity_ids = ''

        if jsn_resp:
            entity_ids = jsn_resp['result']

        return entity_ids, result_offset

    def get_ids(self):
        """
        Returns: all entity ids, that correspond to the set filters

        Companion method to getProjectTagWithOffset. Returns all the results, iterating with offset.
        """
        result_limit = 50
        tag_ids, offset = self.get_ids_with_offset(offset=0, limit=result_limit)

        while offset:
            curr_tag_ids, offset = self.get_ids_with_offset(offset=offset, limit=result_limit)
            tag_ids.extend(curr_tag_ids)

        return tag_ids

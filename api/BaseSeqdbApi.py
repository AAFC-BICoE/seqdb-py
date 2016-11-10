'''
Created on Feb 16, 2016

@author: Oksana Korol

Sequence DB Web Services module. Define basic API calls.
'''
import json
import logging
import requests
import pdb

class UnexpectedContent(requests.exceptions.RequestException):
    error_msg = "Unexpected content of the Web Services response: "

    def __init__(self, response):
        self.response = response
        logging.error(self.error_msg + str(response))

    def __str__(self):
        return repr(self.error_msg + str(self.response))


# CRUD (create, retrieve, update, delete)

class BaseSeqdbApi(object):

    def __init__(self, api_key, base_url):
        self.api_key = api_key
        if not base_url.endswith("/"):
            base_url = base_url + "/"
        self.base_url = base_url


    def retrieve(self, request_url, params=None):
        ''' Submits a request to SeqDB web services
        Kwargs:
            request_url: full request url
        Returns:
            SeqDB API response (still need to format it to JSon, etc.) or
            nothing if resource was not found
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            requests.exceptions.InvalidSchema
        '''
        req_header = { 'apikey': self.api_key }
        
        resp = requests.get(request_url, headers=req_header, params=params)
        logging.debug("Request: {} {} {}".format(resp.request.method, resp.request.url, resp.request.body))
        logging.debug("Response Status Code: {}".format(resp.status_code))
        #resp.content: str {"count":288,"limit":20,"message":"Query completed successfully","offset":0,"result":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],"sortColumn":"regionId","sortOrder":1,"statusCode":200}
        
        if resp.status_code == 404:
            resp = ''
        else:
            # Will raise HTTPError exception if response status was not ok
            try:
                resp.raise_for_status()
            except requests.exceptions.HTTPError as e:
                error_msg = "Retrieve failed. Request body: \n {} \nRequest URL: {}.\nError message: {}.".format(e.request.body, e.request.url, e.response.text)
                e.message = e.message + error_msg
                logging.error(error_msg)
                raise e
      
        return resp


    def retrieveJson(self, request_url, params=None):
        ''' Submits a request to SeqDB web services
            Note: in case of paginated results, will return one page only. 
                    Use WithOffset to get all results.
        Args:
            request_url: part of the API url, specific to the requests (i.e. "/sequences")
        Returns:
            json formatted object
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        resp = self.retrieve("{}{}".format(self.base_url, request_url), params=params)
        if resp:
            jsn_resp =  json.loads(resp.text)
            resp = jsn_resp
            
            if 'result' not in jsn_resp:
                raise UnexpectedContent(response=jsn_resp)
                    
        return resp
    
    
    def retrieveJsonWithOffset(self, request_url, params=None, offset=0, limit=0):
        ''' Submits a request to SeqDB web services with offset to handle paginated results
        Args:
            request_url: part of the API url, specific to the requests (i.e. "/sequences")
            params: str formatted parameters to add to a request (filters and such)
            offset: 0 for the firsst query, number of results returned for the subsequent queries
        Returns:
            json formatted object 
            offset of results. If 0 then all/last set of results have been retrieved, if >= 0,
                then the function has to be called again with this offset to retrieve more results.
                When the request is exhausted, offset will be -1
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        
        ### NOTE: the params need to be in the string format, because dictionary does not allow duplicate keys.
        ###    This restriction will prevent queries with multiple filter names, i.e. when trying to filter
        ###    sequences by specimen number and sequence name keyword: filterName="sequence.name" and 
        ###    filterName="specimen.number"
         
        if params and type(params) is not str:
            raise "BaseSeqdbApi.retrieveJsonWithOffset: Parameters to api have to be in the str format"
        
        if params and "offset=" in params:
            raise "BaseSeqdbApi.retrieveJsonWithOffset: Parameters to api should not contain 'offset=' in this method"

        if offset:
            params ="{}&offset={}&".format(params,offset)
        
        if limit:
            params ="{}&limit={}&".format(params,limit)
                       
        jsn_resp = self.retrieveJson(request_url, params)
        
        result_offset = 0
        
        if jsn_resp:
            if 'metadata' not in jsn_resp:
                raise UnexpectedContent(response=jsn_resp)
                
            if 'resultCount' not in jsn_resp['metadata'] and 'limit' not in jsn_resp['metadata'] and 'offset' not in jsn_resp['metadata']:
                raise UnexpectedContent(response=jsn_resp)
            
            # TODO verify this works for all cases
            #if jsn_resp['metadata']['resultCount'] > 0 and not jsn_resp['result']:
            #    raise UnexpectedContent(response=jsn_resp)

            # Checking for paginated results
            result_total_num = int(jsn_resp['metadata']['resultCount'])
            result_returned_so_far =  int(jsn_resp['metadata']['limit']) + int(jsn_resp['metadata']['offset'])
            if (result_total_num > result_returned_so_far):    
                result_offset = result_returned_so_far
            
        return jsn_resp, result_offset
       
    
    def update(self, request_url, json_data):
        ''' Updates a SeqDB entity
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
        '''
        json_data = str(json_data)
        json_data = json_data.replace("u'", "\"")
        json_data = json_data.replace("'","\"")
        req_header = {
            'apikey': self.api_key, 'Content-Type': 'application/json'}
        resp = requests.put(self.base_url + request_url, headers=req_header, data=json_data)
        logging.debug("Request: {} {} {}".format(resp.request.method, resp.request.url, resp.request.body))
        logging.debug("Response Status Code: {}".format(resp.status_code))

        try:
            resp.raise_for_status()
        except requests.exceptions.HTTPError as e:
            error_msg = "Update failed. Request body: \n {} \nRequest URL: {}.\nError message: {}.".format(e.request.body, e.request.url, e.response.text)
            e.message = e.message + error_msg
            logging.error(error_msg)
            raise e

        return resp


    def create(self, request_url, json_data):
        ''' Creates a SeqDB entity
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
        '''
        req_header = {
            'apikey': self.api_key, 'Content-Type': 'application/json'}
        # (url, data, json)
        resp = requests.post(request_url, headers=req_header, data=json_data)
        logging.debug("Request: {} {} {}".format(resp.request.method, resp.request.url, resp.request.body))
        logging.debug("Response Status Code: {}".format(resp.status_code))

        try:
            resp.raise_for_status()
        except requests.exceptions.HTTPError as e:
            error_msg = "Create failed. Request body: \n {} \nRequest URL: {}.\nError message: {}.".format(e.request.body, e.request.url, e.response.text)
            e.message = e.message + error_msg
            logging.error(error_msg)
            if resp.status_code == 401:
                logging.error("Not sufficient permissions to write information to SeqDB. Please contact SeqDB administrator.")
            raise e
        
        return resp


    def delete(self, request_url):
        ''' Creates a SeqDB entity
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
        '''

        req_header = {'apikey': self.api_key}
        resp = requests.delete(request_url, headers=req_header)
        logging.debug("Request: {} {} {}".format(resp.request.method, resp.request.url, resp.request.body))
        logging.debug("Response Status Code: {}".format(resp.status_code))

        # Will raise HTTPError exception if response status was not ok
        try:
            resp.raise_for_status()
        except requests.exceptions.HTTPError as e:
            error_msg = "Delete failed. Request body: \n {} \nRequest URL: {}.\nError message: {}.".format(e.request.body, e.request.url, e.response.text)
            e.message = e.message + error_msg
            logging.error(error_msg)
            if resp.status_code == 401:
                logging.error("Not sufficient permissions to delete information from SeqDB. Please contact SeqDB administrator.")
            raise e

        return resp  
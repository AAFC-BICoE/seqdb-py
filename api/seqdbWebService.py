'''
Created on Mar 4, 2015

@author: korolo

Sequence DB Web Services module
'''

import base64
import gzip
import json
import logging
import mimetypes
import os

import requests


class UnexpectedContent(requests.exceptions.RequestException):
    error_msg = "Unexpected content of the Web Services response: "

    def __init__(self, response):
        self.response = response
        logging.error(self.error_msg + str(response))

    def __str__(self):
        return repr(self.error_msg + str(self.response))


# CRUD (create, retrieve, update, delete)

class seqdbWebService:

    def __init__(self, api_key, base_url):
        self.api_key = api_key
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
        #resp.content: str {"count":288,"limit":20,"message":"Query completed successfully","offset":0,"result":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],"sortColumn":"regionId","sortOrder":1,"statusCode":200}
        
        if resp.status_code == 404:
            resp = ''
        else:
            # Will raise HTTPError exception if response status was not ok
            try:
                resp.raise_for_status()
            except requests.exceptions.HTTPError as e:
                error_msg = "Retrieve failed. Request body: \n %s \nRequest URL: %s.\nError message: %s." % (e.request.body, e.request.url, e.response.text)
                e.message = e.message + error_msg
                logging.error(error_msg)
                raise e
      

        return resp


    def retrieveJson(self, request_url, params=None):
        ''' Submits a request to SeqDB web services
        Returns:
            json formatted object
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        resp = self.retrieve(request_url, params=params)
        if resp:
            # See if results were paginated, and if yes - retrieve the rest of the results
            jsn_resp =  json.loads(resp.text)
            resp = jsn_resp
            
            #if 'count' and 'sortColum' and 'limit' and 'offset' and 'message' and 'statusCode' and 'sortOrder' not in jsn_resp:
            if 'result' not in jsn_resp.keys():
                raise UnexpectedContent(response=jsn_resp)
            
            '''
            # If these values are present, there is a chance that we can get paginated result
            if 'count' in jsn_resp.keys() and 'limit' in jsn_resp.keys():

                result_num = int(jsn_resp['count'])
                result_per_page =  int(jsn_resp['limit'])
                
                if (result_num > result_per_page):
                    
                    for curr_offset in xrange(result_per_page, result_num, result_per_page):
                        if params:
                            if type(params) is dict:
                                params['offset'] = curr_offset
                            elif type(params) is str:
                                if 'offset' in params:
                                    params = re.sub('offset=[0-9]*', "offset=%s" %curr_offset, params)
                                else:
                                    params = "%s&offset=%s" %(params,curr_offset)
                        else:
                            params={'offset': curr_offset}
                            
                        resp2 = self.retrieve(request_url, params)
                        jsn_resp2 =  json.loads(resp2.text)
                
                        if 'count' and 'limit' and 'result' not in jsn_resp2.keys():
                            raise UnexpectedContent(response=jsn_resp)
                        
                        jsn_resp['result'] = jsn_resp['result'] + jsn_resp2['result'] 
                        
                    jsn_resp['limit'] = jsn_resp['count']
        
            resp = jsn_resp
            '''
        
        return resp
    
    def retrieveJsonWithOffset(self, request_url, params=None, offset=0):
        ''' Submits a request to SeqDB web services with offset to handle paginated results
        Args:
            request_url: part of the API url, specific to the requests (i.e. "/sequences")
            params: str formatted parameters to add to a request (filters and such)
            offset: 0 for the firsst query, number of results returned for the subsequent queries
        Returns:
            json formatted object 
            offset of results. If 0 then all/last set of results have been retrieved, if > 0,
                then the function has to be called again with this offset to retrieve more results
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
            raise "seqdbWebServices.retrieveJsonWithOffset: Parameters to api have to be in the str format"
        
        if params and "offset=" in params:
            raise "seqdbWebServices.retrieveJsonWithOffset: Parameters to api should not contain 'offset=' in this method"

        if offset:
            params ="%s&offset=%s&" %(params,offset)
                       
        jsn_resp = self.retrieveJson(self.base_url + request_url, params)
        
        result_offset = 0
        
        if jsn_resp:
            if 'count' not in jsn_resp and 'limit' not in jsn_resp and 'offset' not in jsn_resp:
                raise UnexpectedContent(response=jsn_resp)
            
            # TODO verify this works for all cases
            #if jsn_resp['count'] > 0 and not jsn_resp['result']:
            #    raise UnexpectedContent(response=jsn_resp)

            # Checking for paginated results
            result_total_num = int(jsn_resp['count'])
            result_returned_so_far =  int(jsn_resp['limit']) + int(jsn_resp['offset'])
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
        req_header = {
            'apikey': self.api_key, 'Content-Type': 'application/json'}
        resp = requests.put(request_url, headers=req_header, data=json_data)

        try:
            resp.raise_for_status()
        except requests.exceptions.HTTPError as e:
            error_msg = "Update failed. Request body: \n %s \nRequest URL: %s.\nError message: %s." % (e.request.body, e.request.url, e.response.text)
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

        try:
            resp.raise_for_status()
        except requests.exceptions.HTTPError as e:
            error_msg = "Create failed. Request body: \n %s \nRequest URL: %s.\nError message: %s." % (e.request.body, e.request.url, e.response.text)
            e.message = e.message + error_msg
            logging.error(error_msg)
            #print error_msg
            #raise requests.exceptions.HTTPError(error_msg)  # this will lose the stacktrace
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

        # Will raise HTTPError exception if response status was not ok
        try:
            resp.raise_for_status()
        except requests.exceptions.HTTPError as e:
            error_msg = "Delete failed. Request body: \n %s \nRequest URL: %s.\nError message: %s." % (e.request.body, e.request.url, e.response.text)
            e.message = e.message + error_msg
            logging.error(error_msg)
            raise e
        

        return resp


    def getSequenceIdsWithOffset(self, specimenNum=None, sequenceName=None, pubRefSeq=None, genBankGI=None, offset=0):
        ''' Returns sequence ids, limited by the specified filter parameters
        Agrs:
            specimenNum: specimen number (identifier) for which sequence IDs will be retrieved
            sequenceName: keyword in the sequence name (i.e. not a direct match)
            pubRefSeq: whether the sequence is a public reference sequence
            offset: nothing if it is a first query, then number of records from which to load the next set of ids
        Returns:
            a list of seqdb sequence ids 
            offset of results. If 0 then all/last set of results have been retrieved, if > 0,
                then the function has to be called again with this offset to retrieve more results
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        
        params = ""
             
        if specimenNum:
            params = "%sfilterName=specimen.number&filterValue=%s&filterOperator=and&filterWildcard=false&" %(params,specimenNum)
    
        if sequenceName:
            params = "%sfilterName=sequence.name&filterValue=%s&filterOperator=and&filterWildcard=true&" %(params,sequenceName)
            
        if pubRefSeq:
            params = params + "filterName=sequence.submittedToInsdc&filterValue=true&filterOperator=and&filterWildcard=true"
        
        if genBankGI:
            params = params + "filterName=sequence.genBankGI&filterValue=%s&filterOperator=and&filterWildcard=false&" %genBankGI
        
        jsn_resp, result_offset = self.retrieveJsonWithOffset(request_url="/sequence", params=params, offset=offset)
        
        sequence_ids = ""
        
        if jsn_resp:            
            sequence_ids = jsn_resp['result']
        
        if genBankGI and len(sequence_ids) > 1:
            raise SystemError(
                "More than one record associated with GenBank GI, which should be unique.")
        
        return sequence_ids, result_offset
       

   
    def getSequenceIdsByRegionWithOffset(self, region_id, offset=0):
        ''' Given a region id, return sequence ids, belonging to this region
        Returns:
            a list of seqdb sequence ids for the created sequences OR empty
            list id created failed
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        
        request_url = "/region/" + str(region_id) + "/sequence"
        jsn_resp, result_offset = self.retrieveJsonWithOffset(request_url=request_url, offset=offset)
        
        sequence_ids = ""
        
        if jsn_resp:            
            sequence_ids = jsn_resp['result']            
        
        return sequence_ids, result_offset



    def getJsonConsensusSequenceIds(self, params=None):
        ''' Gets Json api response object with sequence ids of all consensus sequences.
            Note: be sure to handle pagination, since first api query may not return all 
            possible results.
        Args:
            sequenceName: A name to identify the sequence
            genBankGI: GenBankGI if the sequence is public reference
        Returns:
            Json object with the sequence ids
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        jsn_resp = self.retrieveJson(
            self.base_url + "/consensus", params=params)

        
        if jsn_resp['count'] > 0 and not jsn_resp['result']:
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp


    def getConsensusSequenceIdsWithOffset(self, specimenNum=None, sequenceName=None, pubRefSeq=None, genBankGI=None, offset=0):
        ''' Returns sequence ids, limited by the specified filter parameters
        Agrs:
            specimenNum: specimen number (identifier) for which sequence IDs will be retrieved
            sequenceName: keyword in the sequence name (i.e. not a direct match)
            pubRefSeq: whether the sequence is a public reference sequence
            genBankGI: genBank GI (identifier) for those sequences that are in genBank
            offset: nothing if it is a first query, then number of records from which to load the next set of ids
        Returns:
            a list of seqdb sequence ids 
            offset of results. If 0 then all/last set of results have been retrieved, if > 0,
                then the function has to be called again with this offset to retrieve more results
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        params = ''
        if specimenNum:
            params = params + "filterName=specimen.number&filterValue=%s&filterOperator=and&filterWildcard=false&" %specimenNum
            
        if sequenceName:
            params = params + "filterName=sequence.name&filterValue=%s&filterOperator=and&filterWildcard=true&" %sequenceName
        
        if pubRefSeq:
            params = params + "filterName=sequence.submittedToInsdc&filterValue=true&filterOperator=and&filterWildcard=true"
            
        if genBankGI:
            params = params + "filterName=sequence.genBankGI&filterValue=%s&filterOperator=and&filterWildcard=false&" %genBankGI
        
        #jsn_resp = self.getJsonConsensusSequenceIds(params)
        #jsn_resp = self.retrieveJson(self.base_url + "/consensus", params=params)
        jsn_resp, result_offset = self.retrieveJsonWithOffset(request_url="/consensus", params=params, offset=offset)
        
        
        sequence_ids = ""
        
        if jsn_resp:  
            sequence_ids = jsn_resp['result']
            
        if genBankGI and len(sequence_ids) > 1:
            raise SystemError(
                "More than one record associated with GenBank GI, which should be unique.")
            
        
        return sequence_ids, result_offset


    # TODO verify which payload values are required by the SeqDB WS API
    def createConsensusSequence(
            self, name, sequence, qualities=None,
            seqType="N", readNum=0, additional=None):
        post_data = {
            'consensus': {
                'name': name,
                'seq': sequence,
                'seqType': seqType,
                'readNum': readNum,
                'qualities': qualities,
            },
        }

        if additional is not None:
            post_data['consensus'].update(additional)

        resp = self.create(self.base_url + "/consensus", json.dumps(post_data))
        jsn_resp = resp.json()

        if 'result' and 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result'], jsn_resp['statusCode'], jsn_resp['message']


    def deleteConsensusSequence(self, consensus_id):
        request_url = "/consensus/" + str(consensus_id)
        jsn_resp = self.delete(self.base_url + request_url).json()

        if 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp

    def importChromatSequences(
            self, blob, dest_file_name,
            notes="", trace_file_path=""):
        ''' Imports a binary blob (i.e. chromatogram) to seqdb
        Kwargs:
            dest_file_name name with which the chromatogram will be saved on
                           sedDB side; i.e. expected to have .ab1 extension
        Returns:
            a list of seqdb sequence ids for the created sequences OR empty
            list if creation failed
        Raises:
            UnexpectedContent
        '''

        chromat_b64 = base64.b64encode(blob)

        post_data = {
            "sequenceImportPayload": {
                "base64File": chromat_b64,
                "fileName": dest_file_name,
                "plateType": 1,
                "createLocation": False,
                "traceFilePath": trace_file_path,
                "notes": notes

            }
        }

        resp = self.create(
            self.base_url + '/sequenceImport', json.dumps(post_data))
        jsn_resp = resp.json()

        if 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        result = ""
        if 'result' in jsn_resp.keys():
            result = jsn_resp['result']
        else:
            logging.warn(
                "Importing chromatogram failed with status code: '%s' "
                "and message: '%s'" % (
                    jsn_resp['statusCode'], jsn_resp['message']))

        return result

    def importChromatSequencesFromFile(
            self, chromat_file, notes="",
            trace_file_path="", dest_file_name=""):
        ''' Imports a chromatogram from a file
        Args:
            chromat_file name of the chromatogram file
        Returns:
            list of sequence ids of the sequences that were imported from the
            chromatogram
        Raises:
            IOError
            UnexpectedContent
        '''

        if not os.path.isfile(chromat_file):
            raise IOError("Expecting a file, but got a directory.")

        chromat_file_name = os.path.basename(chromat_file)

        if not dest_file_name:
            # get just the file name, no extension
            dest_file_name = os.path.splitext(
                os.path.basename(chromat_file_name))[0]

        if mimetypes.guess_type(chromat_file_name)[1] == "gzip":

            file_strem = gzip.open(chromat_file, "rb")

            # Remove .gz extension for file name
            chromat_file_name = chromat_file_name[:-3]

        else:
            file_strem = open(chromat_file, "r")

        blob = file_strem.read()

        file_strem.close()

        return self.importChromatSequences(
                blob=blob, dest_file_name=dest_file_name,
                notes=notes, trace_file_path=trace_file_path)

    def deleteSequence(self, seq_id):
        ''' Deletes a SeqDB sequence
        Args:
            seq_id: id of the sequence to be deleted
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        request_url = "sequence/" + str(seq_id)
        jsn_resp = self.delete(self.base_url + request_url).json()

        if 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp

    def bulkDeleteSequence(self, seq_ids):
        ''' Deletes a list of SeqDB sequences
        Args:
            seq_ids: list of seqdb sequence ids to be deleted
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        for seq_id in seq_ids:
            self.deleteSequence(seq_id)

    def updateSeqSource(self, sequenceId, params):
        resp = self.update(
            self.base_url + "/sequence/" + str(sequenceId) + "/seqSource",
            json.dumps(params))
        jsn_resp = resp.json()

        if 'result' and 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result'], jsn_resp['statusCode'], jsn_resp['message']


    def getJsonSeqSource(self, sequenceId):
        jsn_resp = self.retrieveJson(
            self.base_url + "/sequence/" + str(sequenceId) + "/seqSource")

        return jsn_resp


    def getJsonSequence(self, seq_id, params=None):
        jsn_resp = self.retrieveJson(
            self.base_url + "/sequence/" + str(seq_id), params)
        if not jsn_resp['result']:
            raise UnexpectedContent(response=jsn_resp)

        jsn_seq = jsn_resp['result']
        if 'seq' not in jsn_seq.keys() or 'name' not in jsn_seq.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_seq

    def getFastaSeqPlus(self, seq_id):
        ''' Gets sequence from SeqDB and returns it in a fasta format with a
            unique header. Note, the fasta formatting is done here, instead of
            using seqdb fasta web service request
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        jsn_seq = self.getJsonSequence(seq_id)
        # TODO Use BioPython to format
        fasta_seq = '>' + \
            jsn_seq['name'] + '|seqdbId:' + \
            str(seq_id) + '\n' + jsn_seq['seq'] + '\n'
        return fasta_seq

    def getFastaSeq(self, seq_id):
        ''' Gets sequence in fasta format (SeqDB fasta request)
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
        '''
        url = self.base_url + "/sequence/" + str(seq_id) + ".fasta"
        response = self.retrieve(url)
        return response.content



    def getItsRegionIds(self, offset=0):
        ''' Get region IDs of ITS sequences
        Args:
            offset: 0 if it is a first query, then number of records from which to load the next set of ids
        Returns:
            a list of seqdb its region ids 
            offset of results. If 0 then all/last set of results have been retrieved, if > 0,
                then the function has to be called again with this offset to retrieve more results
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        
        its_region_names = ["ssu", "16s", "18s", "its", "5.8s", "lsu", "23s", "25s", "28s", "internal transcribed spacer"]
        its_region_ids = set()
        
        for its_name in its_region_names:
            current_ids, offset = self.getRegionIdsWithOffset(regionName=its_name, offset=offset)
            its_region_ids.update(current_ids)
            while offset:
                current_ids, offset = self.getRegionIdsWithOffset(regionName=its_name, offset=offset)
                its_region_ids.update(current_ids)       
            
        return list(its_region_ids)
        
    '''
    def getRegionIdsByName(self, name, offset=0):
        params = {
            'filterName': 'name',
            'filterValue': name,
            'filterOperator': 'and',
            'filterWildcard': 'true'
        }
        return self.getRegionIdsWithOffset(params, offset)
    '''

    def getRegionIdsWithOffset(self, regionName=None, offset=0):
        ''' Get region IDs
        Args:
            regionName: gene region name to filter results by 
            offset: nothing if it is a first query, then number of records from which to load the next set of ids
        Returns:
            a list of seqdb region ids 
            offset of results. If 0 then all/last set of results have been retrieved, if > 0,
                then the function has to be called again with this offset to retrieve more results
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        params = ""
             
        if regionName:
            params = "%sfilterName=name&filterValue=%s&filterOperator=and&filterWildcard=true&" %(params,regionName)
    
        jsn_resp, result_offset = self.retrieveJsonWithOffset(request_url="/region", params=params, offset=offset)
        
        region_ids = ""
        
        if jsn_resp:
            region_ids = jsn_resp['result']
        
        return region_ids, result_offset
        

    def getRegionName(self, region_id):
        ''' Given a region id, returns region name or empty string if no response
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        jsn_resp = self.retrieveJson(self.base_url + "/region/" + str(region_id))
        if jsn_resp:
            return jsn_resp['result']['name']
        else:
            return ''
        

    def createRegion(self, name, description, group_id=1):
        ''' Creates a region
        Args:
            name: region name
            description: region description
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        # TODO Question requirement for gene region to be associated with a
        # group
        # TODO Don't hard code group
        post_data = {
            "region": {
                "description": description,
                "group": {"id": group_id},
                "name": name}
            }

        resp = self.create(self.base_url + '/region', json.dumps(post_data))
        jsn_resp = resp.json()

        if 'result' and 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result']

    def deleteRegion(self, regionId):
        ''' Deletes a region
        Args:
            regionId: id of a region to be created
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        request_url = "/region/" + str(regionId)
        jsn_resp = self.delete(self.base_url + request_url).json()

        if 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp

        

    def getFeatureTypesWithIds(self):
        ''' 
        Returns:
            a dictionary of Feature types with feature name as keys and featureIds as values
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        feature_types = ''
        #jsn_resp = self.retrieveJson(self.base_url + "/featureType")
        jsn_resp, result_offset = self.retrieveJsonWithOffset(request_url="/featureType")
        

        if jsn_resp:
  
            feature_type_ids = jsn_resp['result']
            # Get the rest of the results, if all where not returned with the first query
            while result_offset:
                jsn_resp, result_offset = self.retrieveJsonWithOffset(request_url="/featureType", 
                                                                      offset=result_offset) 
                feature_type_ids.extend(jsn_resp['result'])
            
            feature_types = {}

            for feat_type_id in feature_type_ids:
                jsn_resp = self.retrieveJson(
                    self.base_url + "/featureType/" + str(feat_type_id))

                if jsn_resp:

                    feature_name = jsn_resp['result']['name']
                    feature_types[feature_name] = feat_type_id

        return feature_types
                

    def createFeatureType(self, featureTypeName, featureTypeDescription=''):
        ''' Creates a FeatureType
        Returns:
            'result' from the json response
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        post_data = {"featureType": {
            "description": featureTypeDescription, "name": featureTypeName}}

        resp = self.create(
            self.base_url + '/featureType', json.dumps(post_data))
        jsn_resp = resp.json()

        if 'result' and 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result']

    def deleteFeatureType(self, featureTypeId):
        ''' Deletes a FeatureType
        Returns:
            json response
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        request_url = "/featureType/" + str(featureTypeId)
        jsn_resp = self.delete(self.base_url + request_url).json()

        if 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp

    def getFeature(self, featureId):
        ''' Retrieves a Feature
        Args:
            featureId: id of a feature to be retrieved
        Returns:
            'result' from the json response OR nothing is feature was not found
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        jsn_resp = self.retrieveJson(
            self.base_url + "/feature/" + str(featureId))

        if jsn_resp:
            return jsn_resp['result']
        else:
            return ''

    def insertFeature(
            self, name, featureTypeId, featureLocations,
            sequenceId, description='', featureDefault=False, parentId=None):
        ''' Creates a Feature
        Args:
            name: name of the feature
        Returns:
            'result' from the json response
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        post_data = {
            "feature": {
                "name": name,
                "featureType": {"id": featureTypeId},
                "featureLocations": featureLocations,
                "description": description,
                "featureDefault": featureDefault
            }
        }
        if parentId is not None:
            post_data["feature"]["parentFeature"] = {"id": parentId}

        resp = self.create(self.base_url + "/sequence/" +
                           str(sequenceId) + "/feature", json.dumps(post_data))

        jsn_resp = resp.json()

        if 'result' and 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result']

    def deleteFeature(self, featureId):
        ''' Deletes a Feature
        Args:
            featureId: id of the feature to be deleted
        Returns:
            json response
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        request_url = "/feature/" + str(featureId)
        jsn_resp = self.delete(self.base_url + request_url).json()

        if 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp

           
    # TODO Not tested
    # TODO create/insert - review for consistency
    def createSpecimen(self):
        ''' Create a Specimen
        Args:
            TBD
        Kargs:
            TBD
        Returns:
            'result' from the json response
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        # TODO Not yet implemented
        logging.warn("createSpecimen not yet implemented")


    # TODO Not tested
    # TODO get/retrieve - review for consistency
    def getSpecimen(self, specimenId):
        ''' Retrieves a Specimen
        Args:
            specimenId: id of a specimen to be retrieved
        Kargs:
            None
        Returns:
            'result' from the json response OR nothing if specimen was not found
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        jsn_resp = self.retrieveJson(
            self.base_url + "/specimen/" + str(specimenId))

        if jsn_resp:
            return jsn_resp['result']
        else:
            return ''


    # TODO Not tested
    def updateSpecimen(self, specimenId, params):
        resp = self.update(
            self.base_url + "/specimen/" + str(specimenId),
            json.dumps(params))
        jsn_resp = resp.json()

        if 'result' and 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result'], jsn_resp['statusCode'], jsn_resp['message']


    def deleteSpecimen(self, specimenId):
        ''' Deletes a Feature
        Args:
            specimenId: id of the specimen to be deleted
        Returns:
            json response
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        request_url = "/specimen/" + str(specimenId)
        jsn_resp = self.delete(self.base_url + request_url).json()

        if 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp


    def getJsonSpecimenIds(self, params=None):
        ''' Gets ids of all specimens
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        jsn_resp = self.retrieveJson(
            self.base_url + "/specimen", params=params)

        
        if jsn_resp['count'] > 0 and not jsn_resp['result']:
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp


    def getJsonSpecimenIdsByOtherIds(self, code, identifier):
        params = {'filterName': 'otherIds',
                  'filterValue': code + identifier,
                  'filterOperator': 'and',
                  'filterWildcard': 'true'}

        return self.getJsonSpecimenIds(params)

    def getJsonSpecimenIdsBySpecimenId(self, code, identifier):
        params = {'filterName': ['biologicalCollection.name', 'number'],
                  'filterValue': [code, identifier],
                  'filterOperator': ['and', 'and'],
                  'filterWildcard': ['false', 'false']}

        jsn_resp = self.getJsonSpecimenIds(params)

        # Collection Name + Specimen Number should be unique
        if jsn_resp['count'] > 1:
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp

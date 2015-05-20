'''
Created on Mar 4, 2015

@author: korolo

Sequence DB Web Services module
'''

import requests, json, logging, base64, os
import urlparse
import gzip, mimetypes


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
            SeqDB API response (still need to format it to JSon, etc.) or nothing if resource was not found 
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
            resp.raise_for_status()
        
        
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
            
            #if 'count' and 'sortColum' and 'limit' and 'offset' and 'message' and 'statusCode' and 'sortOrder' not in jsn_resp:
            if 'result' not in jsn_resp.keys():
                raise UnexpectedContent(response=jsn_resp)
            
            # If these values are present, there is a chance that we can get paginated result
            if 'count' in jsn_resp.keys() and 'limit' in jsn_resp.keys():

                result_num = int(jsn_resp['count'])
                result_per_page =  int(jsn_resp['limit'])
                
                if (result_num > result_per_page):
                    
                    for curr_offset in xrange(result_per_page, result_num, result_per_page):
                        if params:
                            params['offset'] = curr_offset
                        else:
                            params={'offset': curr_offset}
                            
                        resp2 = self.retrieve(request_url, params)
                        jsn_resp2 =  json.loads(resp2.text)
                
                        if 'count' and 'limit' and 'result' not in jsn_resp2.keys():
                            raise UnexpectedContent(response=jsn_resp)
                        
                        jsn_resp['result'] = jsn_resp['result'] + jsn_resp2['result'] 
                        
                    jsn_resp['limit'] = jsn_resp['count']
        
            resp = jsn_resp
        
        return resp
    

    
    def update(self, request_url, json_data):
        ''' Updates a SeqDB entity
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError  
        '''
        req_header = { 'apikey': self.api_key, 'Content-Type': 'application/json' }
        resp = requests.put(request_url, headers=req_header, data=json_data)

        resp.raise_for_status()

        return resp


    
    def create(self, request_url, json_data):
        ''' Creates a SeqDB entity
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError  
        '''
        req_header = { 'apikey': self.api_key, 'Content-Type': 'application/json' }
        resp = requests.post(request_url, headers=req_header, data=json_data)  #(url, data, json)
        
        resp.raise_for_status()
        
        return resp



    def delete(self, request_url):
        ''' Creates a SeqDB entity
        Raises:
            requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout, and requests.exceptions.HTTPError  
        '''
        req_header = { 'apikey': self.api_key }
        resp = requests.delete(request_url, headers=req_header)
        
        # Will raise HTTPError exception if response status was not ok
        resp.raise_for_status()
        
        return resp


    
    def getJsonConsensusSequenceIds(self, params=None):
        ''' Gets sequence ids of all consensus sequences
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent  
        '''
        jsn_resp = self.retrieveJson(self.base_url + "/consensus", params=params)

        
        if jsn_resp['count'] > 0 and not jsn_resp['result']:
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp



    def getJsonConsensusSequenceIdsByName(self, name):
        params = {'filterName':'sequence.name',
                  'filterValue':name,
                  'filterOperator':'and',
                  'filterWildcard':'true'}

        return self.getJsonConsensusSequenceIds(params)



    def getJsonConsensusSequenceIdsByGI(self, gi):
            params = {'filterName':'sequence.genBankGI',
                      'filterValue':gi,
                      'filterOperator':'and',
                      'filterWildcard':'false'
                    }

            jsn_resp = self.getJsonConsensusSequenceIds(params)

            if jsn_resp['count'] > 1:
                raise SystemError("More than one record associated with gi, which should be unique")

            return jsn_resp


    # TODO verify which payload values are required by the SeqDB WS API
    def createConsensusSequence(self, name, sequence, qualities=None, seqType="N", readNum=0, additional=None):
        post_data = { 
            'consensus': {
                'name': name,
                'seq': sequence,
                'seqType': seqType,
                'readNum': readNum,
            },
        }

        if additional != None:
            post_data['consensus'].update(additional)

        resp = self.create(self.base_url + "/consensus", json.dumps(post_data))
        jsn_resp = resp.json()

        if 'result' and 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result'], jsn_resp['statusCode'], jsn_resp['message']


    
    # TODO Not tested
    def deleteConsensusSequence(self, consensus_id):
        request_url = "/consensus/" + str(consensus_id)
        jsn_resp = self.delete(self.base_url + request_url).json()
  
        if 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)
        
        return jsn_resp


    def importChromatSequences(self, blob, dest_file_name, notes="", trace_file_path=""):
        ''' Imports a binary blob (i.e. chromatogram) to seqdb
        Kwargs:
            dest_file_name name with which the chromatogram will be saved on sedDB side; i.e. expected to have .ab1 extension
        Returns:
            a list of seqdb sequence ids for the created sequences OR empty list if creation failed
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
        
        
        resp = self.create(self.base_url + '/sequenceImport', json.dumps(post_data))
        jsn_resp = resp.json()
        
        if 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)
        
        result = ""
        if 'result' in jsn_resp.keys():
            result = jsn_resp['result']
        else:
            logging.warn("Importing chromatogram failed with status code: '%s' and message: '%s'" % (jsn_resp['statusCode'], jsn_resp['message']) )
        
        return result



    def importChromatSequencesFromFile(self, chromat_file, notes="", trace_file_path="", dest_file_name=""):
        ''' Imports a chromatogram from a file
        Args:
            chromat_file name of the chromatogram file
        Returns:
            list of sequence ids of the sequences that were imported from the chromatogram
        Raises:
            IOError
            UnexpectedContent
        '''
        
        if not os.path.isfile(chromat_file):
            raise IOError("Expecting a file, but got a directory.")
        
        chromat_file_name = os.path.basename(chromat_file)
        
        if not dest_file_name:
            # get just the file name, no extension
            dest_file_name = os.path.splitext(os.path.basename(chromat_file_name))[0]
        
        if mimetypes.guess_type(chromat_file_name)[1] == "gzip":
            
            file_strem = gzip.open(chromat_file, "rb")
            
            # Remove .gz extension for file name
            chromat_file_name = chromat_file_name[:-3]

        else:    
            file_strem = open(chromat_file, "r")    
        
        blob = file_strem.read()
        
        file_strem.close()

        return self.importChromatSequences(blob = blob, dest_file_name = dest_file_name, notes=notes, trace_file_path=trace_file_path)



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



    def updateSeqSource(self, seqdb_id, params):
        resp = self.update(self.base_url + "/sequence/" + str(seqdb_id) + "/seqSource", json.dumps(params))
        jsn_resp = resp.json()

        if 'result' and 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result'], jsn_resp['statusCode'], jsn_resp['message']



    def getJsonSeq(self, seq_id):
        jsn_resp = self.retrieveJson(self.base_url + "/sequence/" + str(seq_id))
        if not jsn_resp['result']:
            raise UnexpectedContent(response=jsn_resp)
        
        jsn_seq = jsn_resp['result']
        if 'seq' not in jsn_seq.keys() or 'name' not in jsn_seq.keys():
            raise UnexpectedContent(response=jsn_resp)
        
        return jsn_seq

    def getFastaSeqPlus(self, seq_id):
        ''' Gets sequence from SeqDB and returns it in a fasta format with a unique header. Note, the fasta 
            formatting is done here, instead of using seqdb fasta web service request
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        jsn_seq = self.getJsonSeq(seq_id)
        # TODO Use BioPython to format
        fasta_seq =  '>' + jsn_seq['name'] + '|seqdbId:' + str(seq_id) + '\n' + jsn_seq['seq'] + '\n';
        return fasta_seq



    def getFastaSeq(self, seq_id):
        ''' Gets sequence in fasta format (SeqDB fasta request)
        Raises:
            requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout, and requests.exceptions.HTTPError
        '''
        url = self.base_url + "/sequence/" + str(seq_id) + ".fasta"
        response = self.retrieve(url)
        return response.content

    
    def getItsRegionIds(self):
        ''' Get region IDs of ITS sequences
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''  
        return self.getRegionIdsByName("ITS")
        #self.getRegionIds()



    def getRegionIdsByName(self, name):
        params = { 
                'filterName': 'name', 
                'filterValue': name, 
                'filterOperator': 'and',
                'filterWildcard': 'false' 
                }
        return self.getRegionIds(params)



    def getRegionIds(self, params=None):
        ''' Get region IDs
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''  
        jsn_resp = self.retrieveJson(self.base_url + "/region", params)
        if jsn_resp:
            return jsn_resp['result']
        else:
            return ''



    def createRegion(self, name, description):
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
        # TODO Question requirement for gene region to be associated with a group
        post_data = {"region":{"description":description, "group":{"id":1}, "name":name }}
 
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



    def getSeqIds(self, region_id):
        ''' Given a region id, return sequence ids, belonging to this region
        Returns:
            a list of seqdb sequence ids for the created sequences OR empty list id created failed
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        request_url = "/region/" + str(region_id) + "/sequence"
        jsn_resp = self.retrieveJson(self.base_url + request_url)
        
        if jsn_resp:
            return jsn_resp['result']
        else:
            return ''



    def getFeatureTypesWithIds(self):
        ''' 
        Returns:
            a dictionary of Feature types: name: featureId
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        
        jsn_resp = self.retrieveJson(self.base_url + "/featureType")
        
        if jsn_resp:
            
            feature_type_ids = jsn_resp['result']
            feature_types = {}
            
            for feat_type_id in feature_type_ids:
                jsn_resp = self.retrieveJson(self.base_url + "/featureType/" + str(feat_type_id))
                
                if jsn_resp:
                
                    feature_name = jsn_resp['result']['name']
                    feature_types[feature_name] = feat_type_id 
                
            return feature_types
    
    
        else:
            return ''

    def createFeatureType(self, featureTypeName, featureTypeDescription = ''):
        ''' Creates a FeatureType 
        Returns:
            'result' from the json response
        Raises:
            requests.exceptions.ConnectionError
            requests.exceptions.ReadTimeout
            requests.exceptions.HTTPError
            UnexpectedContent
        '''
        post_data = {"featureType":{"description":featureTypeDescription,"name":featureTypeName }}
        
        resp = self.create(self.base_url + '/featureType', json.dumps(post_data))
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
        jsn_resp = self.retrieveJson(self.base_url + "/feature/" + str(featureId))
        
        if jsn_resp:
            return jsn_resp['result']
        else:
            return ''



    def insertFeature(self, name, featureTypeId, featureLocations, sequenceId, description='', featureDefault=False):
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
                "name":name,
                "featureType":{"id":featureTypeId},
                "featureLocations":featureLocations,
                "description":description,
                "featureDefault":featureDefault
            }
        }
        
        resp = self.create(self.base_url + "/sequence/" + str(sequenceId) + "/feature", json.dumps(post_data))
                     
        #print("name: " + name + "  featureTypeId: " + featureTypeId + "   featureLocations: " + featureLocations+"  description: " + description+ "   featureDefault: " + str(featureDefault) + "\n\n")
        
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
    
    

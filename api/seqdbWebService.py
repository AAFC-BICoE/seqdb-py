'''
Created on Mar 4, 2015

@author: korolo

Sequence DB Web Services module
'''

import requests, json, logging


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

    
        
        
    # Submits a request to SeqDB web services
    # Raises requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout, and requests.exceptions.HTTPError  
    # @request_url": full request url
    # Returns response (still need to format it to JSon, etc.) or nothing if resource was not found 
    def retrieve(self, request_url, params=None):
        
        req_header = { 'apikey': self.api_key }
        resp = requests.get(request_url, headers=req_header, params=params)
        
        if resp.status_code == 404:
            resp = ''
        else:
            # Will raise HTTPError exception if response status was not ok
            resp.raise_for_status()
              
        return resp
    
    
    # Submits a request to SeqDB web services
    # Raises requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout, and requests.exceptions.HTTPError  
    # Returns json formatted object
    def retrieveJson(self, request_url, params=None):
        resp = self.retrieve(request_url, params=params)
        if resp:
            return json.loads(resp.text)
        else:
            return resp
    
    def update(self, request_url, json_data):
        req_header = { 'apikey': self.api_key, 'Content-Type': 'application/json' }
        resp = requests.put(request_url, headers=req_header, data=json_data)

        resp.raise_for_status()

        return resp
    
    def create(self, request_url, json_data):
        req_header = { 'apikey': self.api_key, 'Content-Type': 'application/json' }
        resp = requests.post(request_url, headers=req_header, data=json_data)  #(url, data, json)
        
        resp.raise_for_status()
        
        return resp
    
    def delete(self, request_url):
        req_header = { 'apikey': self.api_key }
        resp = requests.delete(request_url, headers=req_header)
        
        # Will raise HTTPError exception if response status was not ok
        resp.raise_for_status()
        
        return resp
    
    def getJSONConsensusSequenceIds(self, params=None):
        jsn_resp = self.retrieveJson(self.base_url + "/consensus", params=params)

        if 'count' and 'sortColum' and 'limit' and 'offset' and 'message' and 'statusCode' and 'sortOrder' not in jsn_resp:
            raise UnexpectedContent(response=jsn_resp)

        if jsn_resp['count'] > 0 and ('result' not in jsn_resp.keys() or not jsn_resp['result']):
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp

    def getJSONConsensusSequenceIdsByName(self, name):
        params = {'filterName':'sequence.name',
                  'filterValue':name,
                  'filterOperator':'and',
                  'filterWildcard':'true'}

        return self.getJSONConsensusSequenceIds(params)

    def getJSONConsensusSequenceIdsByGI(self, gi):
            params = {'filterName':'sequence.genBankGI',
                      'filterValue':gi,
                      'filterOperator':'and',
                      'filterWildcard':'false'
                    }

            jsn_resp = self.getJSONConsensusSequenceIds(params)

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

    def updateSeqSource(self, seqdb_id, params):
        resp = self.update(self.base_url + "/sequence/" + str(seqdb_id) + "/seqSource", json.dumps(params))
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
            
    def getJSONSeq(self, seq_id):
        jsn_resp = self.retrieveJson(self.base_url + "/sequence/" + str(seq_id))
        if 'result' not in jsn_resp.keys() or not jsn_resp['result']:
            raise UnexpectedContent(response=jsn_resp)
        
        jsn_seq = jsn_resp['result']
        if 'seq' not in jsn_seq.keys() or 'name' not in jsn_seq.keys():
            raise UnexpectedContent(response=jsn_resp)
        
        return jsn_seq

    # Gets sequence from SeqDB and returns it in a fasta format with a unique header.
    # Note, the fast formatting is done here, instead of using seqdb fasta web service request
    # Raises requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout, and requests.exceptions.HTTPError
    def getFastaSeqPlus(self, seq_id):
        jsn_seq = getJSONSeq(seq_id)
        # TODO Use BioPython to format
        fasta_seq =  '>' + jsn_seq['name'] + '|seqdbId:' + str(seq_id) + '\n' + jsn_seq['seq'] + '\n';
        return fasta_seq
       
    # Gets sequence in fasta format (SeqDB fasta request)
    # Raises requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout, and requests.exceptions.HTTPError  
    def getFastaSeq(self, seq_id):
        url = self.base_url + "/sequence/" + str(seq_id) + ".fasta"
        response = self.retrieve(url)
        return response.content
        
    
    # Get region IDs of ITS sequences
    # Raises requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout, and requests.exceptions.HTTPError  
    def getItsRegionIds(self):
        return self.getRegionIdsByName("ITS")

    def getRegionIdsByName(self, name):
        params = { 
                'filterName': 'name', 
                'filterValue': name, 
                'filterOperator': 'and',
                'filterWildcard': 'false' 
                }
        return self.getRegionIds(params)

    def getRegionIds(self, params):
        jsn_resp = self.retrieveJson(self.base_url + "/region", params)
        if jsn_resp:
            if 'result' not in jsn_resp.keys():
                raise UnexpectedContent(response=jsn_resp)
           
            return jsn_resp['result']
        else:
            return ''

    def createRegion(self, name, description):
        # TODO Question requirement for gene region to be associated with a group
        post_data = {"region":{"description":description, "group":{"id":1}, "name":name }}
 
        resp = self.create(self.base_url + '/region', json.dumps(post_data))
        jsn_resp = resp.json()

        if 'result' and 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)

        return jsn_resp['result']



        
    # Given a region id, return sequence ids, belonging to this region
    # api_key and base_url required for ws request
    # Raises requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout, and requests.exceptions.HTTPError  
    def getSeqIds(self, region_id):
        request_url = "/region/" + str(region_id) + "/sequence"
        jsn_resp = self.retrieveJson(self.base_url + request_url)
        
        if jsn_resp:
            if not jsn_resp or 'result' not in jsn_resp.keys():
                raise UnexpectedContent(response=jsn_resp)
            
            return jsn_resp['result']
        else:
            return ''
    
    # Returns a dictionary of Feature types: featureName: featureId    
    # Raises requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout, and requests.exceptions.HTTPError  
    def getFeatureTypesWithIds(self):
        jsn_resp = self.retrieveJson(self.base_url + "/featureType")
        
        if jsn_resp:
            if 'result' not in jsn_resp.keys():
                raise UnexpectedContent(response=jsn_resp)
            
            feature_type_ids = jsn_resp['result']
            feature_types = {}
            
            for feat_type_id in feature_type_ids:
                jsn_resp = self.retrieveJson(self.base_url + "/featureType/" + str(feat_type_id))
                
                if jsn_resp:
                    if 'result' not in jsn_resp.keys():
                        raise UnexpectedContent(response=jsn_resp)
                
                    feature_name = jsn_resp['result']['featureName']
                    feature_types[feature_name] = feat_type_id 
                
            return feature_types
    
    
        else:
            return ''
    
    def createFeatureType(self, featureTypeName, featureTypeDescription = ''):
        post_data = {"featureType":{"featureDescription":featureTypeDescription,"featureName":featureTypeName }}
        
        resp = self.create(self.base_url + '/featureType', json.dumps(post_data))
        jsn_resp = resp.json()
        
        if 'result' and 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)
        
        return jsn_resp['result']


        
    def deleteFeatureType(self, featureTypeId):
        request_url = "/featureType/" + str(featureTypeId)
        jsn_resp = self.delete(self.base_url + request_url).json()
  
        if 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)
        
        return jsn_resp
        

  
    def getFeature(self, featureId):
        jsn_resp = self.retrieveJson(self.base_url + "/feature/" + str(featureId))
        
        if jsn_resp:
            if 'result' not in jsn_resp.keys():
                raise UnexpectedContent(response=jsn_resp)
            
            return jsn_resp['result']
    
        else:
            return ''


    
    def insertFeature(self, name, featureTypeId, featureLocations, sequenceId, description='', featureDefault=False):
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
        request_url = "/feature/" + str(featureId)
        jsn_resp = self.delete(self.base_url + request_url).json()
  
        if 'statusCode' and 'message' not in jsn_resp.keys():
            raise UnexpectedContent(response=jsn_resp)
        
        return jsn_resp
    
    

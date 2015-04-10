'''
Created on Feb 12, 2015

@author: korolo
'''
import unittest, requests, yaml
from api import seqdbWebService

#from seqdb_ws import json_seqdb_request, seqdb_ws_request
    
    # curl -H "apikey: ***REMOVED***" ***REMOVED***/sequence/1
    # curl -H "apikey: ***REMOVED***" ***REMOVED***/region?filterName=name&filterValue=ITS&filterOperator=and&filterWildcard=true
    
class TestSeqdbWebService(unittest.TestCase):

    def setUp(self):
        config = yaml.load(file('config.yaml', 'r'))
        self.seqdbWS = seqdbWebService.seqdbWebService(config['seqdb']['apikey'], config['seqdb']['url'])
        
        
    def tearDown(self):
        pass


    def testConnection(self):
        print "Cleanup done"



class TestSeqdbWebService_NoDataSetup(unittest.TestCase):
    
    def setUp(self):
        config = yaml.load(file('config.yaml', 'r'))
        self.seqdbWS = seqdbWebService.seqdbWebService(config['seqdb']['apikey'], config['seqdb']['url'])
        self.featureTypeIds = []

    
    def tearDown(self):
        for ftId in self.featureTypeIds:
            self.seqdbWS.deleteFeatureType(ftId)


    def testConnection(self):
        self.assertRaises(requests.exceptions.ConnectionError, self.seqdbWS.retrieve, "http://jibberish")
        # TODO: test wrong api key


    def testCreateFeature_valid(self):
        """Test creation of a feature - expected to pass"""
        featureTypeId = self.seqdbWS.createFeatureType("Test", "Test")
        self.featureTypeIds.append(featureTypeId)
        self.assertTrue(featureTypeId, "Feature type ID was not returned after feature type creation.")


    def testCreateMultipleFeatureTypes_valid(self):
        """Test creation of two feature - expected to work because of unique names"""
        featureTypeId = self.seqdbWS.createFeatureType("Test1", "Test")
        self.featureTypeIds.append(featureTypeId)
        featureTypeId = self.seqdbWS.createFeatureType("Test2", "Test")
        self.featureTypeIds.append(featureTypeId)
        self.assertTrue(featureTypeId, "Second feature type ID was not returned after feature type creation.")

    
    def testCreateMultipleFeature_failduplicate(self):
        """Test creation of a feature - expected to fail because of duplicate feature type id"""
        featureTypeId = self.seqdbWS.createFeatureType("Test", "Tesdst")
        self.featureTypeIds.append(featureTypeId)
        
        self.assertRaises(requests.exceptions.HTTPError, self.seqdbWS.createFeatureType, "Test", "Test")
        #featureTypeId = self.seqdbWS.createFeatureType("Test", "Test")
    

class TestSeqdbWebService_FeatureType_Existing(unittest.TestCase):

    
    def setUp(self):
        config = yaml.load(file('config.yaml', 'r'))
        self.seqdbWS = seqdbWebService.seqdbWebService(config['seqdb']['apikey'], config['seqdb']['url'])
        featureTypeId = self.seqdbWS.createFeatureType("TestType", "Test")
        
        self.featureTypeIds = [featureTypeId]
    
    def tearDown(self):
        for ftId in self.featureTypeIds:
            self.seqdbWS.deleteFeatureType(ftId)
    
    def testGetFeatureTypes(self):
        """Test retrieving of a feature - expected to pass"""
        # curl -H "apikey: ***REMOVED***"  "***REMOVED***:8181/seqdb.web-2.5/api/v1/featureType"
        # curl -H "apikey: ***REMOVED***"  "***REMOVED***/api/v1/featureType"
        actual = self.seqdbWS.getFeatureTypesWithIds()
        self.assertTrue(actual, "No feature types returned.")
        self.assertIn("TestType", actual, "No feature type 'Test' found.")

    def testDeleteFeatureType(self):
        """Test deletion of a feature - expected to pass"""
        idToDelete = self.featureTypeIds[0]
        actual = self.seqdbWS.deleteFeatureType(idToDelete)
        self.assertEqual(200, actual['statusCode'], "Could not delete feature type.")
        self.featureTypeIds.remove(idToDelete)



class TestSeqdbWebService_Region_Existing(unittest.TestCase):

    
    def setUp(self):
        config = yaml.load(file('config.yaml', 'r'))
        self.seqdbWS = seqdbWebService.seqdbWebService(config['seqdb']['apikey'], config['seqdb']['url'])
        regionId = self.seqdbWS.createRegion("ITS", "Test")
        self.regionIds = [regionId]
    
    def tearDown(self):
        for rgId in self.regionIds:
            self.seqdbWS.deleteRegion(rgId)


    def testGetItsRegionIds(self):
        # curl -H "apikey: ***REMOVED***"  "***REMOVED***:8181/seqdb.web-2.5/api/v1/region?filterName=name&filterValue=ITS&filterOperator=and&filterWildcard=true"
        # curl -H "apikey: ***REMOVED***"  "***REMOVED***/api/v1/region?filterName=name&filterValue=ITS&filterOperator=and&filterWildcard=true"
        actual = self.seqdbWS.getItsRegionIds()
        self.assertTrue(actual, "No ITS region ids returned.")
        self.assertIn(self.regionIds[0], actual, "Region id is not in the results.")




class TestSeqdbWebService_Sequence_Existing(unittest.TestCase):

    
    def setUp(self):
        config = yaml.load(file('config.yaml', 'r'))
        self.seqdbWS = seqdbWebService.seqdbWebService(config['seqdb']['apikey'], config['seqdb']['url'])
        seqId, errCod, msg = self.seqdbWS.createConsensusSequence("Test", "ACGTCTGATCGATC")
        self.sequenceIds = [seqId]
    
    def tearDown(self):
        for seqId in self.sequenceIds:
            self.seqdbWS.deleteConsensusSequence(seqId)
            
    def testGetFastaSeq(self):
        actual = self.seqdbWS.getFastaSeq(self.sequenceIds[0])
        self.assertTrue(actual, "Fasta sequence is empty.")
        self.assertIn(">", actual, "Fasta does not contain >.")        
    
    def testGetFastaSeqPlus(self):
        # curl -H "apikey: ***REMOVED***"  "***REMOVED***:8181/seqdb.web-2.5/api/v1/sequence/1"
        actual = self.seqdbWS.getFastaSeqPlus(self.sequenceIds[0])
        self.assertTrue(actual, "Fasta sequence is empty.")
        self.assertIn(">", actual, "Fasta does not contain >.")
        self.assertIn("seqdbId:", actual, "Fasta does not contain seqdbId.")
        

    '''
    # TODO: this test requires sequence to be associated with region. The ability to do this via API is not yet implemented,
    #       therefore we can not do a test setup for this test. 
    def testGetSeqIds(self):
        actual = self.seqdbWS.getSeqIds(regionId)
        self.assertTrue(actual, "No Sequence ids returned.")
        self.assertIn(self.sequenceIds[0], actual, "Expected sequence id is not in the results.")
    '''

#### Works till here

class TestSeqdbWebService_Sequence_FeatureType_Feature_Existing(unittest.TestCase):

    def setUp(self):
        config = yaml.load(file('config.yaml', 'r'))
        self.seqdbWS = seqdbWebService.seqdbWebService(config['seqdb']['apikey'], config['seqdb']['url'])
        
        seqId, code, msg = self.seqdbWS.createConsensusSequence("Test", "ACsdgafGTCTGATCGATC")
        self.sequenceIds = [seqId]
        
        featureTypeId = self.seqdbWS.createFeatureType("Test", "Test")
        self.featureTypeIds = [featureTypeId]

        self.sampleFeatureLocations1 = [{"start":3,"end":5,"frame":1,"strand":1}]
        featureId = self.seqdbWS.insertFeature("testName", featureTypeId, self.sampleFeatureLocations1, seqId, "sample description", True)
        self.featureIds = [featureId]
    
    
    def tearDown(self):

        for fId in self.featureIds:
            self.seqdbWS.deleteFeature(fId)

        for seqId in self.sequenceIds:
            self.seqdbWS.deleteConsensusSequence(seqId)
        
        for ftId in self.featureTypeIds:
            self.seqdbWS.deleteFeatureType(ftId)
    
    

    def testGetFeature(self):
        featureId = self.featureIds[0]
        retrieved_feature = self.seqdbWS.getFeature(featureId)
        self.assertTrue(retrieved_feature, "No feature was retrieved.") 
        self.assertEqual("testName", retrieved_feature['name'], "Name of the retrieved feature does not match.")
        self.assertEqual("sample description", retrieved_feature['description'], "Feature description does not match")
        self.assertEqual(self.sampleFeatureLocations1, retrieved_feature['featureLocations'], "")
        self.assertEqual(self.featureTypeIds[0], retrieved_feature['featureType']['id'], "Feature type id does not match.")
        self.assertTrue(retrieved_feature['featureDefault'], "Feature default does not match")
    
    
    def testCreateFeature(self):

        sampleFeatureLocations2 = [{"start":3,"end":5,"frame":1,"strand":1}, {"start":334,"end":454,"frame":2,"strand":1}]
        
        featureId = self.seqdbWS.insertFeature("testName", self.featureTypeIds[0], sampleFeatureLocations2, self.sequenceIds[0], "sample description", True)
        self.assertTrue(featureId, "Feature id was not returned after feature creation.")
        self.featureIds.append(featureId)
            
        
    def testDeleteFeature(self):
        featureId = self.featureIds.pop()
        actual = self.seqdbWS.deleteFeature(featureId)
        self.assertEqual(200, actual['statusCode'], "Could not delete feature type.")
        
        # Get feature again
        retrieved_feature = self.seqdbWS.getFeature(featureId)
        self.assertFalse(retrieved_feature, "Unexpected: Feature was found after being deleted.") 
        
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

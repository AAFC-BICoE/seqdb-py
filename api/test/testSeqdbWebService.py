'''
Created on Feb 12, 2015

@author: korolo
'''
import unittest, requests
from api import seqdbWebService

#from seqdb_ws import json_seqdb_request, seqdb_ws_request

class Test(unittest.TestCase):
    
    # curl -H "apikey: ***REMOVED***" http://localhost:10000/seqdb.web/api/v1/sequence/1
    # curl -H "apikey: ***REMOVED***" http://localhost:10000/seqdb.web/api/v1/sequence/1.fasta
    # curl -H "apikey: ***REMOVED***" http://localhost:10000/seqdb.web/api/v1//region?filterName=name&filterValue=ITS&filterOperator=and&filterWildcard=true
    # curl -H "apikey: ***REMOVED***" http://localhost:8181/seqdb.web-2.5/api/v1/region?filterName=name&filterValue=ITS&filterOperator=and&filterWildcard=true


    def setUp(self):

        #local: "***REMOVED***:8181/seqdb.web-2.5/api/v1"
        #prod: "***REMOVED***/api/v1"
        
        self.seqdbWS = seqdbWebService.seqdbWebService('***REMOVED***', '***REMOVED***:8181/seqdb.web-2.5/api/v1')
    
    def tearDown(self):
        pass

    def testRetrieve(self):
        self.assertRaises(requests.exceptions.ConnectionError, self.seqdbWS.retrieve, "http://jibberish")
        # TODO: test wrong api key



    def testGetItsRegionIds(self):
        # curl -H "apikey: ***REMOVED***"  "***REMOVED***:8181/seqdb.web-2.5/api/v1/region?filterName=name&filterValue=ITS&filterOperator=and&filterWildcard=true"
        # curl -H "apikey: ***REMOVED***"  "***REMOVED***/api/v1/region?filterName=name&filterValue=ITS&filterOperator=and&filterWildcard=true"
        actual = self.seqdbWS.getItsRegionIds()
        self.assertTrue(actual, "No ITS region ids returned.")
        self.assertIn(8, actual, "Region id 8 is not in the results.")
        
                
        
    def testGetSeqIds(self):
        actual = self.seqdbWS.getSeqIds("8")
        self.assertTrue(actual, "No Sequence ids returned.")
        self.assertIn(1685, actual, "Sequence id 1685 is not in the results.")
        
        
        
    def testGetFastaSeq(self):
        actual = self.seqdbWS.getFastaSeq("1")
        self.assertTrue(actual, "Fasta sequence is empty.")
        self.assertIn(">", actual, "Fasta does not contain >.")        


    
    def testGetFastaSeqPlus(self):
        # curl -H "apikey: ***REMOVED***"  "***REMOVED***:8181/seqdb.web-2.5/api/v1/sequence/1"
        actual = self.seqdbWS.getFastaSeqPlus("1")
        self.assertTrue(actual, "Fasta sequence is empty.")
        self.assertIn(">", actual, "Fasta does not contain >.")
        self.assertIn("seqdbId:1", actual, "Fasta does not contain seqdbId.")
        

    
    def testGetFeatureTypes(self):
        # curl -H "apikey: ***REMOVED***"  "***REMOVED***:8181/seqdb.web-2.5/api/v1/featureType"
        # curl -H "apikey: ***REMOVED***"  "***REMOVED***/api/v1/featureType"
        actual = self.seqdbWS.getFeatureTypesWithIds()
        self.assertTrue(actual, "No feature types returned.")
        self.assertIn("Quality Trim", actual, "No feature type 'Quality Trim' found.")
        
    
    def testCreateDeleteFeatureType(self):
        #curl -X POST -H "apikey: ***REMOVED***" -H "Content-Type: application/json" -d '{"featureType":{"featureDescription":"test description 1231q","featureName":"test type 123123"}}' "***REMOVED***:8181/seqdb.web-2.5/api/v1/featureType"
        
        #self.seqdbWS.deleteFeatureType("19")
        
        featureTypeID = self.seqdbWS.createFeatureType("test type", "test type description")
        self.assertTrue(featureTypeID, "Feature type ID was not returned after feature type creation.")
        
        #curl -X DELETE -H "apikey: ***REMOVED***" "***REMOVED***:8181/seqdb.web-2.5/api/v1/featureType/6"
        actual = self.seqdbWS.deleteFeatureType(featureTypeID)
        self.assertEqual(200, actual['statusCode'], "Could not delete feature type.")

    def testCreateGetDeleteFeature(self):
        
        # Create feature
        sampleFeatureLocations1 = [{"start":3,"end":5,"frame":1,"strand":1}]
        sampleFeatureLocations2 = [{"start":3,"end":5,"frame":1,"strand":1}, {"start":334,"end":454,"frame":2,"strand":1}]
        
        featureId = self.seqdbWS.insertFeature("testName", 1, sampleFeatureLocations1, 1, "sample description", True)
        self.assertTrue(featureId, "Feature id was not returned after feature creation.")
        
        # Get feature
        retrieved_feature = self.seqdbWS.getFeature(featureId)
        self.assertTrue(retrieved_feature, "No feature was retrieved.") 
        self.assertEqual("testName", retrieved_feature['name'], "Name of the retrieved feature does not match.")
        self.assertEqual("sample description", retrieved_feature['description'], "Feature description does not match")
        self.assertEqual(sampleFeatureLocations1, retrieved_feature['featureLocations'], "")
        self.assertEqual(1, retrieved_feature['featureType']['id'], "Feature type id does not match.")
        self.assertTrue(retrieved_feature['featureDefault'], "Feature default does not match")
        
        # Delete feature
        actual = self.seqdbWS.deleteFeature(featureId)
        self.assertEqual(200, actual['statusCode'], "Could not delete feature type.")
        
        # Get feature again
        retrieved_feature = self.seqdbWS.getFeature(featureId)
        self.assertFalse(retrieved_feature, "Unexpected: Feature was found after being deleted.") 
                
                
                
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
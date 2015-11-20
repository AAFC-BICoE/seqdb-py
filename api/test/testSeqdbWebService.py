'''
Created on Feb 12, 2015

@author: korolo
'''
import unittest, requests, yaml
from api import seqdbWebService


test_url = '***REMOVED***:2002/seqdb/api/v1/'
test_api_key = '***REMOVED***'
    
class TestSeqdbWebService(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        self.fixture = seqdbWebService.seqdbWebService(api_key=test_api_key, base_url=test_url)

    def setUp(self):
        pass
        
    def tearDown(self):
        pass


    def testConnection(self):
        print "Cleanup done"



class TestSeqdbWebService_NoDataSetup(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        self.fixture = seqdbWebService.seqdbWebService(api_key=test_api_key, base_url=test_url)

    def setUp(self):
        self.featureTypeIds = []
        self.sequenceIds = []

    
    def tearDown(self):
        for ftId in self.featureTypeIds:
            self.fixture.deleteFeatureType(ftId)

        for sId in self.sequenceIds:
            self.fixture.deleteSequence(sId)

    def testConnection(self):
        self.assertRaises(requests.exceptions.ConnectionError, self.fixture.retrieve, "http://jibberish")
        # TODO: test wrong api key


    def testCreateFeature_valid(self):
        """Test creation of a feature - expected to pass"""
        featureTypeId = self.fixture.createFeatureType("Test", "Test")
        self.featureTypeIds.append(featureTypeId)
        self.assertTrue(featureTypeId, "Feature type ID was not returned after feature type creation.")


    def testCreateMultipleFeatureTypes_valid(self):
        """Test creation of two feature - expected to work because of unique names"""
        featureTypeId = self.fixture.createFeatureType("Test1", "Test")
        self.featureTypeIds.append(featureTypeId)
        featureTypeId = self.fixture.createFeatureType("Test2", "Test")
        self.featureTypeIds.append(featureTypeId)
        self.assertTrue(featureTypeId, "Second feature type ID was not returned after feature type creation.")

    
    def testCreateMultipleFeature_failduplicate(self):
        """Test creation of a feature - expected to fail because of duplicate feature type id"""
        featureTypeId = self.fixture.createFeatureType("Test", "Tesdst")
        self.featureTypeIds.append(featureTypeId)
        
        self.assertRaises(requests.exceptions.HTTPError, self.fixture.createFeatureType, "Test", "Test")
   
        
    def testCreateChromatSequence_valid(self):
        """ Test creating a sequence with binary .abi or .ab1 file (chromatogram) """
        seq_id = self.fixture.importChromatSequencesFromFile(chromat_file = "data/GRDI_test_seq.ab1", notes="This is a test upload.", 
                    trace_file_path="test_path",dest_file_name="test.ab1",)
        
        self.sequenceIds.extend(seq_id)
        
        self.assertTrue(seq_id, "Persisting chromatogram did not return an id.")

    def testCreateChromatSequenceGz_valid(self):
        """ Test creating a sequence with a zipped binary file (chromatogram) i.e. ab1.gz """
        
        seq_id = self.fixture.importChromatSequencesFromFile(chromat_file = "data/blob_db.ab1.gz")
        self.sequenceIds.extend(seq_id)
        
        self.assertTrue(seq_id, "Persisting chromatogram did not return an id.")
        
        
    def testCreateChromatSequence_wrongPath(self):
        """ Test creating a sequence with non-existent file (should throw exception)"""
        self.assertRaises(IOError, self.fixture.importChromatSequencesFromFile, "data/")
        self.assertRaises(IOError, self.fixture.importChromatSequencesFromFile, "zzz/non-existent.ab1")
        


    def testDelete(self):
        #curl -X DELETE -H "apikey: ***REMOVED***" "***REMOVED***/sequence/4709479"
        #curl -X GET -H "apikey: ***REMOVED***" "***REMOVED***/sequence/4709479"
        '''
        delete_ids = (
4709635,
4709636,
            )
        self.fixture.bulkDeleteSequence(delete_ids)
        '''
        
        #self.fixture.deleteSequence(4709476)
        #self.fixture.deleteSequence(4709481)
        pass
        


class TestSeqdbWebService_FeatureType_Existing(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        self.fixture = seqdbWebService.seqdbWebService(api_key=test_api_key, base_url=test_url)
    
    def setUp(self):
        featureTypeId = self.fixture.createFeatureType("TestType", "Test")
        
        self.featureTypeIds = [featureTypeId]

    
    def tearDown(self):
        for ftId in self.featureTypeIds:
            self.fixture.deleteFeatureType(ftId)

    
    def testGetFeatureTypes_valid(self):
        """Test retrieval of a feature - expected to pass"""
        actual = self.fixture.getFeatureTypesWithIds()
        self.assertTrue(actual, "No feature types returned.")
        self.assertIn("TestType", actual, "No feature type 'Test' found.")


    def testDeleteFeatureType_valid(self):
        """Test deletion of a feature - expected to pass"""
        idToDelete = self.featureTypeIds[0]
        actual = self.fixture.deleteFeatureType(idToDelete)
        self.assertEqual(200, actual['statusCode'], "Could not delete feature type.")
        self.featureTypeIds.remove(idToDelete)



class TestSeqdbWebService_Region_Existing(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        self.fixture = seqdbWebService.seqdbWebService(api_key=test_api_key, base_url=test_url)
    
    def setUp(self):
        regionId = self.fixture.createRegion("ITS", "Test", 168)
        self.regionIds = [regionId]

    
    def tearDown(self):
        for rgId in self.regionIds:
            self.fixture.deleteRegion(rgId)


    
    def testGetRegionIds_valid(self):
        """Test retrieval of an ITS region ids - expected to pass"""
        
        current_ids, offset = self.fixture.getRegionIdsWithOffset()
        self.assertTrue(current_ids, "No region ids returned.")

        id_found = False
        while offset:
            current_ids, offset = self.fixture.getRegionIdsWithOffset(offset=offset)
            if self.regionIds[0] in current_ids:
                id_found = True
                break
            
        self.assertTrue(id_found, "Region id is not in the results.")
 
    
    def testGetItsRegionIds_valid(self):
        """Test retrieval of an ITS region ids - expected to pass"""
        actual = self.fixture.getItsRegionIds()
        self.assertTrue(actual, "No ITS region ids returned.")
        self.assertIn(self.regionIds[0], actual, "Region id is not in the results.")



class TestSeqdbWebService_Sequence_Existing(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        self.fixture = seqdbWebService.seqdbWebService(api_key=test_api_key, base_url=test_url)
    
    def setUp(self):
        seqId, errCod, msg = self.fixture.createConsensusSequence(name="Test", sequence="ACGTCTGATCGATC")
        self.sequenceIds = [seqId]

    
    def tearDown(self):
        for seqId in self.sequenceIds:
            self.fixture.deleteConsensusSequence(seqId)

        
    def testGetFastaSeq(self):
        """Test retrieval of a sequence in fasta format (seqdb api fasta) - expected to pass"""
        actual = self.fixture.getFastaSeq(self.sequenceIds[0])
        self.assertTrue(actual, "Fasta sequence is empty.")
        self.assertIn(">", actual, "Fasta does not contain >.")        

    
    def testGetFastaSeqPlus(self):
        """Test retrieval of a sequence in fasta format (python fasta formatting) - expected to pass"""
        # curl -H "apikey: ***REMOVED***"  "***REMOVED***:8181/seqdb.web-2.5/api/v1/sequence/1"
        actual = self.fixture.getFastaSeqPlus(self.sequenceIds[0])
        self.assertTrue(actual, "Fasta sequence is empty.")
        self.assertIn(">", actual, "Fasta does not contain >.")
        self.assertIn("seqdbId:", actual, "Fasta does not contain seqdbId.")

    
    # TODO: takes a long time to run: re-design to be more targeted 
    def testGetSeqIds(self):
        current_ids, offset = self.fixture.getSequenceIdsWithOffset(sequenceName="Test")
        self.assertTrue(current_ids, "No Sequence ids returned.")

        id_found = False
        while offset:
            current_ids, offset = self.fixture.getSequenceIdsWithOffset(sequenceName="Test", offset=offset)
            if self.sequenceIds[0] in current_ids:
                id_found = True
                break
            
        self.assertTrue(id_found, "Expected sequence id is not in the results.")
    



class TestSeqdbWebService_Sequence_FeatureType_Feature_Existing(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.fixture = seqdbWebService.seqdbWebService(api_key=test_api_key, base_url=test_url)

    def setUp(self):
        seqId, code, msg = self.fixture.createConsensusSequence("Test", "ACsdgafGTCTGATCGATC")
        self.sequenceIds = [seqId]
        
        featureTypeId = self.fixture.createFeatureType("Test", "Test")
        self.featureTypeIds = [featureTypeId]

        self.sampleFeatureLocations1 = [{"start":3,"end":5,"frame":1,"strand":1}]
        featureId = self.fixture.insertFeature("testName", featureTypeId, self.sampleFeatureLocations1, seqId, "sample description", True)
        self.featureIds = [featureId]
    
    
    def tearDown(self):

        for fId in self.featureIds:
            self.fixture.deleteFeature(fId)

        for seqId in self.sequenceIds:
            self.fixture.deleteConsensusSequence(seqId)
        
        for ftId in self.featureTypeIds:
            self.fixture.deleteFeatureType(ftId)
    

    def testGetFeature(self):
        """Test retrieval of a feature - expected to pass"""
        featureId = self.featureIds[0]
        retrieved_feature = self.fixture.getFeature(featureId)
        self.assertTrue(retrieved_feature, "No feature was retrieved.") 
        self.assertEqual("testName", retrieved_feature['name'], "Name of the retrieved feature does not match.")
        self.assertEqual("sample description", retrieved_feature['description'], "Feature description does not match")
        self.assertEqual(self.sampleFeatureLocations1, retrieved_feature['featureLocations'], "")
        self.assertEqual(self.featureTypeIds[0], retrieved_feature['featureType']['id'], "Feature type id does not match.")
        self.assertTrue(retrieved_feature['featureDefault'], "Feature default does not match")
    
    
    def testCreateFeature(self):
        """Test creation of a feature - expected to pass"""
        sampleFeatureLocations2 = [{"start":3,"end":5,"frame":1,"strand":1}, {"start":334,"end":454,"frame":2,"strand":1}]
        
        featureId = self.fixture.insertFeature("testName", self.featureTypeIds[0], sampleFeatureLocations2, self.sequenceIds[0], "sample description", True)
        self.assertTrue(featureId, "Feature id was not returned after feature creation.")
        self.featureIds.append(featureId)
            
        
    def testDeleteFeature(self):
        """Test deletion of a feature - expected to pass"""
        featureId = self.featureIds.pop()
        actual = self.fixture.deleteFeature(featureId)
        self.assertEqual(200, actual['statusCode'], "Could not delete feature type.")
        
        # Get feature again
        retrieved_feature = self.fixture.getFeature(featureId)
        self.assertFalse(retrieved_feature, "Unexpected: Feature was found after being deleted.") 
        
    
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

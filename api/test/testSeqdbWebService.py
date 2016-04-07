'''
Created on Feb 12, 2015

Note: these tests should be run against a local SeqDB instance only. They do not mock the API. 
    Since the work on SeqDB API is ongoing, they serve as both tests for the request methods 
    as well as tests for the SeqDB API itself.
    
    Therefore, these tests are brittle. They depend not only on the current content of your 
    local SeqDB, but also on the information that you have access to (i.e. groups you're 
    part of). Here, by "you" I mean the user whose API Key is used for testing (seqdb_api_key)

Before running the tests:
    - Make suse that config/config4tests.yaml is created from sample and the modified to point
     to local instance of SeqDB
    - Make sure that local instance of SeqDB is running

@author: korolo
'''
import requests
import unittest 
import yaml

from api import seqdbWebService
from config import config_root


class TestSeqdbWebService(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        config = yaml.load(file(config_root.path() + '/config4tests.yaml', 'r'))        
        self.fixture = seqdbWebService.seqdbWebService(api_key=config['seqdb']['api_key'],
                                                   base_url=config['seqdb']['api_url'])
    
    def testRetrieve(self):
        # Test faulty connection
        self.assertRaises(requests.exceptions.ConnectionError, self.fixture.retrieve, "http://jibberish")
        # TODO: test wrong api key
        
    def testRetrieveJson(self):
        actual = self.fixture.retrieveJson("/consensus/358301")
        self.assertTrue(actual, "retrieveJson: no result was returned.")
        
    def testRetrieveJsonWithOffset(self):
        actual_firstpage, _ = self.fixture.retrieveJsonWithOffset("/consensus?filterName=sequence.name&filterValue=ium&filterWildcard=true")
        self.assertTrue(actual_firstpage, "No result was returned.")
        self.assertIn('metadata', actual_firstpage, "No 'metadata' in the response.")

        actual_secondpage, _ = self.fixture.retrieveJsonWithOffset(request_url="/consensus?filterName=sequence.name&filterValue=ium&filterWildcard=true",
                                                                offset=20)
        self.assertTrue(actual_secondpage, "No result was returned.")
        self.assertIn('metadata', actual_secondpage, "No 'metadata' in the response.")
        
        self.assertNotEqual(actual_firstpage['result'], actual_secondpage['result'], "")
        
    
    ###########################################################################
    # Sequence
    ###########################################################################
    
    # Refactored
    def testGetFastaSequencesWithOffset(self):
        actual = self.fixture.getRawSequencesWithOffset(offset=0, limit=30, sequence_format="fasta", specimenNum=4405)
        self.assertTrue(actual, "No Sequences returned.")
        self.assertIn(">seqdb|27755", actual,"Expecting that fasta return will contain id 27755.")
        self.assertNotIn(">seqdb|358301", actual, "Fasta return is not expected to have sequence 358301, since it is consensus.")
    
    # Refactored    
    def testGetSequenceIds(self):
        actual = self.fixture.getRawSequenceIds(specimenNum=4405)
        self.assertTrue(actual, "No Sequence ids returned.")
        self.assertEqual(22, len(actual),"Expecting 22 sequences associated with this specimen.")
        self.assertIn(27755, actual, "Sequence id 27755 is expected to be associated with specimen 4405.")
        self.assertNotIn(358301, actual, "Sequence id 358301 is not expected to be in results, since it is consensus.")  

    # Refactored
    def testCreateChromatSequence_wrong_path(self):
        self.assertRaises(IOError, self.fixture.importChromatSequencesFromFile, "data/")
        self.assertRaises(IOError, self.fixture.importChromatSequencesFromFile, "zzz/non-existent.ab1")

    # Refactored
    def testCreateDeleteChromatSequence(self):
        """ Test creating a sequence with binary .abi or .ab1 file (chromatogram) """
        seq_id = self.fixture.importChromatSequencesFromFile(chromat_file = "data/GRDI_test_seq.ab1", notes="This is a test upload.", 
                    trace_file_path="test_path",dest_file_name="test.ab1",)
        
        self.assertTrue(seq_id, "Persisting chromatogram did not return an id.")

        # Delete
        delete_jsn_resp = self.fixture.deleteRawSequence(seq_id)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete feature type.")
    
    # Refactored
    def testCreateDeleteChromatSequenceGz_valid(self):
        """ Test creating a sequence with a zipped binary file (chromatogram) i.e. ab1.gz """
        seq_id = self.fixture.importChromatSequencesFromFile(chromat_file = "data/blob_db.ab1.gz")
        self.assertTrue(seq_id, "Persisting chromatogram did not return an id.")
        
        # Delete
        delete_jsn_resp = self.fixture.deleteRawSequence(seq_id)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete feature type.")     
        
    # Refactored
    def testGetFastaSeq(self):
        actual = self.fixture.getFormattedSeq("1", "fasta")
        self.assertTrue(actual, "Fasta sequence is empty.")
        self.assertIn(">", actual, "Fasta does not contain >.")        
    
    # Refactored
    def testGetFastqSeq(self):
        actual = self.fixture.getFormattedSeq("1", "fastq")
        self.assertTrue(actual, "Fastq sequence is empty.")
        self.assertIn("@seqdb", actual, "Fastq does not contain @seqdb.")        
    
    # Refactoring: Not carrying this method forward
    def testGetFastaSeqPlus(self):
        # "http://localhost:2002/seqdb\/api/v1/sequence/1"
        actual = self.fixture.getFastaSeqPlus("1")
        self.assertTrue(actual, "Fasta sequence is empty.")
        self.assertIn(">", actual, "Fasta does not contain >.")
        self.assertIn("seqdbId:1", actual, "Fasta does not contain seqdbId.")        


    ###########################################################################
    # Consensus Sequence
    ###########################################################################
    
    #Refactored
    def testGetConsensusSequenceIds(self):
        actual = self.fixture.getConsensusSequenceIds(specimenNum=4405)
        self.assertTrue(actual, "No Sequence ids returned.")
        self.assertEqual(1, len(actual),"Expecting 1 consensus sequence associated with this specimen.")
        self.assertIn(358301, actual, "Sequence id 358301 is expected to be associated with specimen 4405.")
    #Refactored
    def testCreateGetDeleteSequence(self):
        # create
        seqId, errCod, msg = self.fixture.createConsensusSequence(name="Test", sequence="ACGTCTGATCGATC")
        self.assertTrue(seqId, "Creating consensus sequence did not return an id.")
        self.assertEqual(errCod, 201, "Did not get successful exit code for create consensus sequence.")
        
        # get
        seqIds = self.fixture.getConsensusSequenceIds(sequenceName="Test")
        self.assertTrue(seqIds, "Creating consensus sequence did not return an id.")
        self.assertIn(seqId, seqIds, "Expected sequence id was not in the result.")
        
        # delete
        delete_jsn_resp = self.fixture.deleteConsensusSequence(seqId)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete feature type.")
    
    
    ###########################################################################
    # Determination
    ###########################################################################
    
    #Refactored
    def testCreateGetDeleteDetermination(self):
        test_taxonomy = {'superkingdom': 'Eukaryota', 'genus': 'Phytophthora', 'species': 'Phytophthora ramorum', 'order': 'Peronosporales'}
        det_id = self.fixture.insertSequenceDetermination(28954, test_taxonomy)
        self.assertTrue(det_id, "Creating determination on a sequence did not return a determination id.")
        
        get_determination = self.fixture.getDetermination(det_id)
        self.assertTrue(get_determination, "Did not get determination back.")
        self.assertEquals(False, get_determination['accepted'], "")
        self.assertEquals('Phytophthora', get_determination['taxonomy']['genus'], "")
        
        # delete
        delete_jsn_resp = self.fixture.deleteDetermination(det_id)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete determination.")
    
    #Refactored
    def testGetAcceptedSpecimenDetermination(self):
        actual = self.fixture.getAcceptedSpecimenDetermination(27755)
        self.assertTrue(actual, "Expecting accepted determination, but got none.")
        self.assertEquals("Fungi", actual['taxonomy']['kingdom'], "Expecting kingdom Fungi, but got {}".format(actual['taxonomy']['kingdom']))
        self.assertEquals("arrhenomanes", actual['taxonomy']['species'], "Expecting species arrhenomanes, but got {}".format(actual['taxonomy']['species']))

        actual = self.fixture.getAcceptedSpecimenDetermination(358301)
        self.assertTrue(actual, "Expecting accepted determination, but got none.")
        self.assertEquals("Oomycota", actual['taxonomy']['phylum'], "Expecting phylum Oomycota, but got {}".format(actual['taxonomy']['phylum']))
        self.assertEquals("Pythiales", actual['taxonomy']['taxanomicOrder'], "Expecting order Pythiales, but got {}".format(actual['taxonomy']['taxanomicOrder']))
        
        
    
    ###########################################################################
    # Gene Region
    ###########################################################################
    
    #Refactored
    def testCreateGetDeleteITSRegion(self):
        # create
        create_regId = self.fixture.createRegion("ITS_test", "Test", 168)
        self.assertTrue(create_regId, "Creating region did not return an id.")
        
        # get
        get_regIds = self.fixture.getItsRegionIds()
        self.assertTrue(get_regIds, "No ITS region ids returned.")
        self.assertIn(create_regId, get_regIds, "Expected region is not in the results.")

        # delete
        delete_jsn_resp = self.fixture.deleteRegion(create_regId)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete region.")     
        
        # get
        get_regIds = self.fixture.getItsRegionIds()
        self.assertNotIn(create_regId, get_regIds, "Not expecting this region to be in the results. Region was deleted.")
        
    
    ###########################################################################
    # Feature
    ###########################################################################
    def testCreateGetDeleteFeature(self):
        
        # Create feature
        sampleFeatureLocations1 = [{"start":3,"end":5,"frame":1,"strand":1}]
        sampleFeatureLocations2 = [{"start":3,"end":5,"frame":1,"strand":1}, {"start":334,"end":454,"frame":2,"strand":1}]
        
        featureId = self.fixture.insertFeature("testName", 1, sampleFeatureLocations1, 1, "sample description", True)
        self.assertTrue(featureId, "Feature id was not returned after feature creation.")
        
        # Get feature
        retrieved_feature = self.fixture.getFeature(featureId)
        self.assertTrue(retrieved_feature, "No feature was retrieved.") 
        self.assertEqual("testName", retrieved_feature['name'], "Name of the retrieved feature does not match.")
        self.assertEqual("sample description", retrieved_feature['description'], "Feature description does not match")
        self.assertEqual(sampleFeatureLocations1, retrieved_feature['featureLocations'], "")
        self.assertEqual(1, retrieved_feature['featureType']['id'], "Feature type id does not match.")
        self.assertTrue(retrieved_feature['featureDefault'], "Feature default does not match")
        
        # Delete feature
        delete_jsn_resp = self.fixture.deleteFeature(featureId)
        self.assertEqual(200, delete_jsn_resp['metadata']['statusCode'], "Could not delete feature type.")
        
        # Get feature again
        retrieved_feature = self.fixture.getFeature(featureId)
        self.assertFalse(retrieved_feature, "Unexpected: Feature was found after being deleted.") 
                
    
    
    ###########################################################################
    # Feature Type
    ###########################################################################
    def testGetFeatureTypes(self):
        # "http://localhost:2002/seqdb\/api/v1/featureType"
        actual = self.fixture.getFeatureTypesWithIds()
        self.assertTrue(actual, "No feature types returned.")
        self.assertIn("Quality Trim", actual, "No feature type 'Quality Trim' found.")
        
    
    def testCreateDeleteFeatureType(self):
        #curl -X POST -H "apikey: ***REMOVED***" -H "Content-Type: application/json" -d '{"featureType":{"featureDescription":"test description 1231q","featureName":"test type 123123"}}' "***REMOVED***/featureType"
        
        #self.fixture.deleteFeatureType("19")
        
        featureTypeID = self.fixture.createFeatureType("test type", "test type description")
        self.assertTrue(featureTypeID, "Feature type ID was not returned after feature type creation.")
        
        #curl -X DELETE -H "apikey: ***REMOVED***" "***REMOVED***/featureType/6"
        actual = self.fixture.deleteFeatureType(featureTypeID)
        self.assertEqual(200, actual['metadata']['statusCode'], "Could not delete feature type.")

    def testCreateDeleteFeatureType_multiple_duplicate(self):
        #self.fixture.deleteFeatureType("8")
        # Create
        featureTypeId1 = self.fixture.createFeatureType("Test1", "Test")
        featureTypeId2 = self.fixture.createFeatureType("Test2", "Test")
        self.assertTrue(featureTypeId2, "Second feature type ID was not returned after feature type creation.")
        """Test creation of a feature - expected to fail because of duplicate feature type id"""
        self.assertRaises(requests.exceptions.HTTPError, self.fixture.createFeatureType, "Test1", "Duplicate of Test1")
   
        # Delete
        actual = self.fixture.deleteFeatureType(featureTypeId1)
        self.assertEqual(200, actual['metadata']['statusCode'], "Could not delete feature type.")
        actual = self.fixture.deleteFeatureType(featureTypeId2)
        self.assertEqual(200, actual['metadata']['statusCode'], "Could not delete feature type.")
      

    
    ###########################################################################
    # Project Tag
    ###########################################################################
    
    
    ###########################################################################
    # Specimen
    ###########################################################################

    def testGetSpecimen(self):
        specimen_jsn = self.fixture.getSpecimen("1322")
        self.assertTrue(specimen_jsn, "Specimen was not returned.")
        self.assertEqual(1322, specimen_jsn["id"])
            

                
if __name__ == "__main__":
    unittest.main()
'''
Created on April 29, 2016

@author: Oksana Korol
'''
import unittest
from tools.test.test_pull_seqdb_seqs import TestPullSeqdbSeqs
from tools.test.test_seqdb_config_maker import TestSeqdbConfigMaker
from tools.test.test_push_delete_seqdb_its_feat import TestPushDeleteSeqdbFeatures

#suite = unittest.TestLoader().loadTestsFromTestCase(TestPushDeleteSeqdbFeatures)
#unittest.TextTestRunner(verbosity=2).run(suite)

def test_script_suite():
    test_suite = unittest.TestSuite()
    #test_suite.addTest(unittest.makeSuite(TestPullSeqdbSeqs))
    test_suite.addTest(unittest.makeSuite(TestSeqdbConfigMaker))
    test_suite.addTest(unittest.makeSuite(TestPushDeleteSeqdbFeatures))
    return test_suite

script_suite = test_script_suite()
runner = unittest.TextTestRunner()
runner.run(script_suite)

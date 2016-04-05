'''
Created on Mar 30, 2016

@author: Oksana Korol
'''
import unittest
from api.test.TestSpecimenApi import TestSpecimenApi
from api.test.TestProjectTagApi import TestProjectTagApi
from api.test.TestRawSequenceApi import TestRawSequenceApi

def test_suite():
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.makeSuite(TestSpecimenApi))
    test_suite.addTest(unittest.makeSuite(TestProjectTagApi))
    test_suite.addTest(unittest.makeSuite(TestRawSequenceApi))
    return test_suite

api_suite = test_suite()
runner = unittest.TextTestRunner()
runner.run(api_suite)
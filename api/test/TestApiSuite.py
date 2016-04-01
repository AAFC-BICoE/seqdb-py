'''
Created on Mar 30, 2016

@author: Oksana Korol
'''
import unittest
from api.test import TestSpecimenApi, TestProjectTagApi, TestSequenceApi

def test_suite():
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.makeSuite(TestSpecimenApi.TestSpecimenApi))
    test_suite.addTest(unittest.makeSuite(TestProjectTagApi.TestProjectTagApi))
    test_suite.addTest(unittest.makeSuite(TestSequenceApi.TestSequenceApi))
    return test_suite

api_suite = test_suite()
runner = unittest.TextTestRunner()
runner.run(api_suite)
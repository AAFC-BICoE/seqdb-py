'''
Created on April 29, 2016

@author: Oksana Korol
'''
import unittest
from api.test.TestApiSuite import test_api_suite
from tools.test.TestScriptSuite import test_script_suite

# Below needs fixing. To run all tests, run the api and tools suites separately:
# api/test/TestApiSuite.py
# tools/test/TestScriptSuite.py
"""
api_suite = test_api_suite()
tools_suite = test_script_suite()
runner = unittest.TextTestRunner()
print "***** API Test Suite *****"
runner.run(api_suite)
print "***** Tools Test Suite *****"
runner.run(tools_suite)
"""
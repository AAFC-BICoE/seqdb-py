import os
import sys

# Allow access to api/ package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from seqdb_py.api.test import TestApiSuite as api_suite
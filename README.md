dev branch
==========
The dev branch of this project currently isn't functional.  It is being 
modified to follow the conventions of and to utilize the tools available to 
Python 3.

At the moment, work is focused on making the API functional with Python 3.  The 
tools directory hasn't been modified and likely is not compatible with the rest
of the project.

This message will be amended as work progresses.

_______


SeqDB-Py 
=========

SeqDB-Py is a Python extension of SeqDB project. It allows access to SeqDB 
through its RESful API. The project contains both the API wrapper and the 
command line tools that use it.

How do I get set up?
--------------------

Requirements:
  * Running instance of SeqDB (not currently Open Source)
  * Python 3.6
  * pip 10+ (should come packaged with Python 3.4+)

Running
-------

This section is for those who just want to run (i.e. use) seqdb-py. If you'd 
like to develop it, see instruction below. 

* Create and start your virtual environment:

   > virtualenv venv
   
   > source venv/bin/activate

* Install seqdb-py from git (use one of the options below):

-- Source distribution:
  
   > pip install -e git+https://bitbucket.org/aafc-mbb/seqdb-py.git#egg=seqdb-py

-- Installing a particular version (i.e. Version 1.0 in the example below):

   > pip install -e git+https://bitbucket.org/aafc-mbb/seqdb-py.git@1.0#egg=seqdb-py

* Modify config files for your environment:

Initial distribution will contain <file_name>.yaml.sample files, which need to 
be renamed to <file_name>.yaml and modified according to instructions inside 
the files.

There is a main config file for the package, config/config.yaml, which 
specifies logging information, SeqDB API url, etc. In addition, some tools may 
have specific yaml config files, containing parameters that only those tools 
require. These configuration files will be located alongside the tool.  

* Run desired tool (below are available as entry points):
   > pull_seqdb_its_seqs
   
   > push_seqdb_its_feat
   
   > delete_seqdb_features
   
   > seqdb_gb_insert 
   
   Note: see help on each tool for usage


Example: Running seqdb_gb_insert.py 
  This example assumes you've activated a virtual environment and installed s
  eqdb-py in it.
  
* Create config/config.yaml and tools/seqdb_gb_insert_config.yaml

   > cp config/config.yaml.sample config/config.yaml
   
   > cp tools/seqdb_gb_insert_config.yaml.sample tools/seqdb_gb_insert_config.yaml

* Initialize above files according to instructions inside each .yaml file

* Run script

   > seqdb_gb_insert




Developing
----------

If you would like to contribute to seqdb-py codebase, do the following.

* Clone the repo
    > git clone https://github/AAFC-BICoE/seqdb-py.git

* Switch to the seqdb-py directory
    > cd seqdb-py

* Create virtual environment, source, it and install required packages in it:
   > mkdir env
   
   > python3 -m venv env
   
   > source env/bin/activate
   
   > pip install -r requirements.txt
   
* Run the tests before changing and before committing code.


Tests
-----

Seqdb-py has a testing framework, implemented using python's PyUnit. To keep 
the code maintainable and to minimise frustraction at the merging time, 
please run all tests before you touch the code and before you commit. 
Here is what you need:

* Running local instance of SeqDB. Rationale: since testing seqdb-py 
essentially means testing that SeqDB API calls will work with current SeqDB 
API, you need a current SeqDB applicaion. Also, since the tests create and 
delete features, it is preferable that the SeqDB instance is local. Even 
though the tests are designed not to leave any mess behind, there is no 
guarantee. 

* Admin access to local instance of SeqDB.

* config/config4tests.yaml file. Create it using a supplied template.

* Run the tests for api and tools suite separately (the combined suite is not 
working at the moment):

   >api/test/TestApiSuite.py
   
   >tools/test/TestScriptSuite.py

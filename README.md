# SeqDB-Py 

SeqDB-Py is a Python extension of SeqDB project. It allows access to SeqDB 
through its RESful API. The project contains ~~both~~ the API wrapper ~~and the 
command line tools that use it~~.

The project was previously set up to work with Python2 (tag python2 and previous). The
transition to Python3 has started, but is still a work in progress.  For this reason,
many of its modules do not work.  At the present time, **only the RawSeqenceApi module and
its dependencies** have been tested and are deemed functional.  


## How do I get set up?

Requirements:
  * Running instance of SeqDB (not currently Open Source)
  * Python 3.6+
  * pip 10+

### Running

This section is for those who just want to run (i.e. use) seqdb-py. If you'd 
like to develop, see the instruction below. 

* Create and start your virtual environment:

        mkdir .env
        python3 -m venv .env
        source .env/bin/activate

* Install seqdb-py from git (use one of the options below):

    Source distribution:
  
        pip install git+https://github.com/AAFC-BICoE/seqdb-py.git#egg=seqdb-py

    Installing a particular version (i.e. version 2.0 in the example below):

        pip install -e git+https://github.com/AAFC-BICoE/seqdb-py.git@v2.0#egg=seqdb-py

* Modify config files for your environment:

The initial distribution contains `<file_name>.yaml.sample` files, which need to 
be renamed `<file_name>.yaml` and modified according to instructions inside 
the files.

There is a main config file for the package, config/config.yaml, which 
specifies logging information, SeqDB API url, etc. ~~In addition, some tools may 
have specific yaml config files, containing parameters that only those tools 
require. These configuration files will be located alongside the tool.~~  

### Developing

If you would like to contribute to seqdb-py codebase, do the following:

* Clone the repo
        
        $ git clone https://github/AAFC-BICoE/seqdb-py.git

* Switch to the seqdb-py directory

        $ cd seqdb-py

* Create virtual environment, source, it and install required packages in it:

        seqdb-py$ mkdir .env
        seqdb-py$ python3 -m venv .env
        seqdb-py$ source .env/bin/activate
        seqdb-py$ pip install -r requirements.txt
   
* Run the tests before changing and before committing code.


## Tests

Seqdb-py has a testing framework implemented using python's PyUnit. It is advised 
to run all tests before making changes to the code and before committing.

To run the tests, you will need:

* Running local instance of the SeqDB application: Since testing seqdb-py 
essentially means testing that SeqDB API calls will work with the current SeqDB 
API, a working instance of the SeqDB application is required. Also, since the tests 
create and delete features, it is imperative that these tests not be run on the 
production instance of SeqDB. Even though the tests are designed not to leave any mess 
behind, there is no guarantee. 

* Admin access to local instance of SeqDB (recommended).

* config/config4tests.yaml file. Create it using a supplied template.

* Run the tests for api suite ~~and tools suite separately (the combined suite is not 
working at the moment)~~:

        seqdb-py$ python api/test/TestApiSuite.py   
     >~~$ tools/test/TestScriptSuite.py~~

# README #

Tool to extract chromatograms from GeneSifter LIMS system (GSLE) and load it to SeqDB.

### Background ###
The result of the GSLE workflow is a number of chromatograms (.abi or .ab1 file format). 
These chromatograms need to be loaded to SeqDB and linked back to the sample/specimen 
that were used to create this run. This project is going to be triggered from the GSLE
using their event/handler functionality. Once triggered, the projec will load completed
chromatograms to SeqDB.

## How do I get set up? ##

### Requirements ###
* Running instance of SeqDB (not currently Open Source)
* Running instance of GSLE (not Open Source)
* Python 2.7
* pip (sudo apt-get install python-pip python-dev build-essential)
* virtualenv (sudo apt-get install python-virtualenv)
* Following will be installed using Makefile:
  * python requests module (2.6.0)
  * python psycopg2 module (2.4)
  * python PyYAML module (3.11)
  * See requirements.txt for full / current list of requirements

### Running ###

* Create config.yaml (use provided config.yaml.sample as a template)
* Build project:
  * make
* Run project:
  * run.sh


## Project file description ##

* extract_LIMS_chromats.py: main project file. Contains instructions to extract data from GSLE and load it to SeqDB.
* api/seqdbWebService.py: python-wrapped SeqDB API calls. Used in main project file to communicate with SeqDB.
* config.yaml: configuration file required by main project file. Contains particulars about GSLE, SeqDB, logging, etc.
* config.yaml.sample: sample of configuration above. Use it to create config.yaml for the new project.
* requirements.txt: file containing 3rd party python modules that need to be installed before the project execution.
* README.txt: this file. Information about this project.
* Makefile: describes metadata about this project and is used to build it
* run.sh: run script that invokes main project file with all the necessary environments. (Registered as an event handler in GSLE.) This is a system specific file and is here only to give an example on how this project can be executed. 

## Who do I talk to? ##
* Oksana Korol
* MBB team

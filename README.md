# README #

This repo contains a Python module to interact with the SeqDB WS API, and python scripts that use this module to interact with SeqDB.


## How do I get set up? ##

### Requirements ###
* Running instance of SeqDB (not currently Open Source)
* Python 2.7
* python requests module (2.6.0)
* python biopython module (1.65)
* python PyYAML module (3.11)
* See requirements.txt for full / current list of requirements

### Not required, but recommended ###
* pip (sudo apt-get install python-pip python-dev build-essential)
* virtualenv (sudo apt-get install python-virtualenv)
* make

### Running ###
* Install required modules (preferably in virtual env):
   > pip install -r requirements.txt
* Add project dir to PYTHONPATH
   > PYTHONPATH=“/path/to/project_dir:$PYTHONPATH"
   > export PYTHONPATH
* Run desired tool

### Example: Running seqdb_gb_insert.py from a clean directory using make and a virtual environment
* Create config.yaml file
```
cp config.yaml.sample config.yaml
```
* Initialize config.yaml file
  * Update entrez email address entry
  * Update seqdb url
  * Update seqdb apikey
* Create and activate virtual environment
```
make
source venv/bin/activate
```
* Update PYTHONPATH
```
   PYTHONPATH=“`pwd`:$PYTHONPATH"
   export PYTHONPATH
```
* Run script
```
./bin/seqdb_gb_insert.py
```

### Subsequent use of seqdb_gb_insert.py
* Update PYTHONPATH
```
   PYTHONPATH=“`pwd`:$PYTHONPATH"
   export PYTHONPATH
```
* Activate virtual environment
```
source venv/bin/activate
```
* Run script
```
./bin/seqdb_gb_insert.py
```

## Who do I talk to? ##
* Oksana Korol
* MBB team

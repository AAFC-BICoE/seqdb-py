# README #

This repo contains Python scripts and modules that extend Sequence Database functionality.


## How do I get set up? ##

### Requirements ###
* Python 2.7
* python requests module

### Not required, but recommended ###
* pip (sudo apt-get install python-pip python-dev build-essential)
* virtualenv (sudo apt-get install python-virtualenv)
* make

### Running ###
* Install python requests (preferably in virtual env):
   > pip install -r requirements.txt
* Add project dir to PYTHONPATH
   > PYTHONPATH=“/path/to/project_dir:$PYTHONPATH"
   > export PYTHONPATH

### Running seqdb_gb_insert.py
* cp config.yaml.sample config.yaml
 * Update config.yaml

### Linux Requirements (Needs to be updated / integrated)
* sudo apt-get install python-pip python-dev build-essential
* sudo apt-get install python-virtualenv

### Linux Running (Needs to be updated / integrated)
* sudo make
* source sdbgbload/bin/activate
* cp config.yaml.sample config.yaml
 * Update config.yaml
* ./bin/seqdb_gb_insert.py

## Who do I talk to? ##
* Oksana Korol
* MBB team

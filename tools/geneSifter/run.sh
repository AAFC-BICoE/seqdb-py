#!/bin/bash

#############################################################################
# This is a runner script for a python program, which loads chromatograms   #
# from GeneSifter LIMS (GSLE) and pushed them to SeqDB.                     #
#																			#
# Author: Oksana Korol														#
# Date created: 20.04.2015													#
#############################################################################

## This directory is where python code, its virtual environment and config is. 
## See config.yaml for where logs and chromatograms will be written.
cd /home/geospiza/scripts/uploadChromatsToSeqdb

source /usr/local/geospiza/bin/genv
export LD_LIBRARY_PATH=/usr/local/geospiza/lib/

# Activate virtual environment for the python code to run. The virtual 
# environment should be created with ">make" for this project (see README)
source venv/bin/activate


python extract_LIMS_chromats.py $@


deactivate

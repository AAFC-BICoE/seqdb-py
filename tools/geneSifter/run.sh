#!/bin/bash

#source /usr/local/geospiza/bin/genv
source venv/bin/activate

python /home/geospiza/scripts/uploadChromatsToSqdb/extract_LIMS_chromats.py $@

deactivate
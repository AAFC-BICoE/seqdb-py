#!/bin/bash

date >> /tmp/test.log
echo "New run activated" >> /tmp/test.log
source /usr/local/geospiza/bin/genv
echo "source /usr/local/geospiza/bin/genv  - Done" >> /tmp/test.log
source venv/bin/activate
echo "source venv/bin/activate  - Done" >> /tmp/test.log
export LD_LIBRARY_PATH=/usr/local/geospiza/lib/
echo "export LD_LIBRARY_PATH=/usr/local/geospiza/lib/  - Done. Attempting to run a script" >> /tmp/test.log

python /home/geospiza/scripts/uploadChromatsToSqdb/extract_LIMS_chromats.py $@

echo "Script finished execution." >> /tmp/test.log

deactivate

echo "Deactivated venv. Execution complete." >> /tmp/test.log

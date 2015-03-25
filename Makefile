setup: venv

venv:
	# somehow check python version
	# somehow permit python interpretter to be provided as a parameter to Makefile
	#virtualenv -p ${PYTHON} venv
	virtualenv venv
	venv/bin/pip install -r requirements.txt
	echo "Run 'source venv/bin/activate' before executing"

clean:
	rm -rf venv

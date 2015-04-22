.PHONY: setup
setup: venv

venv:
	# TODO: check python version and fail if it is not python 2
	virtualenv venv
	venv/bin/pip install -r requirements.txt
	echo "Run 'source venv/bin/activate' before executing"

.PHONY: clean-log
clean-log:
	rm *.log *.err

.PHONY: clean
clean: clean-log
	rm -rf venv

pep8.log:
	pep8 --show-source --show-pep8 tools/ > pep8.log || true
	pep8 --show-source --show-pep8 api/ > pep8.log || true

pylint.log:
	pylint tools/ > pylint.log || true
	pylint api/ >> pylint.log || true
    
.PHONY: check-source
check-source: pep8.log pylint.log

.PHONY: check-commit
check-commit:
	@pylint tools/ > pylint.log 2>pylint.err && \
	pylint api/ >> pylint.log 2>>pylint.err  && \
	pep8 --show-source --show-pep8 tools/ > pep8.log 2>pep8.err && \
	pep8 --show-source --show-pep8 api/ > pep8.log 2>>pep8.err && \
	echo "pylint and pep8 passed - okay to commit" || \
	echo "pylint and pep8 failed to pass - not okay to commit"

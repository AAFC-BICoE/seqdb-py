setup: sdbgbload

sdbgbload:
	virtualenv venv
	venv/bin/pip install -r requirements.txt

clean:
	rm -rf sdbgbload

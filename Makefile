setup: sdbgbload

sdbgbload:
	virtualenv sdbgbload
	sdbgbload/bin/pip install -r requirements.txt

clean:
	rm -rf sdbgbload

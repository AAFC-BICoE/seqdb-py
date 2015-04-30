#!/usr/bin/env python
#!/usr/bin/python
'''
Created on Apr 27, 2015

@author: korolo
'''

#from distutils.core import setup
from setuptools import setup

setup(name='seqdb-py-api',
    version='1.0',
    description='SeqDB Python API.',
    author='Oksana Korol',
    author_email='oksana.korol@agr.gc.ca',
    url='https://bitbucket.org/aafc-mbb/seqdb-py',
    #data_files=[('', 'requirements.txt')],
    py_modules = ['seqdbWebService'],
    install_requires = [
        'requests==2.6.0', 
     ]
)

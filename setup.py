#!/usr/bin/env python
#!/usr/bin/python
'''
Created on Apr 27, 2015

@author: korolo
'''

#from distutils.core import setup
from setuptools import setup, find_packages

setup(name='seqdb-py',
    version='1.5',
    description='SeqDB Python Extensions: API and tools',
    author='Oksana Korol',
    author_email='oksana.korol@agr.gc.ca',
    url='https://bitbucket.org/aafc-mbb/seqdb-py',
    
    packages = find_packages(exclude=['test','*.test', '*.test.*']),
    #packages=['', 'api', 'tools'],
    #data_files= [('tools', ['config.yaml.sample',]),('', ['config.yaml.sample',]),],
    include_package_data = True,
    package_data = {
                    '': ['*.yaml.sample', 'README'],
                    },
    install_requires = [
        'numpy==1.11.0',
        'biopython==1.67',
        'PyYAML==3.11',
        'requests==2.10', 
    ],
    entry_points={
    'console_scripts': [
        'seqdb_config_maker = tools.seqdb_config_maker:main',
        'pull_seqdb_seqs = tools.pull_seqdb_seqs:main',
        'push_to_seqdb = tools.push_to_seqdb:main',
        'delete_seqdb_features = tools.delete_seqdb_features:main',
        'seqdb_gb_insert = tools.seqdb_gb_insert:main',
    ],
},
)

#!/usr/bin/env python
#!/usr/bin/python
'''
Created on Apr 27, 2015

@author: korolo
'''

#from distutils.core import setup
from setuptools import setup

setup(name='seqdb-py',
    version='1.1-SNAPSHOT',
    description='SeqDB Python Extensions: API and tools',
    author='Oksana Korol',
    author_email='oksana.korol@agr.gc.ca',
    url='https://bitbucket.org/aafc-mbb/seqdb-py',
    packages=['api', 'tools'],
    data_files= [('', ['config.yaml.sample',]),],
    install_requires = [
        'biopython==1.65',
        'PyYAML==3.11',
        'requests==2.6.0', 
     ],
    entry_points={
    'console_scripts': [
        'seqdb_config_maker = tools.seqdb_config_maker:main',
        'pull_seqdb_its_seqs = tools.pull_seqdb_its_seqs:main',
        'push_seqdb_its_feat = tools.push_seqdb_its_feat:main',
        'delete_seqdb_features = tools.delete_seqdb_features:main',
        'seqdb_gb_insert = tools.seqdb_gb_insert:main',
    ],
},
)

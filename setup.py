#!/usr/bin/env python3
#!/usr/bin/python3

from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()
setup(
    name='seqdb_py',
    version='2.0',
    description='SeqDB Python3 API Extensions',
    long_description=long_description,
    long_description_content_type='text/markdown',
    # url='https://github.com/AAFC-BICoE/seqdb-py',
    license='MIT',
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    packages=find_packages(exclude=['seqdb_py.tools', 'seqdb_py.tools.test']),
    # packages=['seqdb_py'],
    include_package_data=True,
    package_data={
        'config': ['*.yaml.sample', ],
        '': ['README.md', ],
    },
    install_requires=[
        # 'numpy==1.11.0',
        # 'biopython==1.67',
        'PyYAML>=3.13',
        'requests>=2.19',
    ],
)
    #entry_points={
    #    'console_scripts': [
    #        'seqdb_config_maker = tools.seqdb_config_maker:main',
    #        'pull_seqdb_seqs = tools.pull_seqdb_seqs:main',
    #        'push_to_seqdb = tools.push_to_seqdb:main',
    #        'delete_seqdb_features = tools.delete_seqdb_features:main',
    #        'seqdb_gb_insert = tools.seqdb_gb_insert:main',
    #   ],
    # },

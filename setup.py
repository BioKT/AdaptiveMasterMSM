#!/usr/bin/env python

# Setup script for AdaptiveMasterMSM package

import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
		name='AdaptiveMasterMSM',
		version='0.0dev',
		description='Algorithms to do adaptive sampling exploiting MasterMSM',
		url='http://github.com/BioKT/AdaptiveMasterMSM',
		author='David De Sancho and Ion Mitxelena',
		author_email='daviddesancho.at.gmail.com',
		license='GPL',
        	packages=find_packages(),
		keywords= "adaptive sampling",
		long_description=read('README.md'),
		classifiers = ["""\
				Development Status :: 1 - Planning
				Operating System :: POSIX :: Linux
				Operating System :: MacOS
				Programming Language :: Python :: 3.7.6
				Topic :: Scientific/Engineering :: Bio-Informatics
				Topic :: Scientific/Engineering :: Chemistry
				"""]
		)

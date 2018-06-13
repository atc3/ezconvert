#!/usr/bin/env python3
# coding: utf-8

import io
import os

from os import path
from setuptools import setup, find_packages, Command

name = 'ezconvert'

setup(
  name=name,
  version='1.0.0',
  description='Convert, filter, transform tabular data',
  url='https://github.com/blahoink/ezconvert',
  author='Albert Chen',
  author_email='chen.alb@husky.neu.edu',
  license='MIT',
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Environment :: Console',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    
    'License :: OSI Approved :: MIT License',

    'Programming Language :: Python :: 3 :: Only',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',

    'Topic :: Scientific/Engineering :: Bio-Informatics'
  ],
  keywords='data conversion filtering transformation mapping maxquant',
  project_urls={
    'Documentation': 'https://github.com/blahoink/ezconvert',
    'Source': 'https://github.com/blahoink/ezconvert',
    'Tracker': 'https://github.com/blahoink/ezconvert/issues',
    'Lab Page': 'https://web.northeastern.edu/slavovlab/'
  },
  packages=['ezconvert'],
  #libraries=[],
  #setup_requires=[],
  install_requires=[
    'pyyaml>=3.12',
    'numpy>=1.14.3',
    'pandas>=0.22.0'
  ],
  #test_requires=[
  #  'pytest'
  #],
  extras_require={

  },
  package_data={
    'ezconvert': [
      'converters'
    ]
  },
  # data outside the package
  # data_files=[('my_data', ['data/data_file'])],
  
  entry_points={
    'console_scripts': [
      ('ezconvert=ezconvert.convert:convert')
    ]
  }

)
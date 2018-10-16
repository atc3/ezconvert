#!/usr/bin/env python3
# coding: utf-8

import csv
import numpy as np
import pandas as pd

## I/O configuration

# column delimiters for input and output files
input_sep = '\t'
output_sep = '\t'
output_type = '.txt'

# print row names/indices?
write_row_names=False

# print the column titles?
write_header=True

# quoting?
quoting=csv.QUOTE_NONNUMERIC

# leave empty to not print
additional_header = []

# separate by this value
sep_by = 'Raw file'


def __pep_001(df):
  # get PEP, ceil to 1
  pep = df['PEP'].values
  pep[pep > 1] = 1

  return (pep > 0.01)

filters = {
  'remove_decoy': (lambda df: df['Leading razor protein'].str.contains('REV__').values),
  'remove_contaminant': (lambda df: df['Leading razor protein'].str.contains('CON__').values),
  'pep_001': __pep_001
}

transformations = {
  'FirstScan': 'Best MS/MS',
  'Sequence': 'Sequence',
  'RT [min]': 'Retention time',
  'Modifications': 'Modifications'
}
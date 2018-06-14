#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import pandas as pd

## I/O configuration

# column delimiters for input and output files
input_sep = '\t'
output_sep = '\t'
output_type = 'txt'

# print row names/indices?
write_row_names=False

# print the column titles?
write_header=False

# leave empty to not print
additional_header = []

# separate by this value
sep_by = 'Raw file'

def __carrier_quant(df):
  carrier1 = df['Reporter intensity corrected 0']
  carrier2 = df['Reporter intensity corrected 1']

  return ((carrier1 == 0) | (carrier2 == 0))

def __sc_quant(df):
  sc = df['Reporter intensity corrected 4']
  return (sc == 0)

def __fdr_001(df):
  # get PEP, ceil to 1
  pep = df['PEP'].values
  pep[pep > 1] = 1

  # magic!!!
  # basically, we need to cumulatively sum the PEP to get the FDR
  # and then map back the cumulatively summed FDR to its original PEP
  qval = np.cumsum(pep[np.argsort(pep)])[np.argsort(pep).argsort()]

  return (qval < 0.01)

filters = {
  'remove_decoy': (lambda df: df['Leading razor protein'].str.contains('REV__').values),
  'remove_contaminant': (lambda df: df['Leading razor protein'].str.contains('CON__').values),
  'carrier_quant': __carrier_quant,
  'sc_quant': __sc_quant,
  'fdr_001': __fdr_001
}

def sc_to_carrier_ratio(df, df_out):
  sc = 'Reporter intensity corrected 4'
  carrier = 'Reporter intensity corrected 0'
  return df[sc] / df[carrier]

transformations = {
  'PSMId': 'id',
  'Sequence': 'Sequence',
  'RatioValue': sc_to_carrier_ratio,
  'RatioWeight': 'PIF',
  'Accession': (lambda df, df_out: df['Leading razor protein'].str.split('|').apply((
    lambda x: x[1] if len(x) is 3 else x[0])))
}
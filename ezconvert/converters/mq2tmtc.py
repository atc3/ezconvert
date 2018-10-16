#!/usr/bin/env python3
# coding: utf-8

import csv
import numpy as np
import pandas as pd

## I/O configuration

# column delimiters for input and output files
input_sep = '\t'
output_sep = ','
output_type = '_peptides.csv'

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

def __fdr_001(df):
  # get PEP, ceil to 1
  #pep = df['PEP'].values
  pep = df['PEP'].values
  pep[pep > 1] = 1

  # magic!!!
  # basically, we need to cumulatively sum the PEP to get the FDR
  # and then map back the cumulatively summed FDR to its original PEP
  qval = (np.cumsum(pep[np.argsort(pep)]) / np.arange(1, df.shape[0]+1))[np.argsort(np.argsort(pep))]

  return (qval > 0.01)


filters = {
  'remove_decoy': (lambda df: df['Proteins'].str.contains('REV__').values),
  'remove_contaminant': (lambda df: df['Proteins'].str.contains('CON__').values),
  'fdr_001': __fdr_001
}

transformations = {
  # Scan number - used to lookup in mzxml file
  'ScanF': 'Scan number',
  # not sure how ScanL is different.
  #'ScanL': 'Scan number',
  #'Time': 'Retention time',
  'PepID': 'Mod. peptide ID',
  #'MS2 ID': 'id',
  #'Obs M+H':
  #'Obs m/z':
  'z': 'Charge',
  #'Scan Event': 'Scan event number',
  #'Detector Type': 'Mass analyzer',
  #'intensity': 'Precursor Intensity',
  # convert Raw file to raw file ID
  #'SrchID': (lambda df, df_out: df['Raw file'].map({ \
  #    ind: val for val, ind in enumerate(np.sort(df['Raw file'].unique()))})),
  'SrchID': 'Raw file',
  # Pretty much same as SrchID, for now
  #'RunID': (lambda df, df_out: 100 + df['SrchID']),
  # Same as RunID
  #'ScansID': (lambda df, df_out: df['RunID']),
  #'SrchName': 'Raw file',
  #'File': 'Raw file', #??
  #'Rank': 1, # always 1??
  'Theo M+H': 'Mass',
  'Theo m/z': 'm/z',
  # maxquant mass error is still bugged with regard to TMT
  # IDK what is a good default for this. maybe 0? maybe the average of the mass
  # errors of the fragments in the MS/MS spectra? (we can calculate that)
  'PPM': 0,
  #'Dalton': 0,
  #'Xcorr': 'Score',
  #'∆corr': 'Delta score',
  #'Uniq. ∆corr': 'Delta score',
  #
  # Get the leading protein
  'Reference': (lambda df, df_out: df['Proteins'].str.split(';').apply((
    lambda x: x[0] if type(x) == list else x))),
  'Peptide': 'Sequence',
  'LDA Probability': 'PEP'
}
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
quoting=csv.QUOTE_MINIMAL

# leave empty to not print
additional_header = [
  # empty for specid, label, scannr, expmass, calcmass, 
  # and doc features (RT and mass calibration) - 7 blanks total
  'DefaultDirection', '-', '-', '-', '-', '-', '-',
  '0', # mass
  '1', # score
  '1.5', # delta score
  '-0.573', # peptide length
  '0.0335', '0.149', '-0.156', # charge states
  '0', '0' # enzymatic features
  # skip peptide and protein
  '', ''
]


def __calibration_dat_filter(df):
  return pd.isnull(df['Uncalibrated - Calibrated m/z [Da]'])

filters = {
  'remove_no_calibration_data': __calibration_dat_filter
}

def __label(df, df_out):
  # target or decoy
  label = df['Leading razor protein'].str.contains('REV__').values.astype(int)
  label[label==1] = -1
  label[label==0] = 1

  return label

def __scan_num(df, df_out):
  # adjust scan numbers so that the scan numbers from different
  # experiments don't overlap
  max_scan_nums = df.groupby('Raw file')['MS/MS scan number'].apply((lambda x: np.max(x)))
  max_scan_nums = max_scan_nums[np.argsort(max_scan_nums.index.values)]
  max_scan_nums = np.cumsum(max_scan_nums) - max_scan_nums[0]

  scannr = df['MS/MS scan number'] + df['Raw file'].map(max_scan_nums)

  return scannr

def __protein(df, df_out):
  # extract UniProt IDs from protein string
  prots = df['Leading razor protein'].str.split('|').apply((
    lambda x: x[1] if len(x) is 3 else x[0]))
  # mark reverse ones differently still
  #prots.loc[df_out['Label'] == -1] = ('REV_' + prots.loc[df_out['Label'] ==- 1])

  return prots


transformations = {
  # Unique ID for each MS2 spectra
  'SpecId': 'id',
  # Target or Decoy
  'Label': __label,
  # MS/MS scan number
  'ScanNr': __scan_num,
  # ExpMass - measured mass of precursor ion (measured m/z * charge)
  'ExpMass': (lambda df, df_out: df['m/z'] * df['Charge']),
  # CalcMass - theoretical mass of precursor ion
  'CalcMass': 'Mass',

  # --doc features, RT and mass calibration
  'RT':              'Retention time',
  'MassCalibration': 'Uncalibrated - Calibrated m/z [Da]',

  # precursor properties
  'Mass': 'Mass',
  # Score - most important feature
  'Score': 'Score',
  'DeltaScore': 'Delta score',

  # sequence characteristics
  # peptide length
  'PepLen': 'Length',
  # charge states
  'Charge1': (lambda df, df_out: (df['Charge'] == 1).values.astype(int)),
  'Charge2': (lambda df, df_out: (df['Charge'] == 2).values.astype(int)),
  'Charge3': (lambda df, df_out: (df['Charge'] == 3).values.astype(int)),

  # enzymatic performance - trypsin
  # don't have data on n-terminus cleavage...
  'enzC': (lambda df, df_out: df['Sequence'].str.slice(-1).isin(['R', 'K']).values.astype(int)),
  # missed cleavages = number of enzymatic sites 
  'enzInt': 'Missed cleavages',

  # lastly, the peptide and protein
  # need to surround peptide with flanking amino acids
  # we don't have this data, so just append alanines on both sides
  'Peptide': (lambda df, df_out: 'A.' + df['Sequence'] + '.A'),

  # extract UniProt IDs from protein string
  'Protein': __protein
}
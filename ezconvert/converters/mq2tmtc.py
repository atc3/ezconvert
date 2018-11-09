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
#quoting=csv.QUOTE_NONNUMERIC
quoting=csv.QUOTE_NONE

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

# filter out observations w/o any mass error value
def __missing_mass_error(df):
  mass_error = df['Mass error [ppm]'].values
  simple_mass_error = df['Simple mass error [ppm]'].values

  return (pd.isnull(mass_error) & pd.isnull(simple_mass_error))

# filter out mass error above 20 PPM
def __large_mass_error(df):
  mass_error = df['Mass error [ppm]'].values
  simple_mass_error = df['Simple mass error [ppm]'].values

  nan_inds = pd.isnull(mass_error)
  mass_error[nan_inds] = simple_mass_error[nan_inds]

  return np.abs(mass_error) > 20


filters = {
  'remove_decoy': (lambda df: df['Proteins'].str.contains('REV__').values),
  'remove_contaminant': (lambda df: df['Proteins'].str.contains('CON__').values),
  'remove_no_protein': (lambda df: pd.isnull(df['Proteins'])),
  'fdr_001': __fdr_001,
  'missing_mass_error': __missing_mass_error,
  'large_mass_error': __large_mass_error
}

def __mass_error_correction(df, df_out):
  # if Mass error [ppm] is Nan, then replace w simple mass error [ppm]
  # drop all that have ppm > -20 as these clearly got wrong peak
  mass_error = df['Mass error [ppm]'].values
  simple_mass_error = df['Simple mass error [ppm]'].values

  nan_inds = pd.isnull(mass_error)
  mass_error[nan_inds] = simple_mass_error[nan_inds]

  return mass_error

transformations = {
  # Scan number - used to lookup in mzxml file
  'ScanF': 'Scan number',
  # not sure how ScanL is different.
  #'ScanL': 'Scan number',
  'Time': 'Retention time',
  #'PepID': 'Mod. peptide ID',
  'PepID': 'id',
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
  'PPM': __mass_error_correction,
  #'Dalton': 0,
  #'Xcorr': 'Score',
  #'∆corr': 'Delta score',
  #'Uniq. ∆corr': 'Delta score',
  #
  # Preliminary scores - from engines like Sequest
  # 
  #'SP': 0,
  #'Rank/SP': 0,
  #
  # Not sure what these are used for
  # 
  #'# Ions': 0,
  #'# Ions Link': 0,
  #'Ions Matched': 0,
  #'Total Ions': 0,
  # Get the leading protein
  #'Reference': (lambda df, df_out: df['Proteins'].str.split(';').apply((
  #  lambda x: x[0] if type(x) == list else x))),
  'Reference': 'Proteins',
  #'Reference Link': 0,
  #'raw_pep_per_ref': 0,
  #'raw_psm_per_ref': 0,
  #'Redun': 0,
  #'Redun Link': 0,
  'Peptide': 'Sequence',
  #'Modified sequence': 'Modified sequence',
  #'Reverse': 'Reverse',
  #'Peptide Link': 0,
  #'Trimmed Peptide': 'Sequence',
  #'Pept. Length': 'Length',
  #'Trypticity': 0, #fully-tryptic, or something else?
  #'MissedCleav': 0,
  #'Validity': 0, #63?
  #'Solu. z': 0, #2 or 3
  #'Gene Symbol': 0,
  #'Protein MWT(kDa)': 0,
  #'Annotation': 0,
  'LDA Probability': 'PEP'
}
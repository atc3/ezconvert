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
  'remove_acetyl_ox_modifications': (lambda df: df['Modifications'].str.contains('Acetyl|Oxidation').values),
  'fdr_001': __fdr_001,
  'missing_mass_error': __missing_mass_error,
  'large_mass_error': __large_mass_error,
  #'remove_phospho': (lambda df: df['Modified sequence'].str.contains('p').values)
  'only_phospho': (lambda df: ~df['Modified sequence'].str.contains('p').values)
}

# element monoisotopic masses
element_mass = {
    'H': 1.0078250321,
    'C': 12.0,
    'Cx': 13.0033548378,
    'N': 14.0030740052,
    'Nx': 15.0001088984,
    'O': 15.9949146221,
    'F': 18.99840320,
    'P': 30.97376151,
    'S': 31.97207069,
    'Cl': 34.96885271,
    'Na': 22.98976967,
    'K': 38.9637069,
    'Ca': 39.9625912,
    'Fe': 53.9396148
}

# amino acid monoisotopic masses
aa_mass = {
    'A': 71.0371137878,
    'C': 103.00918447779999,
    'D': 115.026943032,
    'E': 129.0425930962,
    'F': 147.0684139162,
    'G': 57.0214637236,
    'H': 137.0589118624,
    'I': 113.0840639804,
    'K': 128.0949630177,
    'L': 113.0840639804,
    'M': 131.0404846062,
    'N': 114.0429274472,
    'P': 97.052763852,
    'Q': 128.0585775114,
    'R': 156.1011110281,
    'S': 87.03202840990001,
    'T': 101.04767847410001,
    'U': 150.95363, # Selenocysteine
    'V': 99.0684139162,
    'W': 186.0793129535,
    'Y': 163.0633285383,
    'X': 113.0840639804 # Xle - Isoleucine or Leucine
}
proton_mass = 1.0072764666
water_mass = 18.0105646863
tmt_mass = 229.162932141

def __peptide_to_mass(seq):
  mass = 0

  # add up amino acid masses
  for i in range(0, len(seq)):
    if seq[i] == '_': continue
    elif seq[i] == 'p': mass = mass + 79.9663304084
    else: mass = mass + aa_mass[seq[i]]

  # add TMT. one for every lysine (K), and one for the n-terminus
  mass = mass + ((seq.count('K') + 1) * tmt_mass)

  # add water
  mass = mass + water_mass

  # fake a neutral loss
  # NL_H3P04 = 97.9768950947
  # mass = mass - NL_H3P04

  return mass

def __predict_mass(df, df_out):
  return df['Modified sequence'].apply(__peptide_to_mass)

def __predict_m_plus_h(df, df_out):
  mass = __predict_mass(df, df_out)
  return mass + (df['Charge'] * proton_mass)

def __predict_mz(df, df_out):
  mass = __predict_m_plus_h(df, df_out)
  return mass / df['Charge']

def __mass_error_correction(df, df_out):
  # if Mass error [ppm] is Nan, then replace w simple mass error [ppm]
  # drop all that have ppm > -20 as these clearly got wrong peak
  mass_error = df['Mass error [ppm]'].values
  simple_mass_error = df['Simple mass error [ppm]'].values

  nan_inds = pd.isnull(mass_error)
  mass_error[nan_inds] = simple_mass_error[nan_inds]

  return mass_error

# duplicate each row three times, and for each duplicate
# subtract a possible phospho neutral loss
def __neutral_losses(df, df_out):
  NL_H3P04 = 97.9768950947
  NL_HPO3  = 79.9663304084
  NL_H5PO5 = 115.987459781

  df_nl_H3PO4 = df_out.copy()
  df_nl_H3PO4['Theo M+H'] = __predict_m_plus_h(df, df_out) - NL_H3P04
  df_nl_H3PO4['Theo m/z'] = df_nl_H3PO4['Theo M+H'] / df_nl_H3PO4['z']

  df_nl_HPO3 = df_out.copy()
  df_nl_HPO3['Theo M+H'] = __predict_m_plus_h(df, df_out) - NL_HPO3
  df_nl_HPO3['Theo m/z'] = df_nl_HPO3['Theo M+H'] / df_nl_HPO3['z']

  df_nl_H5PO5 = df_out.copy()
  df_nl_H5PO5['Theo M+H'] = __predict_m_plus_h(df, df_out) - NL_H5PO5
  df_nl_H5PO5['Theo m/z'] = df_nl_H5PO5['Theo M+H'] / df_nl_H5PO5['z']

  # append all dataframes by row and return
  df_out = df_out.append(df_nl_H3PO4)
  df_out = df_out.append(df_nl_HPO3)
  df_out = df_out.append(df_nl_H5PO5)

  return df_out

transformations = {
  'ScanF': 'Scan number',
  'Reference': 'Proteins',
  'Peptide': 'Sequence',
  'Modified sequence': 'Modified sequence',
  'SrchID': 'Raw file',
  'Theo m/z': __predict_mz,
  'Theo M+H': __predict_m_plus_h,
  'LDA Probability': 'PEP',
  'PepID': 'id',
  'PPM': __mass_error_correction,
  'z': 'Charge'
  # '__neutral_losses': __neutral_losses
  # Scan number - used to lookup in mzxml file
  #'ScanF': 'Scan number',
  # not sure how ScanL is different.
  #'ScanL': 'Scan number',
  #'Time': 'Retention time',
  #'PepID': 'Mod. peptide ID',
  #'PepID': 'id',
  #'MS2 ID': 'id',
  #'Obs M+H':
  #'Obs m/z':
  #'z': 'Charge',
  #'Scan Event': 'Scan event number',
  #'Detector Type': 'Mass analyzer',
  #'intensity': 'Precursor Intensity',
  # convert Raw file to raw file ID
  #'SrchID': (lambda df, df_out: df['Raw file'].map({ \
  #    ind: val for val, ind in enumerate(np.sort(df['Raw file'].unique()))})),
  #'SrchID': 'Raw file',
  # Pretty much same as SrchID, for now
  #'RunID': (lambda df, df_out: 100 + df['SrchID']),
  # Same as RunID
  #'ScansID': (lambda df, df_out: df['RunID']),
  #'SrchName': 'Raw file',
  #'File': 'Raw file', #??
  #'Rank': 1, # always 1??
  #'Theo M+H': 'Mass',
  #'Theo m/z': 'm/z',
  # maxquant mass error is still bugged with regard to TMT
  # IDK what is a good default for this. maybe 0? maybe the average of the mass
  # errors of the fragments in the MS/MS spectra? (we can calculate that)
  #'PPM': __mass_error_correction,
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
  #'Reference': 'Proteins',
  #'Reference Link': 0,
  #'raw_pep_per_ref': 0,
  #'raw_psm_per_ref': 0,
  #'Redun': 0,
  #'Redun Link': 0,
  #'Peptide': 'Sequence',
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
  #'LDA Probability': 'PEP'
}
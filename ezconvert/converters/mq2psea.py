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

def __sqc_set(df):
  return ~df['Raw file'].str.contains('SQC')

def __sc_quant(df):
  dcols = df.columns[df.columns.str.contains('Reporter intensity corrected')]
  return df[dcols[4:10]].apply((lambda x: x == 0)).apply(np.sum, 1) > 0

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
  'remove_decoy': (lambda df: df['Leading razor protein'].str.contains('REV__').values),
  'remove_contaminant': (lambda df: df['Leading razor protein'].str.contains('CON__').values),
  #'sqc_set': __sqc_set,
  'sc_quant': __sc_quant,
  'fdr_001': __fdr_001
}

def __sc_ratios(df, df_out):
  # J - 128C, 129C, 130C
  # U - 129N, 130N, 131N
  
  # extract uniprot accession number
  Proteins = df['Leading razor protein'].str.split("|").apply((lambda x: x[1] if len(x) == 3 else x[0]))

  dcols = df.columns[df.columns.str.contains('Reporter intensity corrected')]
  dmat = df[dcols[4:10]].values

  # normalize by column and then row
  dmatn = dmat / np.apply_along_axis(np.mean, 0, dmat)
  dmatn = (dmatn / np.apply_along_axis(np.mean, 1, dmatn)[:,None])
  
  # get j/u ratios
  j_channels = [0, 2, 4]
  u_channels = [1, 3, 5]
  ratio_mat = np.zeros((dmatn.shape[0], len(j_channels) * len(u_channels)))
  for i, j in enumerate(j_channels):
    for k, u in enumerate(u_channels):
      ratio_mat[:,(i*len(j_channels))+k] = (dmatn[:,j] / dmatn[:,u])

  # collapse ratios by mean for each protein ID
  prot_list = np.sort(Proteins.unique())
  prot_ratios = np.zeros((len(prot_list), len(j_channels) * len(u_channels)))
  for i, p in enumerate(prot_list):
    p_inds = (Proteins == p)
    prot_ratios[i,:] = np.apply_along_axis(np.mean, 0, ratio_mat[p_inds,:])
  
  df_a = pd.DataFrame(prot_ratios)
  df_a = pd.concat([pd.Series(prot_list), df_a], axis=1)

  return df_a

transformations = {
  #'ProteinSymbol': (lambda df, df_out: df['Leading razor protein'].str.split('|').apply((
  #  lambda x: x[1] if len(x) is 3 else x[0])))
  '__sc_ratios': __sc_ratios
}
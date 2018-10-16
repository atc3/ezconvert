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
write_header=False

# quoting?
quoting=csv.QUOTE_MINIMAL

# leave empty to not print
additional_header = []

def __sqc_set(df):
  return ((~df['Raw file'].str.contains('SQC')) | df['Raw file'].str.contains('SQC9'))

def __sc_quant(df):
  dcols = df.columns[df.columns.str.contains('Reporter intensity corrected')]
  #return df[dcols[4:10]].apply((lambda x: x == 0)).apply(np.sum, 1) > 0
  return (np.apply_along_axis(np.sum, 1, df[dcols[4:10]].values == 0) > 0)

def __fdr_001(df):
  # get PEP, ceil to 1
  pep = df['PEP'].values
  pep[pep > 1] = 1

  # magic!!!
  # basically, we need to cumulatively sum the PEP to get the FDR
  # and then map back the cumulatively summed FDR to its original PEP
  qval = (np.cumsum(pep[np.argsort(pep)]) / np.arange(1, df.shape[0]+1))[np.argsort(np.argsort(pep))]

  return (qval > 0.01)

def __new_fdr_001(df):
  # get PEP, ceil to 1
  pep = df['pep_updated'].values
  pep[pep > 1] = 1

  qval = (np.cumsum(pep[np.argsort(pep)]) / np.arange(1, df.shape[0]+1))[np.argsort(np.argsort(pep))]

  return (qval > 0.01)

def __prot_fdr_001(df):
  # filter at protein FDR of 1%
  return (pd.isnull(df['prot_fdr']) | (df['prot_fdr'] > 0.01))

# only select de novo proteins
def __de_novo_proteins(df):
  # extract uniprot accession number
  Proteins = df['Leading razor protein'].str.split("|").apply((lambda x: x[1] if len(x) == 3 else x[0]))
  # by not matching a dash we are ignoring isoforms
  Proteins = Proteins.str.extract('([A-Z0-9_]+)')
  # or, extract gene name (protein symbol)
  #Proteins = df['Leading razor protein'].str.extract(r'([A-Z0-9-]+)_HUMAN', expand=False)

  # grab PEPs
  pep = df['PEP']
  #pep[pep > 1] = 1

  df_a = pd.concat([Proteins, pep], axis=1)
  df_a.columns = ['protein', 'pep']

  prots = Proteins.unique()
  prots = prots[~pd.isnull(prots)]
  prots = prots[(df_a.groupby('protein')['pep'].min() < 0.01)]

  return Proteins.isin(prots)

  """
  # init exclusion vector
  exclude = np.repeat(False, df.shape[0])

  # for each protein ID:
  for p in np.sort(Proteins.unique()):
    p_inds = (Proteins == p)
    # if there is 1 observation that has PEP < 0.01,
    # then exclude all PSMs belonging to this protein
    if np.any(pep[p_inds] < 0.01):
      exclude = (exclude | p_inds)

  return exclude
  """
def __de_novo_peptides(df):
  Peptides = df['Modified sequence'].copy()
  pep = df['PEP'].copy()
  pep[pep > 1] = 1

  df_a = pd.concat([Peptides, pep], axis=1)
  df_a.columns.values = ['peptide', 'pep']

  peptides = Peptides.unique()
  peptides = peptides[~pd.isnull(peptides)]
  peptides = peptides[(df_a.groupby('peptide')['pep'].min() < 0.01)]

  return Peptides.isin(peptides)

filters = {
  'remove_decoy': (lambda df: df['Leading razor protein'].str.contains('REV__').values),
  'remove_contaminant': (lambda df: df['Leading razor protein'].str.contains('CON__').values),
  'sqc_set': __sqc_set,
  'sc_quant': __sc_quant,
  'prot_fdr': __prot_fdr_001,
  #'de_novo_proteins': __de_novo_proteins,
  #'de_novo_peptides': __de_novo_peptides,
  #'new_fdr_001': __new_fdr_001
  'fdr_001': __fdr_001
}

def __sc_ratios(df, df_out):
  # J - 128C, 129C, 130C
  # U - 129N, 130N, 131N
  
  # extract uniprot accession number
  Proteins = df['Leading razor protein'].str.split('|').apply((lambda x: x[1] if len(x) == 3 else x[0]))
  # by not matching a dash we are ignoring isoforms
  Proteins = Proteins.str.extract('([A-Z0-9_]+)')
  # or, extract gene name (protein symbol)
  #Proteins = df['Leading razor protein'].str.extract(r'([A-Z0-9-]+)_HUMAN', expand=False)

  dcols = df.columns[df.columns.str.contains('Reporter intensity corrected')]
  dmat = df[dcols[4:10]].values

  # normalize by column and then row
  dmatn = dmat / np.apply_along_axis(np.mean, 0, dmat)
  dmatn = (dmatn / np.apply_along_axis(np.mean, 1, dmatn)[:, None])
  
  # get j/u ratios
  j_channels = [0, 2, 4]
  u_channels = [1, 3, 5]
  ratio_mat = np.zeros((dmatn.shape[0], len(j_channels) * len(u_channels)))
  for i, j in enumerate(j_channels):
    for k, u in enumerate(u_channels):
      ratio_mat[:,(i*len(j_channels))+k] = (dmatn[:,j] / dmatn[:,u])

  # collapse ratios by mean for each protein ID
  prot_list = Proteins.unique()
  prot_ratios = np.zeros((len(prot_list), len(j_channels) * len(u_channels)))
  for i, p in enumerate(prot_list):
    if pd.isnull(p): continue
    #p_inds = Proteins.str.contains(p).values
    p_inds = (Proteins == p)
    prot_ratios[i,:] = np.apply_along_axis(np.mean, 0, ratio_mat[p_inds,:])
  
  df_a = pd.DataFrame(prot_ratios)

  # shuffle prot list for a null control?
  np.random.shuffle(prot_list)

  df_a = pd.concat([pd.Series(prot_list), df_a], axis=1)

  return df_a

transformations = {
  #'ProteinSymbol': (lambda df, df_out: df['Leading razor protein'].str.split('|').apply((
  #  lambda x: x[1] if len(x) is 3 else x[0])))
  '__sc_ratios': __sc_ratios
}
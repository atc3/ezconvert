#!/usr/bin/env python3
# coding: utf-8

import argparse
import logging
import numpy as np
import os
import pandas as pd
import sys
import yaml

logger = logging.getLogger('root')
  
def convert():
  # load command-line args
  parser = argparse.ArgumentParser()  

  parser.add_argument('-v', '--verbose', action='store_true', default=False,
    help='Run in verbose mode. If piping output from stdout to a file, leave this off to exclude all logging messages.')

  parser.add_argument('--config-file', required=True, 
    type=argparse.FileType('r', encoding='UTF-8'), 
    help='Path to conversion configuration script. See list of converters in converters/ folder')
  input_group = parser.add_mutually_exclusive_group(required=True)
  input_group.add_argument('--input-list', type=argparse.FileType('r', encoding='UTF-8'),
    help='List of input files, in YAML format.')
  input_group.add_argument('-i', '--input', type=argparse.FileType('r', encoding='UTF-8'),
    nargs='+', help='List of input files, separated by spaces.')

  parser.add_argument('-o', '--output', type=str, 
    help='Path to output data. Default: Leave empty to print to stdout')

  args = parser.parse_args()

  # initialize logger
  # set up logger
  for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler) 
   
  logFormatter = logging.Formatter('%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s')
  logger = logging.getLogger('root')

  if args.verbose: logger.setLevel(logging.DEBUG)
  else: logger.setLevel(logging.WARNING)

  """
  if log_to_file:
    fileHandler = logging.FileHandler(log_file_path, mode='w')
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)
  """

  consoleHandler = logging.StreamHandler()
  consoleHandler.setFormatter(logFormatter)
  logger.addHandler(consoleHandler)
  logger.info(' '.join(sys.argv[0:]))

  # load vars from the config file
  logger.info('Loading config file functions from {}.'.format(args.config_file.name))
  exec(compile(open(args.config_file.name, 'rb').read(), args.config_file.name, 'exec'), 
    globals())

  # read inputs, either from the input list or from the command line
  _input = []
  if args.input_list is not None:
    logger.info('Reading in input files from input list {}.'.format(args.input_list.name))
    with open(args.input_list.name, 'r') as f:
      _input = yaml.load(f)
  else:
    logger.info('Reading in input files from command line.')
    _input = [f.name for f in args.input]

  if len(_input) == 0:
    raise Exception('No input files provided, either from the input list or the command line.')


  df = pd.DataFrame()

  # iterate through each input file provided.
  for i, f in enumerate(_input):
    # first expand user or any vars
    f = os.path.expanduser(f)
    f = os.path.expandvars(f)

    logger.info('Reading in input file #{} | {} ...'.format(i+1, f))

    dfa = pd.read_csv(f, sep=input_sep, low_memory=False)

  df = df.append(dfa)

  # filter observations
  logger.info('\nFiltering observations...')

  # before we filter, assign every row an ID
  df['id'] = range(0, df.shape[0])

  # by default, exclude nothing. we'll use binary ORs (|) to
  # gradually add more and more observations to this exclude blacklist
  df['exclude'] = np.repeat(False, df.shape[0])

  # run all the filters specified by the list in the input config file
  # all filter functions are passed df, and the run configuration
  # after each filter, append it onto the exclusion master list with a bitwise OR
  # if the filter function returns None, then just ignore it.
  for i, f in enumerate(filters):
    logger.info('Applying filter #{}: \"{}\"'.format(i+1, f))
    e = filters[f](df)
    if e is not None:
      df['exclude'] = (df['exclude'] | e)

  # apply exclusion filter
  df = df[~df['exclude']].reset_index(drop=True)

  # create output frame
  df_out = pd.DataFrame()

  # apply transformations
  logger.info('\nTransforming data...')

  for i, t in enumerate(transformations):
    logger.info('Applying transformation #{}: \"{}\"'.format(i+1, t))
    trans = transformations[t]
    if type(trans) is str:
      df_out[t] = df[trans]
    elif callable(trans):
      df_out[t] = trans(df, df_out)
    else:
      raise Exception('Invalid transformation type: {}. Please provide either a string or a function'.format(type(trans)))

  # write headers and weights
  
  headers = ''

  # column headers
  for i, col in enumerate(df_out.columns):
    headers += col
    if i != (len(df_out.columns)-1):
      headers += output_sep
  headers += '\n'

  # additional header
  if 'additional_headers' in globals() and len(additional_headers) > 0:
    for i, w in enumerate(additional_headers):
      headers += w
      if i != (len(weights)-1):
        headers += output_sep
    headers += '\n'

  if args.output is None:
    # if none, then just print to stdout
    print(headers, end='')
    print(df_out.to_string(header=False, index=False, sparsify=False))
  else:
    with open(args.output, 'w') as f:
      f.write(headers)
    logger.info('Writing output to {} ...'.format(args.output))
    df_out.to_csv(args.output, sep=output_sep, header=False, index=False, mode='a')

if __name__ == '__main__':
  convert()
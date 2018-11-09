#!/usr/bin/env python3
# coding: utf-8

import argparse
import csv
import logging
import numpy as np
import os
import pandas as pd
import pkg_resources
import sys
import yaml

logger = logging.getLogger('root')

provided_converters = [
  'mq2pin',
  'mq2pcq',
  'mq2psea',
  'mq2elutator_trainer',
  'mq2tmtc'
]

def write_df_to_file(df, headers, out_path):
  with open(out_path, 'w') as f:
    f.write(headers)
  logger.info('Writing output to {} ...'.format(out_path))
  df.to_csv(out_path, sep=output_sep, header=False, 
    index=write_row_names, mode='a', quoting=quoting)
  
def convert():
  # load command-line args
  parser = argparse.ArgumentParser()  

  parser.add_argument('-v', '--verbose', action='store_true', default=False,
    help='Run in verbose mode. If piping output from stdout to a file, leave this off to exclude all logging messages.')

  parser.add_argument('--config-file', required=True, type=str, 
    help='One of these converters: [' + ' '.join(provided_converters) + '], or a path to conversion configuration script. See list of converters in converters/ folder')
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
  config_file = ''
  if args.config_file in provided_converters:
    config_file = pkg_resources.resource_string('ezconvert', '/'.join(('converters', args.config_file + '.py')))
  else:
    logger.info('Loading config file functions from {}.'.format(args.config_file))
    with open(args.config_file, 'rb') as f:
      config_file = f.read()

  exec(compile(config_file, args.config_file, 'exec'), globals())

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

    logger.info('Read {} PSMs'.format(dfa.shape[0]))

    # track input file with input id
    dfa['input_id'] = i

  df = df.append(dfa)

  # filter observations
  logger.info('Filtering observations...')

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
    logger.info('Marked {} observations for filtering'.format(np.sum(e)))
    if e is not None:
      df['exclude'] = (df['exclude'] | e)

  logger.info('{} / {} ({:.2%}) observations pass filters'.format(df.shape[0] - df['exclude'].sum(), df.shape[0], (df.shape[0] - df['exclude'].sum()) / df.shape[0]))

  # apply exclusion filter
  df = df[~df['exclude']].reset_index(drop=True)

  # create output frame
  df_out = pd.DataFrame()

  # apply transformations
  logger.info('Transforming data...')

  for i, t in enumerate(transformations):
    logger.info('Applying transformation #{}: \"{}\"'.format(i+1, t))
    trans = transformations[t]
    # if transformation is a string, then simply copy the old column
    # to the new output one
    if type(trans) is str:
      df_out[t] = df[trans]
    # if transformation is a function, and the transformation name
    # begins with a '__', then apply the transformation over the
    # entire data frame.
    # this is useful for doing matrix maths that spans over multiple
    # columns, or rows, or something that involves more than just one
    # column.
    elif callable(trans) and t[0:2] == '__':
      df_out = trans(df, df_out)
    # if transformation is a function, then call that function to
    # generate the new column for the output
    elif callable(trans):
      df_out[t] = trans(df, df_out)
    # if transformation is a constant number, then just set all values
    # of that name to the specified number
    # don't have to vectorize, pandas will handle that.
    elif type(trans) is int or type(trans) is float:
      df_out[t] = trans
    else:
      raise Exception('Invalid transformation type: {}. Please provide either a string or a function'.format(type(trans)))

  # write headers and weights
  
  headers = ''

  # column headers
  if write_header:
    for i, col in enumerate(df_out.columns):
      if quoting == csv.QUOTE_ALL or quoting == csv.QUOTE_NONNUMERIC:
        headers += ("\"" + col + "\"")
      else: 
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
    print(df_out.to_string(header=False, index=write_row_names, sparsify=False))
  else:
    if 'sep_by' in globals() and type(sep_by) is str:
      # separate output files based on a certain column
      if sep_by in df_out.columns:
        sep_by_vals = df_out[sep_by]
      elif sep_by in df.columns:
        sep_by_vals = df[sep_by]
      else:
        raise Exception('File separator not found in the columns of either the input file or the transformed output file.') 

      # create the output path if necessary
      if not os.path.exists(args.output):
        logger.info('Path for output folder {} does not exist. Creating...'.format(args.output))
        os.makedirs(args.output)

      # get unique categories
      cats = np.unique(sep_by_vals)
      logger.info('Splitting observations into separate files by "' + sep_by + '"')
      logger.info('Categories: [' + ' '.join(cats) + ']')
      # iterate over each category
      for c in cats:
        out_path = os.path.join(args.output, '{}{}'.format(c, output_type))
        logger.info('Saving category file {} to {}'.format(c, out_path))
        df_a = df_out.loc[sep_by_vals == c]
        write_df_to_file(df_a, headers, out_path)

    else:
      # if no separation, then write the entire collated df to file
      logger.info('Saving combined file to {}'.format(args.output))
      write_df_to_file(df_out, headers, args.output)


  logger.info('Done!')

if __name__ == '__main__':
  convert()
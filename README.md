# ezconvert
Convert, filter, and transform tabular data with a simple configuration file

## Introduction

ezconvert is a simple script that combines tabular input files, filters them, and then transforms values.

Easily program your own converter, and then use it from the command line.

Output options are a bit sparse right now. Please suggest any additions in the issues page, or submit a pull request.

## Installation

```
pip install git+https://github.com/blahoink/ezconvert
```

## Usage

```
usage: ezconvert [-h] [-v] --config-file CONFIG_FILE
                  (--input-list INPUT_LIST | -i INPUT [INPUT ...]) [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Run in verbose mode. If piping output from stdout to a
                        file, leave this off to exclude all logging messages.
  --config-file CONFIG_FILE
                        One of these converters: [mq2pin], or a path to
                        conversion configuration script. See list of
                        converters in converters/ folder
  --input-list INPUT_LIST
                        List of input files, in YAML format.
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        List of input files, separated by spaces.
  -o OUTPUT, --output OUTPUT
                        Path to output data. Default: Leave empty to print to
                        stdout
```

Example:

```
ezconvert --config-file mq2pin --input-list ./input_list.yaml
```

```
ezconvert --config-file path/to/converter.py --input-list ./input_list.yaml
```

Clone the repo to get the list of existing converters, or write your own.

### Input List:

The input list is a YAML file, with items in a list like so:

```
- /path/to/input/file/1.csv
- /path/to/input/file/2.csv
- /path/to/input/file/3.csv
```

```
ezconvert --config-file path/to/config/file --input-list path/to/input_list.yaml
```

You can also list these files in the command-line with the "-i" option

```
ezconvert --config-file path/to/config/file -i /path/to/input/file/1.csv /path/to/input/file/2.csv /path/to/input/file/3.csv
```

## Converters

Converters are defined as python scripts, to give this type of configuration file functionality and flexibility

### I/O Configuration
- ```input_sep```: delimiter for the input file
- ```output_sep```: delimiter for the output file
- ```additional_header```: additional string, or list of items to be separated by the output delimiter. This is printed after the column name headers, but before the data.

### Filters

Filters are listed in a dictionary called ```filters```. The key name is arbitrary, but the value is a function which is passed the input data frame. For example:

```
def __calibration_dat_filter(df):
  return pd.isnull(df['Calibration'])

filters = {
  'remove_no_calibration_data': __calibration_dat_filter
}
```

### Transformations

Transformations are listed in a dictionary called ```transformations```. The key name is the column name for the output. The value can either be a string or a function. 

A string value is the name of the column from the input file, and all values will be copied to the output. 

A function is a mapping function, which is passed two variables: ```df```, the input data, and ```df_out```, the output data (as it is currently being built, down the ```transformations``` dictionary).

```

def __label(input, output):
  # target or decoy?
  label = input['Protein'].str.contains('DECOY__').values.astype(int)
  label[label==1] = -1
  label[label==0] = 1

  return label

# values can be string, function, or lambda
transformations = {
  # a string copies values from the input
  'SpecId': 'id',
  # specify a defined function, with input and output as the two arguments
  'Label': __label,
  # specify a lambda function, with input and output as the two arguments
  'ExpMass': (lambda input, output: df['m/z'] * df['Charge']),
  # reference output columns, as long as they are listed before this item
  # (items executed sequentially)
  'ExpMassShift': (lambda input, output: output['ExpMass'] + input['MassShift'])
}
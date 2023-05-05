# FETCH
Pipeline for processing flow cytometry .fcs files to get FETCH score

## Installation: 
`conda env create -n FETCH --file FETCH_env.yml`

`conda activate FETCH`

## Usage (to run an example):
`python FETCH.py -f example -p my_cool_project -s FC114_A2_A02_002.fcs FC114_C1_C01_025.fcs`

## Parameters:

- `-f` or `--folder`: The path to a directory containing your .fcs files
- `-p` or `--project`: The project name
- `-s` or `--skip_renaming`: This script expects fcs files to be formatted as 'FC114_A1_A01_001.fcs'; if just a few files are not formatted like that, pass them as arguments in this function to avoid errors; otherwise edit summarize function too parse your filenames correctly

# FETCH
Pipeline for processing flow cytometry .fcs files to get FETCH score

Installation: 
conda env create -n FETCH --file FETCH_env.yml python=3.8.8
conda activate FETCH

Usage (to run an example):
python FETCH.py -f example -p my_cool_project -s FC114_A2_A02_002.fcs FC114_C1_C01_025.fcs

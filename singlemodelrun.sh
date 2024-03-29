#!/bin/sh

# This script runs the main py script for the ice model for a single buoy.
# edit the buoy no and preprocess options to run the scripts.


alias cd='cd py_script_files'
cd
# changing the working directory to py_script_files.

# buoy variable gives the buoy number for which the model is run
# The options for buoy are 02, 03, 09, 13, 14, 16
export buoy=16

# preprocessing variable tells whether to do preprocessing or not
#1 is for preprocessing and any other value is for not doing preprocessing.
export preprocess=0

python main_script.py $buoy $preprocess

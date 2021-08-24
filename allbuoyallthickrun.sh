#!/bin/sh

# This is a trial script to run through all the buoys.  


alias cd='cd py_script_files'


# changing the working directory to py_script_files.
cd

# buoy variable gives the buoy number for which the model is run
# The options for buoy are 02, 03, 09, 13, 14, 16
#export buoy=16

# preprocessing variable tells whether to do preprocessing or not
#1 is for preprocessing and any other value is for not doing preprocessing.
export preprocess=0

for i in 02 03 09 13 14 16 
do
  echo "Running for... i is set to $i."
  buoy=$i  
  python allthickrun.py $buoy 
done


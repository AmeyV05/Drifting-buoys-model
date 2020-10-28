@echo off 
:: This batch file runs the script for the ice model.

:: buoy variable gives the buoy number for which the model is run
:: The options for buoy are 02, 03, 09, 13, 14, 16
set buoy=16
:: preprocessing variable tells whether to do preprocessing or not
:: 1 is for preprocessing and any other value is for not doing preprocessing.
set preprocess=1

CD py_script_files
python main_script.py %buoy% %preprocess%
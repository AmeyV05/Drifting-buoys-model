# Ice Buoy_data_analysis
A python code to process and analyse buoy data from the SMOS experiment. The code produces simulation output for free drift ice motion and compares it with buoy drift. The generated data is stored in the folder "generated data"

Open the folder py_scripts and run the file "main_script.py"

Please note that you need the following packages:

1. Python 3.6 
2. Netcdf4  
3. Cartopy (plotting on a world map)

It is assumed that you already have numpy and scipy packages installed.


Note that one should be connected to the Deltares P: project directory for the code to successfully run.
One needs access to the location:P:\1230882-emodnet_hrsm\fromAmey\container\Data_from_models 

If you don't have access to the directory please contact the author for obtaining the required data. 

It is highly recommended to copy this folder and run the code if you are connected to Deltares project folder through a remote connection like NetExtender.

Sometimes there is a runtime error while reading the ncfiles from above location.
In such a situation, please copy the folder 'Data_from_models' to the project working directory and run the script again.



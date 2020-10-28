# Model for drifting ice buoys in Barents Sea region.
This project aims to process and analyse buoy data obtained from an experiment to verify a ship route optimization algorithm. The theoretical aspects of the project are published in the EGU2020 proceedings. Following is the link:

https://doi.org/10.5194/egusphere-egu2020-7544 


The python based code produces simulation output for free drift ice motion and compares it with buoy drift. The generated data is stored in the folder "generated_data"


Note that one should be connected to the Deltares P: project directory for the code to successfully run.
One needs access to the folders: 
(Windows hence ' $backslash$ ' and not '/')
P:\1230882-emodnet_hrsm\fromAmey\container\Data_from_models 
P:\1230882-emodnet_hrsm\fromAmey\container\Buoy_data_excel

In case you are connected to Deltares via a VPN then, it is highly recommended to copy these folders in to a directory called 'container' outside the repository "Drifting-ice-buoys" and run the code. 

If you don't have access to the above folders please contact the author for obtaining the required data. 


**Please use release 1.2 onwards. There was an error in velocity calucations which was later corrected**

***Updated version by the name of Nova1.0 is now released***

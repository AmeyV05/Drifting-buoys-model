# main script
## This is the main script you call to run the code. It will do the pre-processing scripts and model runs together. 
## The script gives you an option of pre-processing.
import numpy as np
import runningmodel 
import processingdata 
import readingdata 
import logging
import generalfunc as gf
import datetime
import os
import sys

# system arguments for running main script.
Bnum = sys.argv[1]
Prepro =sys.argv[2]
# 
### creating of log file for the run. 
# Create a folder for generated_data
gendataloc='../../generated_data'
gf.mkdir_p(gendataloc)
#create log file
gf.logcreate(gendataloc)


###Simulation start date.
#We need to provide a start date for the simulation. The  buoy data are available from 15th March, 2014.
# To limit the amount of data generated and increase the speed of the code we pre-process the data also from the given start date.
sdate=96*16   #96 corresponds to 1 day. As the data is in interval of 15mins.
osdte = datetime.datetime(2014, 3, 15)
dt=datetime.timedelta(minutes=int(sdate*15))
sdte=osdte+dt
#16 days start offset is based on the fact that free drift is only available from below 78deg latitude.
#all the buoys reach this latitude or below by this time. 
logging.info("Simulation started for " +str(sdte))

### Selection of buoy
def buoyselect(Bnum,sdate):
	### buoy indexing as each buoy data is available for different end dates.
	# buoy number decides the indexing '03' : 5952
	logging.info("Running for: Buoy_"+Bnum)
	switcher={
		'02':3152,
		'03':5952,
		'09':3160,
		# '07':1192, as data only till 28 MArch
		# '12':3024, buoy stuck in the area. Doesn't move ahead.
		'13':2524,
		'14':2136,
		'16':5952
			}
	indexing=switcher.get(Bnum,"Invalid Buoy number")-sdate
	return(indexing)


#### Preprocessing
# Pre-processing is important for running the code for the first time. 
def preprocess(Bnum,sdate,indexing):
    #Define nc file locations
    if os.path.isdir('../../container/Model_n_buoy_data'):
      fileloc="../../container/Model_n_buoy_data"
    else:
      fileloc="P:/1230882-emodnet_hrsm/fromAmey/container/Model_n_buoy_data"
    FD=readingdata.readingdata(fileloc,Bnum)
    processingdata.processingdata(Bnum,indexing,FD,sdate)
    logging.info("Pre-processing completed.")

def main():
	#This is just for trial
	indexing=buoyselect(Bnum,sdate)
	if Prepro=='1':
		logging.info("Pre-processing started.")
		preprocess(Bnum,sdate,indexing)
	else:
		pass
	runningmodel.run(Bnum,indexing)

if __name__ == '__main__':
	main()


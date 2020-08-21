# main script
## This is the main script you call to run the code. It will do the pre-processing scripts and model runs together. 
## The script gives you an option of pre-processing.
import numpy as np
import icesimulate as isim
import processingdata 
import readingdata 
import logging
import generalfunc as gf
import datetime
import os
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
# all the buoys read this latitude or below by this time. 
logging.info("Simulation started for " +str(sdte))

### Important change dated 20/8/2020
#This was the original idea to create a filtered buoy observations and store in the pos_vel_data excel file.
# now this data set is no longer required and hence this part of code is removed
#numtaps corresponds to filtered buoy observations.
# # numtaps=2*24*2+1
# numtaps=0
# fedge=int(numtaps/2)


### Selection of buoy
# you can either run the code for a particular Buoy or all the buoys together.

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
	Bnum='02'
	indexing=buoyselect(Bnum,sdate)
	preprocess(Bnum,sdate,indexing)

if __name__ == '__main__':
	main()



### Model 
# This part is independent of the pre-processing. 
# But one needs to have the pre-processed data in the form of Pos_VelData.xlsx file to run the model.

# This describes which forces and parameters are included. 
#forcevec=[f,h,Ua, Va, Ut, Vt, Uo, Vo,Pgx,Pgy,Pgxt,Pgyt]  #

# itervec=[0,1,0.8,0.5,0.2,0.1,0,0,0,0,0,0,0,0,0]
# itervec=[1]
# param=[0.5]
# itervec=[1]
# def model(Bnum,indexing,forcevec):
# 	isim.body(Bnum,indexing,numtaps,Cor)	




# # #deciding whether to run the code for all the buoys or 1 buoy. 
# logging.info("Do you wish to run the code for all the buoys or a particular buoy? ")
# bcount=int(input("Type 1, if you want to run for all the buoys else any other number for a particular buoy data: "))

# #checking if all the nc files should be read to create a pos vel excel file?
# logging.info("Do you wish to read and process data from all the NC files? ")
# logging.info("Please note that this is a computationally heavy and time demanding operation.")
# logging.info("So, if you have already created the appropriate position & velocity xlsx file don't perform this operation. ")
# logging.info("Generally, first time users only do this operation.") 
# ncdatcount=int(input("Type 1 if you want to read all the NC files and create .xlsx files. Else, any other number: "))

# for i in range(len(itervec)):
# 	Cor=[1,1,1,1,1,1,1,1,1,1,1,1]
# 	# Cor=[1,1,1,1,0,0,0,0,0,0,1,1]   #for pgt, coriolis, air 
# 	# Cor=[1,1,1,1,1,1,1,1,1,1,0,0] 
# 	if i==0:
# 		Cor[i]=itervec[i]
# 		# Cor[i+1]=itervec[i]
# 		print(i)
# 		print(Cor[i])
# 	elif 0<i<=5: 
# 		Cor[1]=itervec[i]
# 		print(i)
# 		print(Cor[i])
# 	elif i>5:
# 		Cor[i-4]=itervec[i]
# 		Cor[i-3]=itervec[i]
# 		print(i)
# 		print(Cor[i-4])
# 		i=i+1
# 	# Cor[11]=param[i]
# 	if (bcount!=1):
# 	  #Getting buoy number input.
# 	  Bnum=input("Please input the Buoy number from the list: [02, 03, 07, 09, 12, 13, 14, 16]. And press enter: ")
# 	  logging.info("Running for: Buoy_"+Bnum)
	  
# 	  if (ncdatcount==1):
# 	    logging.info("Running the script for reading all the .NC data files.")
# 	    #Define nc file locations
# 	    if os.path.isdir('../../container/Data_from_models'):
# 	      fileloc="../../container/Data_from_models"
# 	    else:
# 	      fileloc="P:/1230882-emodnet_hrsm/fromAmey/container/Data_from_models"
	    

# 	    FD=rdnc.readingalldata(fileloc,Bnum)
# 	    ncdat.mainproc(Bnum,indexing,FD,numtaps,sdate)
# 	  isim.body(Bnum,indexing,numtaps,Cor)
# 	  # emisim.body(Bnum,indexing,numtaps,Cor)
# 	else:
# 	  for i in switcher:
# 	  	# if i=='02' or i=='03' or i=='09':
# 	  	# 	continue
# 	  	# else:
# 	    Bnum=i
# 	    logging.info("Running for: Buoy_"+Bnum)
# 	    indexing=switcher.get(Bnum,"Invalid Buoy number")-sdate-fedge
# 	    if (ncdatcount==1):
# 	      logging.info("Running the script for reading all the .NC data files.")
# 	      #Define nc file locations
# 	      if os.path.isdir('../../container/Data_from_models'):
# 	        fileloc="../../container/Data_from_models"
# 	      else:
# 	        fileloc="P:/1230882-emodnet_hrsm/fromAmey/container/Data_from_models"
# 	      FD=rdnc.readingalldata(fileloc,Bnum)
# 	      ncdat.mainproc(Bnum,indexing,FD,numtaps,sdate)
# 	    isim.body(Bnum,indexing,numtaps,Cor)
# 		    # emisim.body(Bnum,indexing,numtaps,Cor)



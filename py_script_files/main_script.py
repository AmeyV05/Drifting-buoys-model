# main script file.
import numpy as np
import icesimulate as isim
import enkfmain as enkf
import processdata as ncdat
import readallNC as rdnc
import os
import logging
import sys
import optCw as oCw
import iteratescript as itscrpt 
import generalfunc as gf
# buoy number decides the indexing
# '03' : 5952

for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

gendata='../generated_data'
gf.mkdir_p(gendata)
# logging changes to create a log file and also show log in the terminal.
logging.basicConfig(filename=gendata+'/out.log', filemode='w', level=logging.DEBUG)
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# add the handler to the root logger
logging.getLogger('').addHandler(console)


#indexing for different dates and times based on available buoy data.
switcher={
  '02':3152,
  '03':5952,
  '09':3160,
  # '07':1192, as data only till 28 MArch
  '12':3024,
  '13':2524,
  '14':2136,
  '16':5952
  }

#numtaps corresponds to filtered buoy observations.
numtaps=2*24*2+1
#numtaps=0
fedge=int(numtaps/2)
sdate=96*16
#96 means 1 day

# #deciding whether to run the code for all the buoys or 1 buoy. 
logging.info("Do you wish to run the code for all the buoys or a particular buoy? ")
bcount=int(input("Type 1, if you want to run for all the buoys else any other number for a particular buoy data: "))

#checking if all the nc files should be read to create a pos vel excel file?
logging.info("Do you wish to read and process data from all the NC files? ")
logging.info("Please note that this is computationally heavy and time demanding operation.")
logging.info("So, if you have already created the appropriate position & velocity xlsx file don't perform this operation. ")
logging.info("Generally, first time users only do this operation.") 
ncdatcount=int(input("Type 1 if you want to read all the NC files and create .xlsx files. Else, any other number: "))

# This describes which forces and parameters are included. 
#cor=[f,h,Ua, Va, Ut, Vt, Uo, Vo,Pgx,Pgy,Pgxt,Pgyt]
Cor=[1,1,1,1,1,1,1,1,1,1,1,1]

if (bcount!=1):
  #Getting buoy number input.
  Bnum=input("Please input the Buoy number from the list: [02, 03, 07, 09, 12, 13, 14, 16]. And press enter: ")
  logging.info("Running for: Buoy_"+Bnum)
  indexing=switcher.get(Bnum,"Invalid Buoy number")-sdate-fedge
  if (ncdatcount==1):
    logging.info("Running the script for reading all the .NC data files.")
    #Define nc file locations
    if os.path.isdir('../../container/Data_from_models'):
      fileloc="../../container/Data_from_models"
    else:
      fileloc="P:/1230882-emodnet_hrsm/fromAmey/container/Data_from_models"
    FD=rdnc.readingalldata(fileloc)
    ncdat.mainproc(Bnum,indexing,FD,numtaps,sdate)
  isim.body(Bnum,indexing,numtaps,Cor)

else:
  for i in switcher:
    Bnum=i
    logging.info("Running for: Buoy_"+Bnum)
    indexing=switcher.get(Bnum,"Invalid Buoy number")-sdate-fedge
    if (ncdatcount==1):
      logging.info("Running the script for reading all the .NC data files.")
      #Define nc file locations
      if os.path.isdir('../../container/Data_from_models'):
        fileloc="../../container/Data_from_models"
      else:
        fileloc="P:/1230882-emodnet_hrsm/fromAmey/container/Data_from_models"
      FD=rdnc.readingalldata(fileloc)
      ncdat.mainproc(Bnum,indexing,FD,numtaps,sdate)
    isim.body(Bnum,indexing,numtaps,Cor)




## To run ensemble filter scripts and optimization
# itscrpt.main(Bnum,indexing,Cor)
# enkf.body(Bnum,indexing,Cor)
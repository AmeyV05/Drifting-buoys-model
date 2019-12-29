# main script file.

import numpy as np
import body_program as BP
import readallNC as rdnc
import os
# buoy number decides the indexing
# '03' : 5856

def bindexing(Bnum,switcher):
 return switcher.get(Bnum,"Invalide Buoy number")

#

#nc file locations

#nc file locations
if os.path.isdir('../Data_from_models'):
 fileloc="../Data_from_models"
else:
 fileloc="P:/1230882-emodnet_hrsm/fromAmey/container/Data_from_models"
def readingalldata(fileloc):
 #####GTSM data
 print("Reading GTSM data")
 file=fileloc+"/gtsm_truncated.nc"
 [Xt,Yt,ut,vt,Tt]=rdnc.read_GTSM_map(file)
 TD={'Xt': Xt,'Yt': Yt,'ut': ut,'vt': vt,'Tt': Tt}
 #### ERA 5 winds
 print("Reading ERA5 wind data.")
 #ERA 5 winds
 file=fileloc+"/era5_wind_201403_05.nc"
 [Xa,Ya,u10,v10,Ta]=rdnc.read_wind(file)
 AD={'Xa': Xa,'Ya': Ya,'u10': u10,'v10': v10,'Ta': Ta}
 ## Ocean currents
 print("Reading ocean currents data.")
 file=fileloc+"/Ocean_currents_buoy.nc"
 [Xo,Yo,uo,vo,To]=rdnc.read_ocean(file)
 OD={'Xo': Xo,'Yo': Yo,'uo': uo,'vo': vo,'To': To}
 return(TD,AD,OD)
switcher={
  '02':3152,
  '03':5952,
  '09':3160,
  '07':1192,
  '12':3024,
  '13':2524,
  '14':2136,
  '16':5952
  }
#Bnum="15"
print("Do you wish to run the code for all the buoys or a particular buoy? ")
print("Type 1, if you want to run for all the buoys else any other number for a particular buoy data.")
count=int(input())
Cor=int(input("Do you want Coriolis force? If yes, type 1 else any other number: "))
[TD,AD,OD]=readingalldata(fileloc)
if (count==1):
 for i in switcher:
  Bnum=i
  indexing=switcher[i]-96
  BP.main(Bnum,indexing,Cor,TD,AD,OD)
else:
 print("Please input the Buoy number from the list: [02, 03, 07, 09, 12, 13, 14, 16]. And press enter.")
 Bnum=input()
 indexing=switcher.get(Bnum,"Invalid Buoy number")-96
 BP.main(Bnum,indexing,Cor,TD,AD,OD)



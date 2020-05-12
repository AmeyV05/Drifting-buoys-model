#buoy_data_processing
#Using pandas for excel and csv 
# Please note that velocity vector is 1 less than the time vectors as it is a forward difference method.
import pandas as pd 
import numpy as np
import generalfunc as gf
np.set_printoptions(threshold=np.inf)


def bdata(file):
 D=pd.read_csv(file)
 Xib=np.array(D['Lon'])
 Yib=np.array(D['Lat'])
 Tib=np.array(D['Date(GMT)'])
 D.head()
 #velocity computation using the buoy data
 dt=15*60
 Uib=[]
 Vib=[]
 for j in range(len(Xib)-1):
  dlon=Xib[j+1]-Xib[j]
  dlat=Yib[j+1]-Yib[j]
  lat=Yib[j]
  (dlatm,dlonm)=gf.latlon2meters(lat,dlat,dlon)
  Uib=np.append(Uib,(dlonm)/dt)
  Vib=np.append(Vib,(dlatm)/dt)
 return(Xib,Yib,Uib,Vib,Tib)

def main():
 file ='../buoy_data_excel/procesd_buoy_loc_data_Lars/BUOY_02.csv'
 [Xib,Yib,Uib,Vib,Tib]=bdata(file)
 print(Tib)

if __name__ == '__main__':
	main()


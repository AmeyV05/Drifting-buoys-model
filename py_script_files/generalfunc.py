#general useful functions 
# change in lat long to meters
# change in meters to  change in lat long
# Fourier Transforms

import numpy as np
import datetime
import pandas as pd
import math
deg2rad= np.pi/180.0
# A conversion function to convert the latitudes and longitudes to meteres.
# note that this is a simple conversion formula and not very accurate
def latlon2meters(lat,dlat,dlon):
 #lat and lon are changes in lat and lon in deg
 dlatm=dlat*111.32e03
 dlonm=dlon*40075e03*np.cos(deg2rad*lat)/360
 return(dlatm,dlonm)

#conversion of meters to latitude longitude and coriolis calculator
# conversion of meters to latitude longitude.

def m2latlon(lat,dlonm,dlatm):
 dlat=dlatm/(111.32e03)
 dlon=dlonm/(111.32e03*np.cos(deg2rad*lat))
 return(dlon,dlat)


def num2datetimehrs(y,m,d,sindex,eindex,tvec):
 x=datetime.datetime(y, m, d)
 y=[] #tvec modified to date time.
 #index=336  #indexing for going from 1st march to 15th march (only for winds)
 for i in range(sindex, eindex):
  h=tvec[i]
  dt=datetime.timedelta(hours=int(h))
  y1=dt+x
  y=np.append(y,y1.strftime("%d/%m/%Y %H:%M:%S"))
 #print(y)
 #print(len(time))
 return(y)

def num2datetimesecs(y,m,d,sindex,eindex,tvec):
 x=datetime.datetime(y, m, d)
 y=[] #tvec modified to date time.
 #index=336  #indexing for going from 1st march to 15th march (only for winds)
 for i in range(sindex, eindex):
  s=tvec[i]
  dt=datetime.timedelta(seconds=int(s))
  y1=dt+x
  y=np.append(y,y1.strftime("%d/%m/%Y %H:%M:%S"))
 #print(y)
 #print(len(time))
 return(y)

def mkdir_p(mypath):
#Creates a directory. equivalent to using mkdir -p on the command line'''
 from errno import EEXIST
 from os import makedirs,path
 try:
  makedirs(mypath)
  print("Buoy directory created")
 except OSError as exc: # Python >2.5

  if exc.errno == EEXIST and path.isdir(mypath):
   print("Buoy directory already exists.")
   pass
  else: raise

def veltransform(Uib,Vib,Ut,Vt,Ua,Va,Uo,Vo,path):
 Utvec=np.column_stack((Ut,Vt))
 #Utvec=Utvec[:-1] #don't consider the last element as vel vector is 1 less than time
 Uibvec=np.column_stack((Uib,Vib))
 Uavec=np.column_stack((Ua,Va))
 Uavec=np.repeat(Uavec,4,axis=0)
 Uovec=np.column_stack((Uo,Vo))
 Uovec=np.repeat(Uovec,4,axis=0)
 Uwvec=Utvec+Uovec
 return(Utvec,Uibvec,Uavec,Uovec,Uwvec)

def save2excel(Tib,Utvec,Uibvec,Uisvec,Uovec,Uwvec,Uavec,path):
 # Create some Pandas dataframes from some data.
 dfd = pd.DataFrame({'Date(GMT)': Tib})
 dfibu = pd.DataFrame({'U': Uibvec[:,0]})
 dfibv = pd.DataFrame({'V': Uibvec[:,1]})
 dfisu = pd.DataFrame({'U': Uisvec[:,0]})
 dfisv = pd.DataFrame({'V': Uisvec[:,1]})
 dftu = pd.DataFrame({'U': Utvec[:,0]})
 dftv = pd.DataFrame({'V': Utvec[:,1]})
 dfau = pd.DataFrame({'U': Uavec[:,0]})
 dfav = pd.DataFrame({'V': Uavec[:,1]})
 dfou = pd.DataFrame({'U': Uovec[:,0]})
 dfov = pd.DataFrame({'V': Uovec[:,1]})
 dfwu = pd.DataFrame({'U': Uwvec[:,0]})
 dfwv = pd.DataFrame({'V': Uwvec[:,1]})
 #df4 = pd.DataFrame({'Data': [41, 42, 43, 44]})

 # Create a Pandas Excel writer using XlsxWriter as the engine.
 writer = pd.ExcelWriter(path+'/Vel_data.xlsx', engine='xlsxwriter')

 # Position the dataframes in the worksheet.
 dfd.to_excel(writer, startcol=0,startrow=1,index=False)  # Default position, cell A1.
 #buoy
 dfibu.to_excel(writer, startcol=1,startrow=1,index=False)
 dfibv.to_excel(writer, startcol=2,startrow=1,index=False)
 #simulated
 dfisu.to_excel(writer, startcol=3,startrow=1,index=False)
 dfisv.to_excel(writer, startcol=4,startrow=1,index=False)
 #GTSM
 dftu.to_excel(writer, startcol=5,startrow=1,index=False)
 dftv.to_excel(writer, startcol=6,startrow=1,index=False)
 #Ocean
 dfou.to_excel(writer, startcol=7,startrow=1,index=False)
 dfov.to_excel(writer, startcol=8,startrow=1,index=False)
 #Ocean+Tidal
 dfwu.to_excel(writer, startcol=9,startrow=1,index=False)
 dfwv.to_excel(writer, startcol=10,startrow=1,index=False)
 #ERA 5 winds
 dfau.to_excel(writer, startcol=11,startrow=1,index=False)
 dfav.to_excel(writer, startcol=12,startrow=1,index=False)
 # Close the Pandas Excel writer and output the Excel file.
 #writer.save()
 workbook  = writer.book
 worksheet = writer.sheets['Sheet1']
 merge_format = workbook.add_format({
    'bold': 1,
    'border': 1,
    'align': 'center',
    'valign': 'vcenter',
     })
 worksheet.merge_range("B1:C1","Buoy", merge_format)
 worksheet.merge_range("D1:E1","Simulated", merge_format)
 worksheet.merge_range("F1:G1","GTSM", merge_format)
 worksheet.merge_range("H1:I1","Ocean currents", merge_format)
 worksheet.merge_range("J1:K1","Ocean+Tides", merge_format)
 worksheet.merge_range("L1:M1","ERA5", merge_format)
 workbook.close()
def index2time(indexing):
 tothours=indexing/4	
 ndays=math.floor(tothours/24)
 nhrs=tothours-(ndays*24)
 return(ndays,nhrs)

def  main():
 Tib=[2,3,4]
 Utvec=np.zeros((4,2))
 Uibvec=np.zeros((4,2))
 Uisvec=np.zeros((4,2))
 Uovec=np.zeros((4,2))
 Uwvec=np.zeros((4,2))
 Uavec=np.zeros((4,2))
 path="C:/Users/vasulkar/OneDrive - Stichting Deltares/Documents/Research_data/buoy_analysis/py_script_files"
 save2excel(Tib,Utvec,Uibvec,Uisvec,Uovec,Uwvec,Uavec,path)

if __name__ == '__main__':
	main()
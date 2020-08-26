#general useful functions 
# change in lat long to meters
# change in meters to  change in lat long
# Fourier Transforms
import scipy.ndimage as ndimage
import scipy.signal as signal
import numpy as np
import datetime
import pandas as pd
import math
import shutil
import os
import settings
import logging
deg2rad= np.pi/180.0

## function to start creating a log file named out.log

def logcreate(gendataloc):
  #remove any redundancies with previous log handlers.
  for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
  # logging changes to create a log file and also show log in the terminal.
  logging.basicConfig(filename=gendataloc+'/out.log', filemode='w', level=logging.DEBUG)
  # define a Handler which writes INFO messages or higher to the sys.stderr
  console = logging.StreamHandler()
  console.setLevel(logging.INFO)
  # add the handler to the root logger
  logging.getLogger('').addHandler(console)
  logging.disable(logging.DEBUG)

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

def rotmatrix(phi):
  #this function creates a rotation matrix phi is in rad.
  R=np.array([[np.cos(phi),np.sin(phi)],[-np.sin(phi),np.cos(phi)]])
  return(R)

def num2datetimesecs(y,m,d,sindex,eindex,tvec):
 x=datetime.datetime(y, m, d)
 xtot=[] #tvec modified to date time.
 #index=336  #indexing for going from 1st march to 15th march (only for winds)
 for i in range(sindex, eindex+1):
  s=tvec[i]
  dt=datetime.timedelta(seconds=int(s))
  y1=dt+x
  xtot=np.append(xtot,y1.strftime("%d/%m/%Y %H:%M:%S"))
 #print(y)
 #print(len(time))
 return(xtot)

def mkdir_p(mypath):
#Creates a directory. equivalent to using mkdir -p on the command line'''
 from errno import EEXIST
 from os import makedirs,path
 basename=path.basename(mypath)
 try:
  makedirs(mypath)

  print(basename+" directory created")
 except OSError as exc: # Python >2.5

  if exc.errno == EEXIST and path.isdir(mypath):
   print(basename+" directory already exists.")
   pass
  else: raise

def index2time(indexing):
  #This function gives you the number of days and hours based on the indexing value.
  # it is suited to give values based on 15 mins time step. 
  tothours=indexing/4	
  ndays=math.floor(tothours/24)
  nhrs=tothours-(ndays*24)
  return(ndays,nhrs)


def icethicktyp(forcevec,trate):
  s=settings.settings()
  forcevecn=[]
  if forcevec[1]!='v':
    forcevecn=forcevec
    trate=0
    logging.info("Running with constant ice thickness.")
  else:
    forcevecn=forcevec[:]
    forcevecn[1]=1
  return(forcevecn,trate)
  
def forcedetail(Cor,trate,hs):
  #this function gives an output string array
  #which describes the forces and parameters of the simulations.
  f='Yes' if Cor[0]==1 else 'No' #Coriolis 
  h=str(hs)+str(Cor[1])+"_"+str(trate) if Cor[1]=='v' else str(Cor[1]) #ice thickness
  Uax='Yes' if Cor[2]==1 else 'No' #wind stresss
  Uay=Uax
  Utx='Yes' if Cor[4]==1 else 'No' #Tidal stress
  Uty=Utx
  Uox='Yes' if Cor[6]==1 else 'No' #ocean
  Uoy=Uox
  Pgxo='Yes' if Cor[8]==1 else 'No' #pressure gradients
  Pgyo=Pgxo
  Pgxt='Yes' if Cor[10]==1 else 'No' #pressure gradients tides
  Pgyt=Pgxt
  Cornam=[f,h,Uax,Uay,Utx,Uty,Uox,Uoy,Pgxo,Pgyo,Pgxt,Pgyt]
  #Folder name for storing the data.
  folname='h'+str(Cor[1])+'f'+str(Cor[0])+ \
          'A'+str(Cor[2])+'T'+str(Cor[4])+ \
          'O'+str(Cor[6])+'Po'+str(Cor[8])+'Pt'+str(Cor[10]) if Cor[1]!='v' else \
          'h'+str(hs)+str(Cor[1])+str(trate)[1:4]+'f'+str(Cor[0])+ \
          'A'+str(Cor[2])+'T'+str(Cor[4])+ \
          'O'+str(Cor[6])+'Po'+str(Cor[8])+'Pt'+str(Cor[10])
  return(Cornam,folname)

#function to copy the log file to the folder of simulations.
def logcopy(path):
  srfile='../../generated_data/out.log'
  shutil.copy(srfile,path+'/out.log')
  print("Log file available in: "+path)

## low pass filter
def LPfilter (numtaps,px):     
  window=signal.hamming(numtaps)
  window=window/sum(window)
  pxfilter=signal.convolve(px, window, mode='valid') 
  edge=int(numtaps/2)
  pxres=px[edge:-edge]-pxfilter
  px=px[edge:-edge] #shortening px for later use.
  return (pxres,px,pxfilter)
#subroutine for FFT 
def FFTsig(y,N,ts):
  Yamp=abs(np.fft.fft(y))/N
  Ypha=np.angle(np.fft.fft(y),deg=True)
  if (N%2==1):
      N1=int((N-1)/2)
      Yampplt=Yamp[:N1]
      Yphaplt=Ypha[:N1]
  else:
      N1=int((N/2))
      Yampplt=Yamp[:N1]
      Yphaplt=Ypha[:N1]
  #sampling frequency
  fs=1/ts #ts is sampling time in seconds
  f_nq=fs/2 #Nyquist frequency (Fmax)
  fvec=np.linspace(0,f_nq,N1)
  tvec=1/fvec[1:]
  tvec=np.append(np.inf,tvec)
  return(Yampplt,Yphaplt,fvec,tvec)

def thinrate(ho,trate):
  dtobs=15*60
  hn=ho+trate*dtobs
  if hn<0.1:
    logging.info("Error: Ice thickness decreasing to less than 0.1 ")
    logging.info("Using constant ice thickness of 0.1m")
    hn=0.1
  return(hn)


def  main(): # doesn't work anymore 
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
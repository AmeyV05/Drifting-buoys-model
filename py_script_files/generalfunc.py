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
 tothours=indexing/4	
 ndays=math.floor(tothours/24)
 nhrs=tothours-(ndays*24)
 return(ndays,nhrs)

#subroutine for FFT 
def FFT_signal(y,N,ts):
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


def TidCorio_comput(tvec,tide,Cordeg): 
    deg1=Cordeg['deg1'];deg2=Cordeg['deg2'] 
    #coriolis frequ
    phi1=np.deg2rad(deg1);phi2=np.deg2rad(deg2)
    coriolis=(2*2*np.pi)/(86164.09/3600)
    corperd1=2*np.pi/(coriolis*np.sin(phi1))
    corperd2=2*np.pi/(coriolis*np.sin(phi2))
    arg_deg1=(np.abs(tvec-corperd1)).argmin()
    arg_deg2=(np.abs(tvec-corperd2)).argmin()
    arg_cor={'argdeg1':arg_deg1,'argdeg2':arg_deg2}
    arg_m2 = (np.abs(tvec-tide['M2'])).argmin()
    arg_s2 = (np.abs(tvec-tide['S2'])).argmin()
    arg_mu2 = (np.abs(tvec-tide['MU2'])).argmin()
    arg_o1 = (np.abs(tvec-tide['O1'])).argmin()
    arg_k1 = (np.abs(tvec-tide['K1'])).argmin()
    arg_m4 = (np.abs(tvec-tide['M4'])).argmin()
    arg_tide={'M2':arg_m2,'S2':arg_s2,'MU2':arg_mu2,'O1':arg_o1,'K1':arg_k1,'M4':arg_m4}
    return(arg_tide,arg_cor)

def LPfilter (numtaps,p_x):     
    window=signal.hamming(numtaps)
    window=window/sum(window)
    p_xfilter=signal.convolve(p_x, window, mode='valid') 
    edge=int(numtaps/2)
    p_xresiduum=p_x[edge:-edge]-p_xfilter
    p_x=p_x[edge:-edge]
    return (p_xresiduum,p_x,p_xfilter)


#Error statistics function 

def errstats(Xobs,Yobs,Xsim,Ysim,tmplierinv):
 Xsim=Xsim[::tmplierinv]
 Ysim=Ysim[::tmplierinv]
 dx=(Xobs - Xsim);dy=(Yobs - Ysim);da=np.sqrt((dx**2+dy**2))
 #mean error
 merrx=dx.mean();merry=dy.mean();merra=da.mean()
 merr=[merrx,merry,merra]
 #rms error
 rmsx=np.sqrt((dx ** 2).mean())
 rmsy=np.sqrt((dy ** 2).mean())
 rmsa=np.sqrt((da ** 2).mean())
 rms=[rmsx,rmsy,rmsa]
 #weighted error
 N=len(dx)
 s=0;
 for i in range(N): s+=(N-i)*dx[i]
 werrx=(1/(2*N*(N+1)))*s
 s=0; 
 for i in range(N): s+=(N-i)*dy[i]
 werry=(1/(2*N*(N+1)))*s
 for i in range(N): s+=(N-i)*da[i]
 werra=(1/(2*N*(N+1)))*s
 werr=[werrx,werry,werra]
 return(merr,rms,werr)

def Cordesfunc(Cor,trate,hs):
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

#getting FT by subracting the mean drift and obtaining tidal components for M2 and coriolis 
def FTremMD(numtaps,Xib,Xis,tmplierinv):
  Xis=Xis[::tmplierinv]
  #filtering with lowpas filer
  [Xbfil,X1,xfilter]=LPfilter (numtaps,Xib)
  [Xsfil,X1,xfilter]=LPfilter (numtaps,Xis)
  Nft=len(Xbfil);dt=15*60.0 #time difference in observations
  [xbamres,xbphres,fvec,tvec]=FFT_signal(Xbfil,Nft,dt) 
  [xsamres,xsphres,fvec,tvec]=FFT_signal(Xsfil,Nft,dt)
  tvec=tvec/3600 #Period in hours
  #computation of tidal and coriolis frequency arguments
  tide={'M2':12.421,'S2':12.,'MU2':12.871,'O1':25.819,'K1':23.934,'M4':6.2103}
  Cordeg={'deg1':74.7,'deg2':79}  #latitude for coriolis
  [arg_tide,arg_cor]= TidCorio_comput(tvec,tide,Cordeg)
  errftam=errorFT(xbamres,xsamres,arg_tide)
  errftph=errorFT(xbphres,xsphres,arg_tide)
  xbres=np.row_stack((xbamres,xbphres))
  xsres=np.row_stack((xsamres,xsphres))
  tidamb=[xbamres[arg_tide['M2']],xbamres[arg_tide['S2']],
          xbamres[arg_tide['MU2']],xbamres[arg_tide['O1']],
          xbamres[arg_tide['K1']],xbamres[arg_tide['M4']]]
  tidams=[xsamres[arg_tide['M2']],xsamres[arg_tide['S2']],
          xsamres[arg_tide['MU2']],xsamres[arg_tide['O1']],
          xsamres[arg_tide['K1']],xsamres[arg_tide['M4']]]  
  tidphb=[xbphres[arg_tide['M2']],xbphres[arg_tide['S2']],
          xbphres[arg_tide['MU2']],xbphres[arg_tide['O1']],
          xbphres[arg_tide['K1']],xbphres[arg_tide['M4']]]
  tidphs=[xsphres[arg_tide['M2']],xsphres[arg_tide['S2']],
          xsphres[arg_tide['MU2']],xsphres[arg_tide['O1']],
          xsphres[arg_tide['K1']],xsphres[arg_tide['M4']]]  
  tidb=np.row_stack((tidamb,tidphb))
  tids=np.row_stack((tidams,tidphs))
  return(tvec,xbres,xsres,Xbfil,Xsfil,arg_tide,arg_cor,errftam,errftph,tidb,tids)

def errorFT(xbres,xsres,arg_tide):
  errM2=xbres[arg_tide['M2']]-xsres[arg_tide['M2']]
  errS2=xbres[arg_tide['S2']]-xsres[arg_tide['S2']]
  errMU2=xbres[arg_tide['MU2']]-xsres[arg_tide['MU2']]
  errO1=xbres[arg_tide['O1']]-xsres[arg_tide['O1']]
  errK1=xbres[arg_tide['K1']]-xsres[arg_tide['K1']]
  errM4=xbres[arg_tide['M4']]-xsres[arg_tide['M4']]
  errftvec=[errM2,errS2,errMU2,errO1,errK1,errM4]
  return(errftvec)

def thinrate(ho,trate):
  dtobs=15*60
  hn=ho+trate*dtobs
  if hn<0.05:
    logging.info("Error: Ice thickness decreasing to less than 0.1 ")
    logging.info("Using constant ice thickness of 0.1m")
    hn=0.05
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
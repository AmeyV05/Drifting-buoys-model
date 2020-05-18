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
 xtot=[] #tvec modified to date time.
 #index=336  #indexing for going from 1st march to 15th march (only for winds)
 for i in range(sindex, eindex):
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

def veltransform(Uib,Vib,Ut,Vt,Ua,Va,Uo,Vo,Pgx,Pgy,Pgxt,Pgyt,Hi):
 Utvec=np.column_stack((Ut,Vt))
 Pgtvec=np.column_stack((Pgxt,Pgyt))
 #Utvec=Utvec[:-1] #don't consider the last element as vel vector is 1 less than time
 Uibvec=np.column_stack((Uib,Vib))
 Uavec=np.column_stack((Ua,Va))
 Uavec=np.repeat(Uavec,4,axis=0)
 Uovec=np.column_stack((Uo,Vo))
 Uovec=np.repeat(Uovec,4,axis=0)
 Uwvec=Utvec+Uovec
 Pgvec=np.column_stack((Pgx,Pgy))
 Pgvec=np.repeat(Pgvec,4,axis=0)
 Hivec=np.repeat(Hi,4*24,axis=0)
 add=np.repeat(Hivec[-1],4*12,axis=0)
 Hivec=np.concatenate((Hivec,add),axis=0)
 return(Utvec,Uibvec,Uavec,Uovec,Uwvec,Pgvec,Pgtvec,Hivec)

def save2excel(Tib,Utvec,Uibvec,Uovec,Uwvec,Pgvec,Pgtvec,Hivec,Uavec,Xib,Yib,Xibf,Yibf,path):
 # Create some Pandas dataframes from some data.
 dfd = pd.DataFrame({'Date(GMT)': Tib[:-1]})
 dfibu = pd.DataFrame({'Uib': Uibvec[:,0]})
 dfibv = pd.DataFrame({'Vib': Uibvec[:,1]})
 dftu = pd.DataFrame({'Ut': Utvec[:,0]})
 dftv = pd.DataFrame({'Vt': Utvec[:,1]})
 dfau = pd.DataFrame({'Ua': Uavec[:,0]})
 dfav = pd.DataFrame({'Va': Uavec[:,1]})
 dfou = pd.DataFrame({'Uo': Uovec[:,0]})
 dfov = pd.DataFrame({'Vo': Uovec[:,1]})
 dfwu = pd.DataFrame({'Uw': Uwvec[:,0]})
 dfwv = pd.DataFrame({'Vw': Uwvec[:,1]})
 dfpgx = pd.DataFrame({'Pgx': Pgvec[:,0]})
 dfpgy = pd.DataFrame({'Pgy': Pgvec[:,1]})
 dfpgxt = pd.DataFrame({'Pgxt': Pgtvec[:,0]})
 dfpgyt = pd.DataFrame({'Pgyt': Pgtvec[:,1]})
 dfhi=pd.DataFrame({'Hi':Hivec})
 dfibx = pd.DataFrame({'Xib': Xib[:-1]})
 dfiby = pd.DataFrame({'Yib': Yib[:-1]})
 dfibxf = pd.DataFrame({'Xibf': Xibf[:-1]})
 dfibyf = pd.DataFrame({'Yibf': Yibf[:-1]})

 # Create a Pandas Excel writer using XlsxWriter as the engine.
 writer = pd.ExcelWriter(path+'/Pos_Vel_data.xlsx', engine='xlsxwriter')

 # Position the dataframes in the worksheet.
 dfd.to_excel(writer, startcol=0,startrow=0,index=False)  # Default position, cell A1.
 #buoy
 dfibu.to_excel(writer, startcol=1,startrow=0,index=False)
 dfibv.to_excel(writer, startcol=2,startrow=0,index=False)
 #GTSM
 dftu.to_excel(writer, startcol=3,startrow=0,index=False)
 dftv.to_excel(writer, startcol=4,startrow=0,index=False)
 #Ocean
 dfou.to_excel(writer, startcol=5,startrow=0,index=False)
 dfov.to_excel(writer, startcol=6,startrow=0,index=False)
 #Ocean+Tidal
 dfwu.to_excel(writer, startcol=7,startrow=0,index=False)
 dfwv.to_excel(writer, startcol=8,startrow=0,index=False)
 #ERA 5 winds
 dfau.to_excel(writer, startcol=9,startrow=0,index=False)
 dfav.to_excel(writer, startcol=10,startrow=0,index=False)
 # Pressure gradients
 dfpgx.to_excel(writer, startcol=11,startrow=0,index=False)
 dfpgy.to_excel(writer, startcol=12,startrow=0,index=False)
 # Pressure gradients
 dfpgxt.to_excel(writer, startcol=13,startrow=0,index=False)
 dfpgyt.to_excel(writer, startcol=14,startrow=0,index=False)
 # Ice thickness
 dfhi.to_excel(writer, startcol=15,startrow=0,index=False)
 #buoy locs
 dfibx.to_excel(writer, startcol=16,startrow=0,index=False)
 dfiby.to_excel(writer, startcol=17,startrow=0,index=False)
  #buoy locs filtered
 dfibxf.to_excel(writer, startcol=18,startrow=0,index=False)
 dfibyf.to_excel(writer, startcol=19,startrow=0,index=False)
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
 # worksheet.merge_range("B1:C1","Buoy", merge_format)
 # worksheet.merge_range("D1:E1","Simulated", merge_format)
 # worksheet.merge_range("F1:G1","GTSM", merge_format)
 # worksheet.merge_range("H1:I1","Ocean currents", merge_format)
 # worksheet.merge_range("J1:K1","Ocean+Tides", merge_format)
 # worksheet.merge_range("L1:M1","ERA5", merge_format)
 workbook.close()

def index2time(indexing):
 tothours=indexing/4	
 ndays=math.floor(tothours/24)
 nhrs=tothours-(ndays*24)
 return(ndays,nhrs)

#subroutine for FFT 
def FFT_signal(y,N,ts):
    Y=abs(np.fft.fft(y))/N
    if (N%2==1):
        N1=int((N-1)/2)
        Y_plot=Y[:N1]
    else:
        N1=int((N/2))
        Y_plot=Y[:N1]
    #sampling frequency
    fs=1/ts #ts is sampling time in seconds
    f_nq=fs/2 #Nyquist frequency (Fmax)
    fvec=np.linspace(0,f_nq,N1)
    tvec=1/fvec[1:]
    tvec=np.append(np.inf,tvec)
    return(Y_plot,fvec,tvec)


def TidCorio_comput(tvec,M2,deg):  
    #coriolis frequ
    phi=np.deg2rad(deg)
    coriolis=(2*2*np.pi*np.sin(phi))/(86164.09/3600)
    cori_period=2*np.pi/coriolis
    arg_m2 = (np.abs(tvec-M2)).argmin()
    arg_cor=(np.abs(tvec-cori_period)).argmin()
    return(arg_m2,arg_cor,cori_period)

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

def Cordesfunc(Cor):
  #this function gives an output string array
  #which describes the forces and parameters of the simulations.
  f='Yes' if Cor[0]==1 else 'No' #Coriolis 
  h=str(Cor[1]) #ice thickness
  Uax='Yes' if Cor[2]==1 else 'No' #wind stresss
  Uay=Uax
  Utx='Yes' if Cor[4]==1 else 'No' #Tidal stress
  Uty=Utx
  Uox='Yes' if Cor[6]==1 else 'No' #ocean
  Uoy=Uox
  Pgxo='Yes' if Cor[8]==1 else 'No' #pressure gradients
  Pgyo=Pgxo
  Pgxt='Yes' if Cor[10]==1 else 'No' #pressure gradients
  Pgyt=Pgxt
  Cornam=[f,h,Uax,Uay,Utx,Uty,Uox,Uoy,Pgxo,Pgyo,Pgxt,Pgyt]
  #Folder name for storing the data.
  folname='h'+str(Cor[1])+'f'+str(Cor[0])+ \
          'A'+str(Cor[2])+'T'+str(Cor[4])+ \
          'O'+str(Cor[6])+'P'+str(Cor[8])
  return(Cornam,folname)
# def logcopy(path):
#   srfile='../generated_data/out.log'
#   shutil.copy(srfile,path+'/out.log')
#   print("Log file available in: "+path)

#getting FT by subracting the mean drift and obtaining tidal components for M2 and coriolis 
def FTremMD(numtaps,Xib,Xis,tmplierinv):
  Xis=Xis[::tmplierinv]
  #filtering with lowpas filer
  [Xbres,X1,xfilter]=LPfilter (numtaps,Xib)
  [Xsres,X1,xfilter]=LPfilter (numtaps,Xis)
  Nft=len(Xbres);dt=15*60.0 #time difference in observations
  [xbres,fvec,tvec]=FFT_signal(Xbres,Nft,dt)
  [xsres,fvec,tvec]=FFT_signal(Xsres,Nft,dt)
  tvec=tvec/3600 #Period in hours
  #computation of tidal and coriolis frequency arguments
  M2=12.421 #M2 tidal frequency period
  deg=74.7;deg1=79.  #latitude for coriolis
  [arg_m2,arg_74,cori_period]= TidCorio_comput(tvec,M2,deg)
  [arg_m2,arg_79,cori_period1]= TidCorio_comput(tvec,M2,deg1)
  return(tvec,xbres,xsres,arg_m2,arg_74,arg_79)

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
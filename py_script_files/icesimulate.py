import numpy as np
import generalfunc as gf
import genplots as gp
import logging
import icemodel as im
import pandas as pd
import settings 
import openpyxl


rad2deg=  180.0/np.pi
deg2rad=np.pi/180.0
mod='RungeKutta2'
# mod='ExplicitEuler'


#this program runs the total simulation for the buoy and 
#plots the simulated buoy positions and observations in the folder
# It also computes the error statistics which are printed in the log file.

def simulate(s,Bnum,indexing,Cor):
 [x0,times,PD]=im.initialisation(Bnum,indexing,s)
 Uavec=PD['Uavec']; Utvec=PD['Utvec']; Uovec=PD['Uovec'];Pgvec=PD['Pgvec'];Pgtvec=PD['Pgtvec']
 h=s['h']
 Hivec=PD['Hivec']
 omega=s['omega']
 tmplier=s['tmplier']
 Xis=lon=x0[0]
 Yis=lat=x0[1]
 x=x0 
 x[0]=x[1]=0 #initial values in SI unit assuming the starting latitude as 0m.
 results=np.zeros((len(x),len(times)))
 ti=0
 # loop
 for t in times:
  # if Cor==1: #Old logic for coriolis
      # f = 2.0*omega*np.sin(lat*deg2rad)
  # else:
  #  f=0
  results[:,ti]=x[:]
  f = 2.0*omega*np.sin(lat*deg2rad)
  tiv=int(np.floor(ti*tmplier))
  Ua=Uavec[tiv,0];Va=Uavec[tiv,1]
  Ut=Utvec[tiv,0];Vt=Utvec[tiv,1]
  Uo=Uovec[tiv,0];Vo=Uovec[tiv,1]
  Pgx=Pgvec[tiv,0];Pgy=Pgvec[tiv,1]
  Pgxt=Pgtvec[tiv,0];Pgyt=Pgtvec[tiv,1]
  #for ice thickness from some model (variable)
  # if Hivec[tiv]>0.2:
  #   h=Hivec[tiv]
  # else:
  #   h=0.2
  # h=0.5
  consts=[f,h,Ua, Va, Ut, Vt, Uo, Vo,Pgx,Pgy,Pgxt,Pgyt]
  consts=np.multiply(Cor,consts)
  if (mod=='ExplicitEuler'):
   xn=im.expeumodel(x,consts,s)
  else:
   xn=im.rk2model(x,consts,s) 
  x  = xn
  #coriolis calculator
  dx=xn[0]-results[0,ti]
  dy=xn[1]-results[1,ti]
  #print dx,dy
  (dlon,dlat)=gf.m2latlon(lat,dx,dy)
  lat=lat+dlat
  lon=lon+dlon
  Xis=np.append(Xis,lon)
  Yis=np.append(Yis,lat)
  ti+=1

 #simulated ice positions and velocity #in meters
 xis=results[0,:]; yis=results[1,:]
 Uis=results[2,:]; Vis=results[3,:]
 # Nu=results[4,:]; Nv=results[5,:]
 Uisvec=np.column_stack((Uis,Vis))
 return(Xis,Yis,Uisvec,results,times,PD)
 
def body(Bnum,indexing,numtaps,Cor):
  prefix="BUOY_"
  bname=prefix+Bnum
  path = "../generated_data/"+bname
  (Cornam,folname)=gf.Cordesfunc(Cor)
  path=path+'/'+folname
  gf.mkdir_p(path)
  fedge=int(numtaps/2)
  logging.info("Model simulation started.")
  s=settings.settings()
  # s['tmplier']=1/30
  # s['dt'] = 15*s['tmplier']*s['minut
  if (mod=='ExplicitEuler'):
    logging.info("Using Explicit Euler for simulation.")
  else:
    logging.info("Using Runge Kutta 2 for simulation.")
  [Xis,Yis,Uisvec,results,times,PD]=simulate(s,Bnum,indexing,Cor)

  tmplierinv=int(1/s['tmplier']) 
  logging.info("Model Simulations done.")
  Xib=PD['Xib'][:-fedge+1]; Yib=PD['Yib'][:-fedge+1]
  Xibf=PD['Xibf']; Yibf=PD['Yibf']
  Uibvec=PD['Uibvec']
  logging.info("Plotting started")
  gp.plticevel(Uisvec,Uibvec,tmplierinv,path)
  gp.plticepos(Xib,Yib,Xis,Yis,path)
  logging.info("Plotting completed. Files available in:" +path)
  #error statistics
  
  #pos error stats
  (merr,rms,werr)=gf.errstats(Xib,Yib,Xis,Yis,tmplierinv)
  logging.info("Mean error in position is: "+str(merr))
  logging.info("RMS error in position is: "+str(rms))
  logging.info("Weighted mean error in position is: "+str(werr))
  errpos=np.column_stack((merr,rms,werr))
  #vel error stats
  Uib=Uibvec[:,0][:-fedge];Vib=Uibvec[:,1][:-fedge]
  Uis=Uisvec[:,0];Vis=Uisvec[:,1]
  (merr,rms,werr)=gf.errstats(Uib,Vib,Uis,Vis,tmplierinv)
  logging.info("Mean error in velocity is: "+str(merr))
  logging.info("RMS error in velocity is: "+str(rms))
  logging.info("Weighted mean error in velocity is: "+str(werr))
  errvel=np.column_stack((merr,rms,werr))
  #creating excel file
  simpost2excel(path,bname,Xis,Yis,Cornam,errpos,errvel)
  logging.info("Processing completed for buoy: "+ bname )
  # gf.logcopy(path)
  logging.shutdown()

def poserrormodel(Xib,Yib,Xis,Yis,tmplierinv):
 Xis=Xis[::tmplierinv]
 Yis=Yis[::tmplierinv]
 print(len(Xis))
 print(len(Xib))
 dX=np.sqrt(((Xib - Xis) ** 2).mean())
 dY=np.sqrt(((Yib - Yis) ** 2).mean())
 rmspos=[dX,dY]
 return(rmspos)

def simpost2excel(path,bname,Xis,Yis,Cornam,errpos,errvel):
 fnp=['Coriolis','Ice thickness','x Air Velocity','y Air Velocity', 
     'x Tidal Velocity','y Tidal Velocity','x Ocean Velocity', 
     'y Ocean Velocity','x-PGs Ocean','y-PGs Ocean','x-PGs tides','y-PGs tides']
 writer = pd.ExcelWriter(path+'/simdata.xlsx', engine='xlsxwriter') 
 dfxis= pd.DataFrame({'Xis': Xis})
 dfyis= pd.DataFrame({'Yis': Yis})
 dfb=pd.DataFrame({bname:[]})
 dfe=pd.DataFrame({'Error statistics':[]})
 dff=pd.DataFrame({'Model F&P':fnp})
 dffv=pd.DataFrame({'Values':Cornam})
 dfep=pd.DataFrame({'Position': ['Mean Error','RMS Error','Weighted Mean Error']})
 dfepx=pd.DataFrame({'x': errpos[0,:]})
 dfepy=pd.DataFrame({'y': errpos[1,:]})
 dfev=pd.DataFrame({'Velocity': ['Mean Error','RMS Error','Weighted Mean Error']})
 dfevx=pd.DataFrame({'x': errvel[0,:]})
 dfevy=pd.DataFrame({'y': errvel[1,:]})
 dfb.to_excel(writer,'Sheet1',startcol=0,startrow=0,index=False)
 dfxis.to_excel(writer,'Sheet1', startcol=0,startrow=2,index=False)
 dfyis.to_excel(writer,'Sheet1', startcol=1,startrow=2,index=False)
 dfe.to_excel(writer,'Sheet1',startcol=6,startrow=4,index=False)
 dff.to_excel(writer,'Sheet1', startcol=3,startrow=2,index=False)
 dffv.to_excel(writer,'Sheet1',startcol=4,startrow=2,index=False)
 dfep.to_excel(writer,'Sheet1',startcol=6,startrow=5,index=False)
 dfepx.to_excel(writer,'Sheet1',startcol=7,startrow=5,index=False)
 dfepy.to_excel(writer,'Sheet1',startcol=8,startrow=5,index=False)
 dfev.to_excel(writer,'Sheet1',startcol=6,startrow=9,index=False)
 dfevx.to_excel(writer,'Sheet1',startcol=7,startrow=9,index=False)
 dfevy.to_excel(writer,'Sheet1',startcol=8,startrow=9,index=False)
 workbook  = writer.book
 worksheet = writer.sheets['Sheet1']
 merge_format = workbook.add_format({
    'bold': 1,
    'border': 1,
    'align': 'center',
    'valign': 'vcenter',
     })
 workbook.close()





#code snippet for running convergence script for ice modelling.

  # tmvec=np.array([1,1/3,1/5,1/15])
  # rmsposvec=[]
  # for i in range(len(tmvec)):
  #  s['tmplier']=tmvec[i]
  #  s['dt'] = 15*s['tmplier']*s['minutes']
  #  [Xis,Yis,Uisvec,results,PD]=simulate(s,Bnum,indexing,Cor)
  #  print(len(Xis))
  #  tmplier30=int(30*s['tmplier']) # 30 or 1 30 because standard model is 30sec 1 because obs is 15mins
  #  Xib=PD['Xib']                         
  #  Yib=PD['Yib']
  #  # print(Xis)
     # rmspos=poserrormodel(Xisc,Yisc,Xis,Yis,tmplier30)
  # rmsposvec=np.append(rmsposvec,rmspos)
  # print(rmsposvec)
    # # gp.convplt(rmsposvec,tmvec,s,path,mod)
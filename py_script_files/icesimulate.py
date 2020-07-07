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
# mod='geostromodel'

#this program runs the total simulation for the buoy and 
#plots the simulated buoy positions and observations in the folder
# It also computes the error statistics which are printed in the log file.

def simulate(s,Bnum,indexing,Cor):
 [x0,times,PD]=im.initialisation(Bnum,indexing,s)
 Uavec=PD['Uavec']; Utvec=PD['Utvec']; Uovec=PD['Uovec'];Pgvec=PD['Pgvec'];Pgtvec=PD['Pgtvec']
 ho=s['h'];trate=s['trate'];Yib=PD['Yib']
 Hivec=PD['Hivec']
 omega=s['omega']
 tmplier=s['tmplier']
 Xis=lon=x0[0]
 Yis=lat=x0[1]
 x=x0 
 x[0]=x[1]=0 #initial values in SI unit assuming the starting latitude as 0m.
 results=np.zeros((len(x),len(times)));hvec=[]
 ti=0
 # loop
 for t in times:
  # if Cor==1: #Old logic for coriolis
      # f = 2.0*omega*np.sin(lat*deg2rad)
  # else:
  #  f=0
  results[:,ti]=x[:]
  tiv=int(np.floor(ti*tmplier))
  f = 2.0*omega*np.sin(Yib[tiv]*deg2rad)
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
  #ice thickness with thinning rate (trate) trate=0 for constat thickness 
  if ((ti*tmplier)%1==0):
    h=gf.thinrate(ho,trate)
    hvec=np.append(hvec,h*Cor[1]);ho=h
  consts=[f,h,Ua, Va, Ut, Vt, Uo, Vo,Pgx,Pgy,Pgxt,Pgyt]
  consts=np.multiply(Cor,consts)
  if (mod=='ExplicitEuler'):
   xn=im.expeumodel(x,consts,s)
  elif (mod=='RungeKutta2'):
   xn=im.rk2model(x,consts,s) 
  else:
   xn=im.geostromodel(x,consts,s)
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
 PD['hvec']=hvec
 return(Xis,Yis,Uisvec,results,times,PD)
 
def body(Bnum,indexing,numtaps,Cor):
  prefix="BUOY_"
  bname=prefix+Bnum
  path = "../../generated_data/"+bname
  s=settings.settings()
  h=s['h'];trate=s['trate']
  (Cornam,folname)=gf.Cordesfunc(Cor,trate,h)
  path=path+'/'+folname
  gf.mkdir_p(path)
  fedge=int(numtaps/2)
  logging.info("Model simulation started.")
  Corn=[]
  if Cor[1]!='v':
    Corn=Cor
    s['trate']=0
    logging.info("Running with constant ice thickness.")
  else:
    Corn=Cor[:]
    Corn[1]=1
  if (mod=='ExplicitEuler'):
    logging.info("Using Explicit Euler for simulation.")
  else:
    logging.info("Using Runge Kutta 2 for simulation.")

  [Xis,Yis,Uisvec,results,times,PD]=simulate(s,Bnum,indexing,Corn)

  tmplierinv=int(1/s['tmplier']) 
  logging.info("Model Simulations done.")
  Xib=PD['Xib'][:-fedge+1]; Yib=PD['Yib'][:-fedge+1]
  Xibf=PD['Xibf']; Yibf=PD['Yibf']
  Uibvec=PD['Uibvec'];hvec=PD['hvec'];Tib=PD['Tib']
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
  # #low pass plots
  # [Xsres,Xis1,xisfil]=gf.LPfilter (numtaps,Xis)
  # gp.pltfilsig(Xis1,xisfil,Xsres,path,'sim_lon_filt')
  # [Ysres,Yis1,yisfil]=gf.LPfilter (numtaps,Yis)
  # gp.pltfilsig(Yis1,yisfil,Ysres,path,'sim_lat_filt')
  # Fourier Transforms  
  # Longitude  
  (tvec,xbres,xsres,Xbfil,Xsfil,arg_tide,arg_cor,errftamlon,errftphlon,tidblon,tidslon)=gf.FTremMD(numtaps,Xib,Xis,tmplierinv)
  (tvec,ubres,usres,Ubfil,Usfil,arg_tide,arg_cor,errftamlonv,errftphlonv,tidbu,tidsu)=gf.FTremMD(numtaps,Uib,Uis,tmplierinv)
  gp.pltFT(path,"Longitude",xbres,xsres,tvec,arg_tide,arg_cor)
  
  # Latitude
  (tvec,ybres,ysres,Ybfil,Ysfil,arg_tide,arg_cor,errftamlat,errftphlat,tidblat,tidslat)=gf.FTremMD(numtaps,Yib,Yis,tmplierinv)
  (tvec,vbres,vsres,Vbfil,Vsfil,arg_tide,arg_cor,errftamlonv,errftphlonv,tidbv,tidsv)=gf.FTremMD(numtaps,Vib,Vis,tmplierinv)
  gp.pltFT(path,"Latitude",ybres,ysres,tvec,arg_tide,arg_cor)
  logging.info("Fourier Transforms plotted.")

  # error computations
  err_ft=np.row_stack((errftamlon,errftamlat,errftphlon,errftphlat))
  #tide computations
  tidb=np.row_stack((tidblon[0,:],tidblat[0,:],tidblon[1,:],tidblat[1,:]))
  tids=np.row_stack((tidslon[0,:],tidslat[0,:],tidslon[1,:],tidslat[1,:]))
  #creating excel file
  # Xdata={'Xis':Xis,'Xbfil':Xbfil,'Xsfil':Xsfil,'Usfil':Usfil,'Ubfil':Ubfil}
  # Ydata={'Yis':Yis,'Ybfil':Ybfil,'Ysfil':Ysfil,'Vsfil':Vsfil,'Vbfil':Vbfil}
  Xdata={'Xis':Xis,'Xbfil':Xbfil,'Xsfil':Xsfil,'Usfil':Uis[::tmplierinv],'Ubfil':Ubfil}
  Ydata={'Yis':Yis,'Ybfil':Ybfil,'Ysfil':Ysfil,'Vsfil':Vis[::tmplierinv],'Vbfil':Vbfil}

  simpost2excel(path,bname,Xdata,Ydata,hvec,Cornam,errpos,errvel,err_ft,tidb,tids)
  logging.info("Excel data file created for simulated data.")
  logging.info("Processing completed for buoy: "+ bname )
  gf.logcopy(path)
  logging.shutdown()


def simpost2excel(path,bname,Xdata,Ydata,hvec,Cornam,errpos,errvel,err_ft,tidb,tids):
 fnp=['Coriolis','Ice thickness','x Air Velocity','y Air Velocity', 
     'x Tidal Velocity','y Tidal Velocity','x Ocean Velocity', 
     'y Ocean Velocity','x-PGs Ocean','y-PGs Ocean','x-PGs tides','y-PGs tides']
 writer = pd.ExcelWriter(path+'/simdata_'+bname+'.xlsx', engine='xlsxwriter') 
 dfxis= pd.DataFrame({'Xis': Xdata['Xis']})
 dfyis= pd.DataFrame({'Yis': Ydata['Yis']})
 dfxsf= pd.DataFrame({'Xsfil': Xdata['Xsfil']})
 dfysf= pd.DataFrame({'Ysfil': Ydata['Ysfil']})
 dfxbf= pd.DataFrame({'Xbfil': Xdata['Xbfil']})
 dfybf= pd.DataFrame({'Ybfil': Ydata['Ybfil']})
 dfubf= pd.DataFrame({'Ubfil': Xdata['Ubfil']})
 dfvbf= pd.DataFrame({'Vbfil': Ydata['Vbfil']})
 dfusf= pd.DataFrame({'Usfil': Xdata['Usfil']})
 dfvsf= pd.DataFrame({'Vsfil': Ydata['Vsfil']})
 dfh=pd.DataFrame({'Ice thickness':hvec})
 dfb=pd.DataFrame({bname:[]})
 dfe=pd.DataFrame({'Error statistics':[]})
 dff=pd.DataFrame({'Model F&P':fnp})
 dffv=pd.DataFrame({'Values':Cornam})
 dfep=pd.DataFrame({'Position': ['Mean Error','RMS Error','Weighted Mean Error']})
 dfepx=pd.DataFrame({'x': errpos[0,:]})
 dfepy=pd.DataFrame({'y': errpos[1,:]})
 dfepa=pd.DataFrame({'absolute': errpos[2,:]})
 dfev=pd.DataFrame({'Velocity': ['Mean Error','RMS Error','Weighted Mean Error']})
 dfevx=pd.DataFrame({'x': errvel[0,:]})
 dfevy=pd.DataFrame({'y': errvel[1,:]})
 dfeva=pd.DataFrame({'absolute': errvel[2,:]})
 dfeft=pd.DataFrame({'Tide': ['M2','S2','MU2','O1','K1','M4']})
 dfeftlo=pd.DataFrame({'Err_Lon_amp': err_ft[0,:]})
 dfeftla=pd.DataFrame({'Err_Lat_amp': err_ft[1,:]})
 dfeftphlo=pd.DataFrame({'Err_Lon_ph': err_ft[2,:]})
 dfeftphla=pd.DataFrame({'Err_Lat_ph': err_ft[3,:]})
 dftidbalo=pd.DataFrame({'B_Lon_amp': tidb[0,:]})
 dftidbala=pd.DataFrame({'B_Lat_amp': tidb[1,:]})
 dftidbplo=pd.DataFrame({'B_Lon_ph': tidb[2,:]})
 dftidbpla=pd.DataFrame({'B_Lat_ph': tidb[3,:]})
 dftidsalo=pd.DataFrame({'M_Lon_amp': tids[0,:]})
 dftidsala=pd.DataFrame({'M_Lat_amp': tids[1,:]})
 dftidsplo=pd.DataFrame({'M_Lon_ph': tids[2,:]})
 dftidspla=pd.DataFrame({'M_Lat_ph': tids[3,:]})
 dfb.to_excel(writer,'Sheet1',startcol=0,startrow=0,index=False)
 dfxis.to_excel(writer,'Sheet1', startcol=0,startrow=2,index=False)
 dfyis.to_excel(writer,'Sheet1', startcol=1,startrow=2,index=False)
 dfxsf.to_excel(writer,'Sheet1', startcol=2,startrow=2,index=False)
 dfysf.to_excel(writer,'Sheet1', startcol=3,startrow=2,index=False) 
 dfxbf.to_excel(writer,'Sheet1', startcol=4,startrow=2,index=False)
 dfybf.to_excel(writer,'Sheet1', startcol=5,startrow=2,index=False)
 dfubf.to_excel(writer,'Sheet1', startcol=6,startrow=2,index=False)
 dfvbf.to_excel(writer,'Sheet1', startcol=7,startrow=2,index=False) 
 dfusf.to_excel(writer,'Sheet1', startcol=8,startrow=2,index=False)
 dfvsf.to_excel(writer,'Sheet1', startcol=9,startrow=2,index=False)
 dfh.to_excel(writer,'Sheet1', startcol=10,startrow=2,index=False)
 dfe.to_excel(writer,'Sheet1',startcol=17,startrow=4,index=False)
 dff.to_excel(writer,'Sheet1', startcol=14,startrow=2,index=False)
 dffv.to_excel(writer,'Sheet1',startcol=15,startrow=2,index=False)
 dfep.to_excel(writer,'Sheet1',startcol=17,startrow=5,index=False)
 dfepx.to_excel(writer,'Sheet1',startcol=18,startrow=5,index=False)
 dfepy.to_excel(writer,'Sheet1',startcol=19,startrow=5,index=False)
 dfepa.to_excel(writer,'Sheet1',startcol=20,startrow=5,index=False)
 dfev.to_excel(writer,'Sheet1',startcol=17,startrow=9,index=False)
 dfevx.to_excel(writer,'Sheet1',startcol=18,startrow=9,index=False)
 dfevy.to_excel(writer,'Sheet1',startcol=19,startrow=9,index=False)
 dfeva.to_excel(writer,'Sheet1',startcol=20,startrow=9,index=False)
 dfeft.to_excel(writer,'Sheet1',startcol=17,startrow=14,index=False)
 dfeftlo.to_excel(writer,'Sheet1',startcol=18,startrow=14,index=False)
 dfeftla.to_excel(writer,'Sheet1',startcol=19,startrow=14,index=False)
 dfeftphlo.to_excel(writer,'Sheet1',startcol=20,startrow=14,index=False)
 dfeftphla.to_excel(writer,'Sheet1',startcol=21,startrow=14,index=False)
 dftidbalo.to_excel(writer,'Sheet1',startcol=23,startrow=14,index=False)
 dftidbala.to_excel(writer,'Sheet1',startcol=24,startrow=14,index=False)
 dftidbplo.to_excel(writer,'Sheet1',startcol=25,startrow=14,index=False)
 dftidbpla.to_excel(writer,'Sheet1',startcol=26,startrow=14,index=False)
 dftidsalo.to_excel(writer,'Sheet1',startcol=27,startrow=14,index=False)
 dftidsala.to_excel(writer,'Sheet1',startcol=28,startrow=14,index=False)
 dftidsplo.to_excel(writer,'Sheet1',startcol=29,startrow=14,index=False)
 dftidspla.to_excel(writer,'Sheet1',startcol=30,startrow=14,index=False)
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
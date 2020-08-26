## This script runs the model according to the defined parameters
#The model results are then used to compute error statistics and filtering.
# The statistics, model results etc are stored in excel files in the model run folders (eg. h1f1A1T1O1Po1Pt1)
import numpy as np
import generalfunc as gf
import genplots as gp
import logging
import model
import pandas as pd
import settings 
import openpyxl


rad2deg=  180.0/np.pi
deg2rad=np.pi/180.0



def body(Bnum,indexing,forcevec):
  prefix="BUOY_"
  bname=prefix+Bnum
  path = "../../generated_data/"+bname
  PD=model.readposveldata(path)  
  s=settings.settings()
  h=s['h'];trate=s['trate']
  (forcenam,folname)=gf.forcedetail(forcevec,trate,h)
  #creation of the folder for storing the simulated data.
  path=path+'/'+folname
  gf.mkdir_p(path)
  logging.info("Model simulation started.")
  ## description of consant ice thickness or not in the simulation.
  [forcevecn,trate]=gf.icethicktyp(forcevec,trate)
  s['trate']=trate 
  PD=model.simulate(s,Bnum,indexing,forcevecn,PD)
  logging.info("Model Simulations done.")
  Xib=PD['Xib']; Yib=PD['Yib']   
  Uibvec=PD['Uibvec'];hvec=PD['hvec'];Tib=PD['Tib']
  Xis=PD['Xis']; Yis=PD['Yis'];Uisvec=PD['Uisvec']


  ## We will do the plotting part later in postprocessing.
  logging.info("Plotting started")
  # gp.plticevel(Uisvec,Uibvec,tmplierinv,path)
  gp.plticepos(Xib,Yib,Xis,Yis,path)
  logging.info("Plotting completed. Files available in:" +path)

  ## Statistic computation.
  (XD,YD)=statscompute(s,Xib,Yib,Xis,Yis)
  (UD,VD)=statscompute(s,Uibvec[:,0],Uibvec[:,1],Uisvec[:,0],Uisvec[:,1])
  logging.info("Statistics computation done.")
  simdata2excel(path,Bnum,forcenam,XD,YD,UD,VD,PD) 
  logging.info("Excel data file created for simulated data.")
  logging.info("Processing completed for buoy: "+ bname )
  gf.logcopy(path)
  logging.shutdown()




## Computation of statistics and filtering 
#Error statistics function 

def errstats(Xobs,Yobs,Xsim,Ysim):
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

## Filtering and FT

# Filtering followed by Fourier transform
# This gives am and ph of high pass filtered solution.
def filternFT(s,numtaps,Xo,Xs):
  #filtering with lowpas filer
  [Xores,Xos,Xofil]=gf.LPfilter (numtaps,Xo)
  [Xsres,Xss,Xsfil]=gf.LPfilter (numtaps,Xs)
  Nft=len(Xores);dt=15*60.0 #time difference in observations
  [Xoresam,Xoresph,fvec,tvec]=gf.FFTsig(Xores,Nft,dt) 
  [Xsresam,Xsresph,fvec,tvec]=gf.FFTsig(Xsres,Nft,dt)
  Xoft=np.row_stack((Xoresam,Xoresph))
  Xsft=np.row_stack((Xsresam,Xsresph))
  tvec=tvec/3600 #Period in hours
  return(tvec,Xoft,Xsft,Xofil,Xsfil)

def tidcorargcompute(s,tvec): 
  tidedict=s['tidedict'];cordict=s['cordict']
  corargdict={}
  for k in cordict.keys():
    phi=np.deg2rad(cordict[k])
    f=2*s['omega']*np.sin(phi)
    T=2*np.pi/(f*3600)
    arg=(np.abs(tvec-T)).argmin()
    corargdict[k]=arg
  tideargdict={}
  for k in tidedict.keys():
    arg=(np.abs(tvec-tidedict[k])).argmin()
    tideargdict[k]=arg
  return(corargdict,tideargdict)

def tidcorvalcompute(corarg,tidearg,Xoresam,Xsresam):
  tideamo={};tideams={};erram={}
  for k in tidearg.keys():
    tideamo[k]=Xoresam[tidearg[k]]
    tideams[k]=Xsresam[tidearg[k]]
    erram[k]=Xoresam[tidearg[k]]-Xsresam[tidearg[k]]
  return(tideamo,tideams,erram)


#Xo is observed and Xs is simulated.
def statscompute(s,Xo,Yo,Xs,Ys): 
## computation of error statistics  
  (merr,rms,werr)=errstats(Xo,Yo,Xs,Ys)
  logging.info("Mean error in position is: "+str(merr))
  logging.info("RMS error in position is: "+str(rms))
  logging.info("Weighted mean error in position is: "+str(werr))
  errorvec=np.column_stack((merr,rms,werr))
## Filtering of the paths.
  numtaps=2*24*2+1
  fedge=int(numtaps/2)
# Filtering and computation of Fourier transforms.
#   # Longitude  
  [tvec,Xoft,Xsft,Xofil,Xsfil]=filternFT(s,numtaps,Xo,Xs)
  #computation of tidal and coriolis frequency arguments
  [corarg,tidearg]= tidcorargcompute(s,tvec)
  [longtamo,longtams,longerram]=tidcorvalcompute(corarg,tidearg,Xoft[0,:],Xsft[0,:])
  [longtpho,longtphs,longerrph]=tidcorvalcompute(corarg,tidearg,Xoft[1,:],Xsft[1,:])
#   # Latitude 
  [tvec,Yoft,Ysft,Yofil,Ysfil]=filternFT(s,numtaps,Yo,Ys)
  #computation of tidal and coriolis frequency arguments
  [corarg,tidearg]= tidcorargcompute(s,tvec)
  [lattamo,lattams,laterram]=tidcorvalcompute(corarg,tidearg,Yoft[0,:],Ysft[0,:])
  [lattpho,lattphs,laterrph]=tidcorvalcompute(corarg,tidearg,Yoft[1,:],Ysft[1,:])  
  #dictionary to store all data.
  # X dict
  XD={'Xoft':Xoft,'Xofil':Xofil,'Xsft':Xsft,'Xsfil':Xsfil,
      'tamo':longtamo,'tams':longtams,'tame':longerram,
      'tpho':longtpho,'tphs':longtphs,'tphe':longerrph,'err':errorvec}
  YD={'Yoft':Yoft,'Yofil':Yofil,'Ysft':Ysft,'Ysfil':Ysfil,
    'tamo':lattamo,'tams':lattams,'tame':laterram,
    'tpho':lattpho,'tphs':lattphs,'tphe':laterrph}
  return (XD,YD)

## creating excel file of the data


# This function creates an excel file with all the simulation data. 
# then this data can be used for post processing.

def simdata2excel(path,Bnum,forcenam,XD,YD,UD,VD,PD):
  prefix="BUOY_"
  bname=prefix+Bnum
  fnp=['Coriolis','Ice thickness','x Air Velocity','y Air Velocity', 
      'x Tidal Velocity','y Tidal Velocity','x Ocean Velocity', 
      'y Ocean Velocity','x-PGs Ocean','y-PGs Ocean','x-PGs tides','y-PGs tides']
  writer = pd.ExcelWriter(path+'/Simdata_'+bname+'.xlsx', engine='xlsxwriter') 
  dft=pd.DataFrame({'Date(GMT)': PD['Tib']})
  dfxis= pd.DataFrame({'Xis': PD['Xis']})
  dfyis= pd.DataFrame({'Yis': PD['Yis']})
  dfxib= pd.DataFrame({'Xib': PD['Xib']})
  dfyib= pd.DataFrame({'Yib': PD['Yib']})
  dfxsf= pd.DataFrame({'Xsfil': XD['Xsfil']})
  dfysf= pd.DataFrame({'Ysfil': YD['Ysfil']})
  dfxbf= pd.DataFrame({'Xbfil': XD['Xofil']})
  dfybf= pd.DataFrame({'Ybfil': YD['Yofil']})
  dfus= pd.DataFrame({'Uis': PD['Uisvec'][:,0]})
  dfvs= pd.DataFrame({'Vis': PD['Uisvec'][:,1]})  
  dfub= pd.DataFrame({'Uib': PD['Uibvec'][:,0]})
  dfvb= pd.DataFrame({'Vib': PD['Uibvec'][:,1]})  
  dfusf= pd.DataFrame({'Usfil': UD['Xsfil']})
  dfvsf= pd.DataFrame({'Vsfil': VD['Ysfil']})
  dfubf= pd.DataFrame({'Ubfil': UD['Xofil']})
  dfvbf= pd.DataFrame({'Vbfil': VD['Yofil']})
  dfh=pd.DataFrame({'Ice thickness':PD['hvec']})
  dfe=pd.DataFrame({'Error statistics':[]})
  dff=pd.DataFrame({'Model F&P':fnp})
  dffv=pd.DataFrame({'Values':forcenam})
  dfep=pd.DataFrame({'Position': ['Mean Error','RMS Error','Weighted Mean Error']})
  dfepx=pd.DataFrame({'x': XD['err'][0,:]})
  dfepy=pd.DataFrame({'y': XD['err'][1,:]})
  dfepa=pd.DataFrame({'absolute': XD['err'][2,:]})
  dfev=pd.DataFrame({'Velocity': ['Mean Error','RMS Error','Weighted Mean Error']})
  dfevx=pd.DataFrame({'x': UD['err'][0,:]})
  dfevy=pd.DataFrame({'y': UD['err'][1,:]})
  dfeva=pd.DataFrame({'absolute': UD['err'][2,:]})
  dfeft=pd.DataFrame({'Tide': ['M2','S2','MU2','O1','K1','M4']})
  dftamolo=pd.DataFrame({'OAmpLon': np.array(list(XD['tamo'].values()))})
  dftamslo=pd.DataFrame({'SAmpLon': np.array(list(XD['tams'].values()))})
  dftpholo=pd.DataFrame({'OPhLon': np.array(list(XD['tpho'].values()))})
  dftphslo=pd.DataFrame({'SPhLon': np.array(list(XD['tphs'].values()))})
  dftamola=pd.DataFrame({'OAmpLat': np.array(list(YD['tamo'].values()))})
  dftamsla=pd.DataFrame({'SAmpLat': np.array(list(YD['tams'].values()))})
  dftphola=pd.DataFrame({'OPhLat': np.array(list(YD['tpho'].values()))})
  dftphsla=pd.DataFrame({'SPhLat': np.array(list(YD['tphs'].values()))})
  # writing 
  dft.to_excel(writer,'Sheet1', startcol=0,startrow=0,index=False)
  dfxis.to_excel(writer,'Sheet1', startcol=1,startrow=0,index=False)
  dfyis.to_excel(writer,'Sheet1', startcol=2,startrow=0,index=False)
  dfxib.to_excel(writer,'Sheet1', startcol=3,startrow=0,index=False)
  dfyib.to_excel(writer,'Sheet1', startcol=4,startrow=0,index=False)
  dfxsf.to_excel(writer,'Sheet1', startcol=5,startrow=0,index=False)
  dfysf.to_excel(writer,'Sheet1', startcol=6,startrow=0,index=False) 
  dfxbf.to_excel(writer,'Sheet1', startcol=7,startrow=0,index=False)
  dfybf.to_excel(writer,'Sheet1', startcol=8,startrow=0,index=False)
  dfus.to_excel(writer,'Sheet1', startcol=9,startrow=0,index=False)
  dfvs.to_excel(writer,'Sheet1', startcol=10,startrow=0,index=False)
  dfub.to_excel(writer,'Sheet1', startcol=11,startrow=0,index=False)
  dfvb.to_excel(writer,'Sheet1', startcol=12,startrow=0,index=False)
  dfusf.to_excel(writer,'Sheet1', startcol=13,startrow=0,index=False)
  dfvsf.to_excel(writer,'Sheet1', startcol=14,startrow=0,index=False) 
  dfubf.to_excel(writer,'Sheet1', startcol=15,startrow=0,index=False)
  dfvbf.to_excel(writer,'Sheet1', startcol=16,startrow=0,index=False)
  dfh.to_excel(writer,'Sheet1', startcol=17,startrow=0,index=False)
  dff.to_excel(writer,'Sheet1', startcol=19,startrow=0,index=False)
  dffv.to_excel(writer,'Sheet1',startcol=20,startrow=0,index=False)
  dfe.to_excel(writer,'Sheet1',startcol=22,startrow=2,index=False)
  dfep.to_excel(writer,'Sheet1',startcol=22,startrow=3,index=False)
  dfepx.to_excel(writer,'Sheet1',startcol=23,startrow=3,index=False)
  dfepy.to_excel(writer,'Sheet1',startcol=24,startrow=3,index=False)
  dfepa.to_excel(writer,'Sheet1',startcol=25,startrow=3,index=False)
  dfev.to_excel(writer,'Sheet1',startcol=22,startrow=7,index=False)
  dfevx.to_excel(writer,'Sheet1',startcol=23,startrow=7,index=False)
  dfevy.to_excel(writer,'Sheet1',startcol=24,startrow=7,index=False)
  dfeva.to_excel(writer,'Sheet1',startcol=25,startrow=7,index=False)
  dfeft.to_excel(writer,'Sheet1',startcol=22,startrow=12,index=False)
  dftamolo.to_excel(writer,'Sheet1',startcol=23,startrow=12,index=False)
  dftamslo.to_excel(writer,'Sheet1',startcol=24,startrow=12,index=False)
  dftpholo.to_excel(writer,'Sheet1',startcol=25,startrow=12,index=False)
  dftphslo.to_excel(writer,'Sheet1',startcol=26,startrow=12,index=False)
  dftamola.to_excel(writer,'Sheet1',startcol=27,startrow=12,index=False)
  dftamsla.to_excel(writer,'Sheet1',startcol=28,startrow=12,index=False)
  dftphola.to_excel(writer,'Sheet1',startcol=29,startrow=12,index=False)
  dftphsla.to_excel(writer,'Sheet1',startcol=30,startrow=12,index=False)
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
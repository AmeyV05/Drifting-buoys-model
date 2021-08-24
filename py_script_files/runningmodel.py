## This script runs the model according to the defined parameters
#The model results are then used to compute error statistics and filtering.
# The statistics, model results etc are stored in excel files in the model run folders (eg. h1f1A1T1O1Po1Pt1)
import numpy as np
import pandas as pd
import generalfunc as gf
import plots
import logging
import model
import settings 


def run(Bnum,indexing):
  prefix="BUOY_"
  bname=prefix+Bnum
  path = "../../generated_data/"+bname
  PD=model.readposveldata(path)  
  s=settings.settings()
  h=settings.h;trate=s['trate'];forcevec=settings.forcevec
  (forcenam,folname)=gf.forcedetail(forcevec,trate,h)
  #creation of the folder for storing the simulated data.
  path=path+'/'+folname
  gf.mkdir_p(path)
  ## description of consant ice thickness or not in the simulation.
  [forcevecn,trate]=gf.icethicktyp(forcevec,trate)
  s['trate']=trate 
  logging.info("Model simulation started.")
  PD=model.simulate(s,Bnum,indexing,forcevecn,PD)
  logging.info("Model Simulations done.")
  Xib=PD['Xib']; Yib=PD['Yib']   
  Uibvec=PD['Uibvec'];hvec=PD['hvec'];Tib=PD['Tib']
  Xis=PD['Xis']; Yis=PD['Yis'];Uisvec=PD['Uisvec']
  ## We will do elaborate plotting part later in postprocessing.
  logging.info("Plotting started")
  plots.plticepos(Xib,Yib,Xis,Yis,path,Bnum)
  logging.info("Plotting completed for the buoy path.")
  ## Statistic computation.
  (XD,YD)=statscompute(s,Xib,Yib,Xis,Yis)
  (UD,VD)=statscompute(s,Uibvec[:,0],Uibvec[:,1],Uisvec[:,0],Uisvec[:,1])
  logging.info("Statistics computation done.")
  # simdata2excel(path,Bnum,forcenam,XD,YD,UD,VD,PD) 
  simdata2excl(s,path,Bnum,forcenam,XD,YD,UD,VD,PD)
  logging.info("Excel data file created for simulated data.")
  logging.info("Processing completed for buoy: "+ bname )
  gf.logcopy(path)
  logging.shutdown()


## Functions
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
  [corargX,tideargX]= tidcorargcompute(s,tvec)
  [longtamo,longtams,longerram]=tidcorvalcompute(corargX,tideargX,Xoft[0,:],Xsft[0,:])
  [longtpho,longtphs,longerrph]=tidcorvalcompute(corargX,tideargX,Xoft[1,:],Xsft[1,:])
#   # Latitude 
  [tvec,Yoft,Ysft,Yofil,Ysfil]=filternFT(s,numtaps,Yo,Ys)
  #computation of tidal and coriolis frequency arguments
  [corargY,tideargY]= tidcorargcompute(s,tvec)
  [lattamo,lattams,laterram]=tidcorvalcompute(corargY,tideargY,Yoft[0,:],Ysft[0,:])
  [lattpho,lattphs,laterrph]=tidcorvalcompute(corargY,tideargY,Yoft[1,:],Ysft[1,:])  
  #dictionary to store all data.
  # X dict
  XD={'Xoft':Xoft,'Xofil':Xofil,'Xsft':Xsft,'Xsfil':Xsfil,
      'tamo':longtamo,'tams':longtams,'tame':longerram,
      'tpho':longtpho,'tphs':longtphs,'tphe':longerrph,'err':errorvec,'corarg':corargX,'tidearg':tideargX,'tft':tvec}
      # Y dict
  YD={'Yoft':Yoft,'Yofil':Yofil,'Ysft':Ysft,'Ysfil':Ysfil,
    'tamo':lattamo,'tams':lattams,'tame':laterram,
    'tpho':lattpho,'tphs':lattphs,'tphe':laterrph,'corarg':corargY,'tidearg':tideargY}
  return (XD,YD)

## creating excel file of the data


# This function creates an excel file with all the simulation data. 
# then this data can be used for post processing.

def valuewriter(writer,dname,dvalue,sheetname,i,j):
  dft=pd.DataFrame({dname:dvalue})
  dft.to_excel(writer,sheetname, startcol=i,startrow=j,index=False)



def simdata2excl(s,path,Bnum,forcenam,XD,YD,UD,VD,PD):
  prefix="BUOY_"
  bname=prefix+Bnum
  fnp=['Coriolis','Ice thickness','x Air Velocity','y Air Velocity', 
      'x Tidal Velocity','y Tidal Velocity','x Ocean Velocity', 
      'y Ocean Velocity','x-PGs Ocean','y-PGs Ocean','x-PGs tides','y-PGs tides']
  dnamevecs1=['Date(GMT)','Xis','Yis','Xib','Yib','Xsfil','Ysfil','Xbfil','Ybfil',
            'Uis','Vis','Uib','Vib','Usfil','Vsfil','Ubfil','Vbfil','Ice thickness',
            'Model F&P','Values','Error statistics','Position','x','y','absolute','Velocity','x','y','absolute']
  dvaluevecs1=[PD['Tib'],PD['Xis'],PD['Yis'],PD['Xib'],PD['Yib'],XD['Xsfil'],YD['Ysfil'],XD['Xofil'],YD['Yofil'],
            PD['Uisvec'][:,0],PD['Uisvec'][:,1],PD['Uibvec'][:,0],PD['Uibvec'][:,1],UD['Xsfil'],VD['Ysfil'],UD['Xofil'],VD['Yofil'],
            PD['hvec'],fnp,forcenam,[],['Mean Error','RMS Error','Weighted Mean Error'],XD['err'][0,:],XD['err'][1,:],XD['err'][2,:],
            ['Mean Error','RMS Error','Weighted Mean Error'],UD['err'][0,:],UD['err'][1,:], UD['err'][2,:]]
  colvecs1=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,22,22,23,24,25,22,23,24,25]
  rowvecs1=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,3,3,3,3,7,7,7,7]
  writer = pd.ExcelWriter(path+'/Simdata_'+bname+'.xlsx',engine='xlsxwriter') 
  for i in range(len(dnamevecs1)):
    valuewriter(writer,dnamevecs1[i],dvaluevecs1[i],'Model_Data',colvecs1[i],rowvecs1[i])

  dnamevecs2=['Xsam','Ysam','Xsph','Ysph','Xbam','Ybam','Xbph','Ybph',
              'Usam','Vsam','Usph','Vsph','Ubam','Vbam','Ubph','Vbph','Tft',
              'Tide','BAmpLon','SAmpLon','BPhLon','SPhLon',
              'BAmpLat','SAmpLat','BPhLat','SPhLat','tideargLon','tideargLat','corarglon','corarglat',
              'BAmpU','SAmpU','BPhU','SPhU',
              'BAmpV','SAmpV','BPhV','SPhV','tideargU','tideargV','corargU','corargV']
  dvaluevecs2=[XD['Xsft'][0,:],YD['Ysft'][0,:],XD['Xsft'][1,:],YD['Ysft'][1,:],
              XD['Xoft'][0,:],YD['Yoft'][0,:],XD['Xoft'][1,:],YD['Yoft'][1,:],
              UD['Xsft'][0,:],VD['Ysft'][0,:],UD['Xsft'][1,:],VD['Ysft'][1,:],
              UD['Xoft'][0,:],VD['Yoft'][0,:],UD['Xoft'][1,:],VD['Yoft'][1,:],XD['tft'],
              np.array(list(s['tidedict'].keys())),np.array(list(XD['tamo'].values())),
              np.array(list(XD['tams'].values())),np.array(list(XD['tpho'].values())),
              np.array(list(XD['tphs'].values())),np.array(list(YD['tamo'].values())),
              np.array(list(YD['tams'].values())),np.array(list(YD['tpho'].values())),
              np.array(list(YD['tphs'].values())),np.array(list(XD['tidearg'].values())),
              np.array(list(YD['tidearg'].values())),np.array(list(XD['corarg'].values())),
              np.array(list(YD['corarg'].values())),np.array(list(UD['tamo'].values())),
              np.array(list(UD['tams'].values())),np.array(list(UD['tpho'].values())),
              np.array(list(UD['tphs'].values())),np.array(list(VD['tamo'].values())),
              np.array(list(VD['tams'].values())),np.array(list(VD['tpho'].values())),
              np.array(list(VD['tphs'].values())),np.array(list(UD['tidearg'].values())),
              np.array(list(VD['tidearg'].values())),np.array(list(UD['corarg'].values())),
              np.array(list(VD['corarg'].values()))]
  colvecs2=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,
            30,31,32,33,34,35,36,37,38,39,40,41]
  rowvecs2=(np.zeros(len(colvecs2),dtype=int))
  for i in range(len(dnamevecs2)):
    valuewriter(writer,dnamevecs2[i],dvaluevecs2[i],'FT_Data',colvecs2[i],rowvecs2[i])
  # writer.save()
  workbook  = writer.book
  # worksheet1 = writer.sheets['Model_Data']
  # worksheet2 = writer.sheets['FT_Data']
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
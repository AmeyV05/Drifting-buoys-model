import cartopy.crs as ccrs
import cartopy.feature as cpf 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import scipy.io
from scipy.interpolate import interp2d
from scipy.interpolate import Rbf
import Bdataprocess as bdat
import drifticemodel as dimodel
import os
import generalfunc as gf
import genplots as gp
import logging

def main(Bnum,indexing,Cor,TD,AD,OD):
    prefix="BUOY_"
    bname=prefix+Bnum
    logging.info("Processing started for buoy: "+ bname)
    #Buoy data analysis
    # define the name of the directory to be created
    path = "../generated_data/"+bname
    gf.mkdir_p(path)
    file ='../buoy_data_excel/procesd_buoy_loc_data_Lars/'+bname+'.csv'
    #file ='../buoy_data_excel/CliSAP_Boje_16_modified.csv'
    [Xib,Yib,Uib,Vib,Tib]=bdat.bdata(file)
    Xib=Xib[96:]
    Yib=Yib[96:]
    Uib=Uib[96:]
    Vib=Vib[96:]
    Tib=Tib[96:]
    # print(len(Xib))
    logging.info("Running simulation for buoy number: " + Bnum)
    # obtaining GTSM tidal data in vector form from dictonary TD
    Xt=TD['Xt']
    Yt=TD['Yt']
    ut=TD['ut']
    vt=TD['vt']
    Tt=TD['Tt']
    sindex=0
    eindex=sindex+indexing
    arrlen=eindex-sindex   
    Tt=gf.num2datetimesecs(2014,3,1,sindex,eindex,Tt)
    logging.info("Interpolating GTSM data to buoy locations....")
    # processing of gtsm data
    #interpolating the u and v velocities to the buoy locations.
    Buoy_Xi=Xib #buoy longitude in hourly basis (x position)
    Buoy_Yi=Yib
    Ut=[]
    Vt=[]
    for i in range(arrlen):
        ut_time=ut[i,:]
        vt_time=vt[i,:]
        xb=Buoy_Xi[i]
        yb=Buoy_Yi[i]
        Xtb=[]
        Ytb=[]
        Utb=[]
        Vtb=[]
        for j in range(len(Xt)):
         if (Xt[j]<=xb+0.2 and Xt[j]>=xb-0.2 and Yt[j]>=yb-0.2 and Yt[j]<=yb+0.2):
          Xtb=np.append(Xtb,Xt[j])
          Ytb=np.append(Ytb,Yt[j])
          Utb=np.append(Utb,ut_time[j])
          Vtb=np.append(Vtb,vt_time[j])
        #uinterp=interp2d(Xtb,Ytb,Utb,kind='linear')
        #vinterp=interp2d(Xtb,Ytb,Vtb,kind='linear')
        uinterp=Rbf(Xtb,Ytb,Utb)
        vinterp=Rbf(Xtb,Ytb,Vtb)
        #creating winds at buoy locations
        Ut=np.append(Ut,uinterp(xb,yb))
        Vt=np.append(Vt,vinterp(xb,yb))
        #print(i)

    logging.info("Interpolation successfully done. Processing complete for GTSM data.")

    sindex=360 #converting to 15th march 00
    eindex=sindex+int((indexing/4))
    arrlen=eindex-sindex  
    # obtaining GTSM tidal data in vector form from dictonary TD
    Xa=AD['Xa']
    Ya=AD['Ya']
    u10=AD['u10']
    v10=AD['v10']
    Ta=AD['Ta']
    Ta=gf.num2datetimehrs(1900,1,1,sindex,eindex,Ta)
    u10=u10[sindex:eindex,:,:]
    v10=v10[sindex:eindex,:,:]
    logging.info("Interpolating wind data to buoy locations....")
    # processing of wind data
    #interpolating the u and v velocities to the buoy locations. Note that the ERA5 winds are hourly. 
    Buoy_Xi=Xib[::4] #buoy longitude in hourly basis (x position)
    Buoy_Yi=Yib[::4]
    Ua=[]
    Va=[]
    for i in range(arrlen):
        u10_time=u10[i,:,:]
        v10_time=v10[i,:,:]
        
        uinterp=interp2d(Xa,Ya,u10_time,kind='linear')
        vinterp=interp2d(Xa,Ya,v10_time,kind='linear')
        #creating winds at buoy locations
        Ua=np.append(Ua,uinterp(Buoy_Xi[i],Buoy_Yi[i]))
        Va=np.append(Va,vinterp(Buoy_Xi[i],Buoy_Yi[i]))

    logging.info("Interpolation successfully done.")
    logging.info("Processing complete for wind data.")

    #Ocean_currents_CMEMS

    sindex=0
    eindex=int(indexing/4)
    arrlen=eindex-sindex
    Xo=OD['Xo']
    Yo=OD['Yo']
    uo=OD['uo']
    vo=OD['vo']
    To=OD['To']
    To=gf.num2datetimehrs(1950,1,1,sindex,eindex,To)
    uo=uo[sindex:eindex,0,:,:]
    vo=vo[sindex:eindex,0,:,:]
    #processing of ocean data 

    logging.info("Interpolating ocean currents data to buoy locations....")
    Buoy_Xi=Xib[::4] #buoy longitude in hourly basis (x position)
    Buoy_Yi=Yib[::4]
    Uo=[]
    Vo=[]
    for i in range(arrlen):
        uo_time=uo[i,:,:]
        vo_time=vo[i,:,:]
        
        uinterp=interp2d(Xo,Yo,uo_time,kind='linear')
        vinterp=interp2d(Xo,Yo,vo_time,kind='linear')
        #creating winds at buoy locations
        Uo=np.append(Uo,uinterp(Buoy_Xi[i],Buoy_Yi[i]))
        Vo=np.append(Vo,vinterp(Buoy_Xi[i],Buoy_Yi[i]))

    logging.info("Interpolation successfully done.")
    logging.info("Processing complete for ocean currents data.")
    # Vel_transformations
    [Utvec,Uibvec,Uavec,Uovec,Uwvec]=gf.veltransform(Uib,Vib,Ut,Vt,Ua,Va,Uo,Vo,path)
    logging.info("Velocity vectors successfully merged. Model simulation started.")
    #simulations
    #simulation parameters
    seconds= 1.0 #si unit
    minutes= 60.0*seconds
    hours  = 60.0*minutes  
    days   = 24.0*hours
    tstart = 0.0 #arbitrary initial time
    tstep  = 15*minutes
    #tstop  = 61.0*days+23.75*hours
    [ndays,nhrs]=gf.index2time(indexing)
    tstop=days*ndays+nhrs*hours
    times = np.arange(tstart,tstop,tstep)
    [Xis,Yis,Uisvec]=dimodel.simulation(Uavec,Uwvec,Xib[0],Yib[0],Uib[0],Vib[0],times,Cor)

    logging.info("Model Simulations done.")
    gf.save2excel(Tib,Utvec,Uibvec,Uisvec,Uovec,Uwvec,Uavec,path)
    logging.info("Velocity vectors are available in the following loc:"+ path)
    logging.info("Plotting started")
    gp.plticevel(Uisvec,Uibvec,path)
    gp.plticepos(Xib,Yib,Xis,Yis,path)
    logging.info("Plotting completed. Files available in:" +path)
    logging.info("Processing completed for buoy: "+ bname )

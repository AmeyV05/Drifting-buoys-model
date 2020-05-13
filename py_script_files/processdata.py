# processing all nc data to create vectors a for a particular buoy
import pandas as pd
import numpy as np
import scipy.io
from scipy.interpolate import interp2d
from scipy.interpolate import Rbf
import buoydatprocess as bdat
import os
import generalfunc as gf
import logging

def processing(Bnum,indexing,TD,AD,OD,ID,numtaps,sdate):
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
    Xib=Xib[sdate:]
    Yib=Yib[sdate:]
    Uib=Uib[sdate:]
    Vib=Vib[sdate:]
    Tib=Tib[sdate:]
    # #filtered data
    fedge=int(numtaps/2)
    if fedge==0:
        Xibf=Xib
        Yibf=Yib
        logging.info("Buoy data not filtered.")
    else:
        [Xib_res,Xibmod,Xibf]=gf.LPfilter (numtaps,Xib)
        [Yib_res,Yibmod,Yibf]=gf.LPfilter (numtaps,Yib)
        Xibf=np.append(Xib[0:fedge],Xibf)
        Yibf=np.append(Yib[0:fedge],Yibf)

    # obtaining GTSM tidal data in vector form from dictonary TD
    Xt=TD['Xt']; Yt=TD['Yt']
    ut=TD['ut']; vt=TD['vt']
    Tt=TD['Tt']; ssht=TD['ssht']
    sindex=0
    eindex=sindex+indexing
    arrlen=eindex-sindex   
    Tt=gf.num2datetimesecs(2014,3,1,sindex,eindex,Tt)
    logging.info("Interpolating GTSM data to buoy locations....")
    # processing of gtsm data
    #interpolating the u and v velocities to the buoy locations.
    Buoy_Xi=Xib #buoy longitude in hourly basis (x position)
    Buoy_Yi=Yib
    Ut=[];Vt=[];SSHt=[]
    for i in range(arrlen):
        ut_time=ut[i,:];vt_time=vt[i,:];ssht_time=ssht[i,:]
        xk=Buoy_Xi[i];yk=Buoy_Yi[i]
        box_x1=xk-0.45;box_x2=xk+0.45
        box_y1=yk-0.15;box_y2=yk+0.15
        j1=np.where(Xt>=box_x1,Xt,0);j2=np.where(Xt<=box_x2,j1,0)
        i1=np.where(Yt>=box_y1,Yt,0);i2=np.where(Yt<=box_y2,i1,0)
        Xtn=[];Ytn=[]
        ut_ti=[];vt_ti=[];ssht_ti=[]
        for j in range(len(j2)):
            if j2[j]!=0 and i2[j]!=0:
                Xtn=np.append(Xtn,j2[j]);Ytn=np.append(Ytn,i2[j])
                ut_ti=np.append(ut_ti,ut_time[j])
                vt_ti=np.append(vt_ti,vt_time[j])
                ssht_ti=np.append(ssht_ti,ssht_time[j])
        uinterp=Rbf(Xtn,Ytn,ut_ti,kind='linear')
        vinterp=Rbf(Xtn,Ytn,vt_ti,kind='linear')
        sshinterp=Rbf(Xtn,Ytn,ssht_ti,kind='linear')
        #creating tides at buoy locations
        Ut=np.append(Ut,uinterp(Buoy_Xi[i],Buoy_Yi[i]))
        Vt=np.append(Vt,vinterp(Buoy_Xi[i],Buoy_Yi[i]))
        SSHt=np.append(SSHt,sshinterp(Buoy_Xi[i],Buoy_Yi[i]))
    # computing the ssh gradients on the buoy locations. 
    logging.info("Computing ssh gradients for tidal sshs.")
    Pgxt=[];Pgyt=[];dlatvec=dlonvec=dzettvec=[]
    Pgxt=np.append(Pgxt,0);Pgyt=np.append(Pgyt,0)
    for i in range(arrlen-1):
        dzett=SSHt[i+1]-SSHt[i]
        dzettvec=np.append(dzettvec,dzett)
        lat=Buoy_Yi[i]
        dlat=Buoy_Yi[i+1]-Buoy_Yi[i]
        dlatvec=np.append(dlatvec,dlat)
        dlon=Buoy_Xi[i+1]-Buoy_Xi[i]
        dlonvec=np.append(dlonvec,dlon)
        (dy,dx)=gf.latlon2meters(lat,dlat,dlon)
        dy1=abs(dy);dx1=abs(dx)
        if dy1<=1e-25 and dx1>=1e-25 :
            dzdy=0;dzdx=dzett/dx
        elif dy1>=1e-25 and dx1<=1e-25 :
            dzdy=dzett/dy;dzdx=0
        elif dy1<=1e-25 and dx1<=1e-25 :
            dy=0;dzdx=0
        else:
            dzdx=dzett/dx;dzdy=dzett/dy
        Pgxt=np.append(Pgxt,dzdx);Pgyt=np.append(Pgyt,dzdy)
    logging.info("Tidal ssh gradients computation done.")        
    logging.info("Interpolation successfully done. Processing complete for GTSM data.")
   
    sindex=336+int(sdate/4) #converting to 15th march 00
    eindex=sindex+int((indexing/4))
    arrlen=eindex-sindex  
    # obtaining wind data in vector form from dictonary TD
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

    sindex=int(sdate/4)
    eindex=sindex+int(indexing/4)
    arrlen=eindex-sindex
    Xo=OD['Xo']
    Yo=OD['Yo']
    uo=OD['uo']
    vo=OD['vo']
    ssh=OD['ssh']
    To=OD['To']
    To=gf.num2datetimehrs(1950,1,1,sindex,eindex,To)
    uo=uo[sindex:eindex,0,:,:]
    vo=vo[sindex:eindex,0,:,:]
    ssh=ssh[sindex:eindex,0,:,:]
    # processing sea surface heights to create pressure gradients

    #processing of ocean data 
    logging.info("Calculating pressure gradients from sea surface heights.")

    Buoy_Xi=Xib[::4] #buoy longitude in hourly basis (x position)
    Buoy_Yi=Yib[::4]
    Uo=[]
    Vo=[]
    Pgx=[]
    Pgy=[]
    for i in range(arrlen):
    # evaluation for sea surface heights
    # calculating the area of influence for buoy location.
        xk=Buoy_Xi[i];yk=Buoy_Yi[i]
        box_x1=xk-1.5;box_x2=xk+1.5
        box_y1=yk-0.5;box_y2=yk+0.5
        j1=np.array(np.where(Xo<=box_x1));j2=np.array(np.where(Xo<=box_x2))
        i1=np.array(np.where(Yo<=box_y1));i2=np.array(np.where(Yo<=box_y2))
        Xon=Xo[j1[0,-1]:j2[0,-1]];Yon=Yo[i1[0,-1]:i2[0,-1]]
        dzdxmat=np.zeros((len(Yon),len(Xon)))
        for j in range(len(Yon)-1):
            dzdxvec=[]
            dzdxvec=np.append(dzdxvec,0)
            for k in range(len(Xon)-1):
                dzet=ssh[i,j,k+1]-ssh[i,j,k]
                lat=Yo[j]
                dlat=0
                dlon=Xo[k+1]-Xo[k]
                (dy,dx)=gf.latlon2meters(lat,dlat,dlon)
                dzdxvec=np.append(dzdxvec,dzet/dx)
            dzdxmat[j,:]=dzdxvec
        xinterp=interp2d(Xon,Yon,dzdxmat,kind='linear')
        Pgx=np.append(Pgx,xinterp(Buoy_Xi[i],Buoy_Yi[i]))
        dzdymat=np.zeros((len(Yon),len(Xon)))
        for j in range(len(Xon)-1):
            dzdyvec=[]
            dzdyvec=np.append(dzdyvec,0)
            for k in range(len(Yon)-1):
                dzet=ssh[i,k+1,j]-ssh[i,k,j]
                lat=Yo[k]
                dlat=Yo[k+1]-Yo[k]
                dlon=0
                (dy,dx)=gf.latlon2meters(lat,dlat,dlon)
                dzdyvec=np.append(dzdyvec,dzet/dy)
            dzdymat[:,j]=dzdyvec
        yinterp=interp2d(Xon,Yon,dzdymat,kind='linear')
        Pgy=np.append(Pgy,yinterp(Buoy_Xi[i],Buoy_Yi[i]))
    # for ocean currents    
        uo_time=uo[i,:,:]
        vo_time=vo[i,:,:]
        uo_ti=uo_time[i1[0,-1]:i2[0,-1],j1[0,-1]:j2[0,-1]]
        vo_ti=vo_time[i1[0,-1]:i2[0,-1],j1[0,-1]:j2[0,-1]]
        uinterp=interp2d(Xon,Yon,uo_ti,kind='linear')
        vinterp=interp2d(Xon,Yon,vo_ti,kind='linear')
        #creating winds at buoy locations
        Uo=np.append(Uo,uinterp(Buoy_Xi[i],Buoy_Yi[i]))
        Vo=np.append(Vo,vinterp(Buoy_Xi[i],Buoy_Yi[i]))
    logging.info("Done calculating and interpolating the pressure gradients to buoy locations.")
    logging.info("Interpolating ocean currents to buoy locations.")
    logging.info("Interpolation successfully done.")
    logging.info("Processing complete for ocean currents data.")

    # processing for ice thickness data

    sindex=int(sdate/(4*24))
    eindex=sindex+int((indexing/(4*24)))
    arrlen=eindex-sindex 
    # obtaining wind data in vector form from dictonary TD
    Xe=ID['Xe']
    Ye=ID['Ye']
    hi=ID['hi']
    Te=ID['Te']
    Te=gf.num2datetimehrs(1950,1,1,sindex,eindex,Te)
    hi=hi[sindex:eindex,:,:]
    logging.info("Interpolating wind data to buoy locations....")
    # processing of wind data
    #interpolating the u and v velocities to the buoy locations. Note that the ERA5 winds are hourly. 
    Xibn=Xib[12*4:];Yibn=Yib[12*4:]
    Buoy_Xi=Xibn[::4*24] #buoy longitude in hourly basis (x position)
    Buoy_Yi=Yibn[::4*24]
    Hi=[]
    for i in range(arrlen):
        hi_t=hi[i,:,:]
        hi_t=np.where(hi_t==-32767.,0,hi_t)
        hinterp=interp2d(Xe,Ye,hi_t,kind='linear')
        #creating winds at buoy locations
        hii=hinterp(Buoy_Xi[i],Buoy_Yi[i])
        if hii==0:
            Hi=np.append(Hi,Hi[i-1])
        else:
            Hi=np.append(Hi,hii)

    logging.info("Interpolation successfully done for ice thickness")
    logging.info("Processing complete for ice thickness")

    # Vel_transformations
    [Utvec,Uibvec,Uavec,Uovec,Uwvec,Pgvec,Pgtvec,Hivec]=gf.veltransform(Uib,Vib,Ut,Vt,Ua,Va,Uo,Vo,Pgx,Pgy,Pgxt,Pgyt,Hi)
    PD={'Utvec':Utvec,
        'Uibvec':Uibvec,
        'Uavec':Uavec,
        'Uovec':Uovec,
        'Uwvec':Uwvec,
        'Pgvec':Pgvec,
        'Pgtvec':Pgtvec,
        'Hivec':Hivec,
        'Tib':Tib,
        'Xib':Xib,
        'Yib':Yib,
        'Xibf':Xibf,
        'Yibf':Yibf,}
    return(PD)


def mainproc(Bnum,indexing,FD,numtaps,sdate):
    TD=FD['TD']
    AD=FD['AD']
    OD=FD['OD']
    ID=FD['ID']
    prefix="BUOY_"
    bname=prefix+Bnum
    path = "../generated_data/"+bname
    PD=processing(Bnum,indexing,TD,AD,OD,ID,numtaps,sdate)
    Utvec=PD['Utvec']
    Uavec=PD['Uavec']
    Uwvec=PD['Uwvec']
    Uibvec=PD['Uibvec']
    Uovec=PD['Uovec'] 
    Pgvec=PD['Pgvec'];Pgtvec=PD['Pgtvec']
    Hivec=PD['Hivec']
    Xib=PD['Xib']
    Tib=PD['Tib'] 
    Yib=PD['Yib']
    Xibf=PD['Xibf']
    Yibf=PD['Yibf']
    logging.info("Velocity vectors successfully merged.")
    gf.save2excel(Tib,Utvec,Uibvec,Uovec,Uwvec,Pgvec,Pgtvec,Hivec,Uavec,Xib,Yib,Xibf,Yibf,path)
    logging.info("All processed data is available in the following loc:"+ path)
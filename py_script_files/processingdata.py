## This function processes the data from all the models.
# It particularly reads the buoy data and other corresponding model data 
# then interpolates the data to buoy locations and outputs those vectors to an excel file. 
import pandas as pd
import numpy as np
import scipy.io
from scipy.interpolate import interp2d
from scipy.interpolate import Rbf
import generalfunc as gf
import logging
import os

def processall(Bnum,indexing,FD,sdate):
    prefix="BUOY_"
    bname=prefix+Bnum
    ##Buoy data processing
    BD=FD['BD'] #dictionary where buoy_data is saved.
    Xib=BD['Xib'][sdate:]; Yib=BD['Yib'][sdate:]
    Uib=BD['Uib'][sdate:]; Vib=BD['Vib'][sdate:]
    Tib=BD['Tib'][sdate:]
    ### This part was there for creating a buoy dataset which has filtered tides.
    # It is no more required and hence obsolete. But still kept here to have the functionality in future.
    # #filtered data
    # fedge=int(numtaps/2)
    # if fedge==0:
    #     Xibf=Xib
    #     Yibf=Yib
    #     logging.info("Buoy data not filtered.")
    # else:
    #     [Xib_res,Xibmod,Xibf]=gf.LPfilter (numtaps,Xib)
    #     [Yib_res,Yibmod,Yibf]=gf.LPfilter (numtaps,Yib)
    #     Xibf=np.append(Xib[0:fedge],Xibf)
    #     Yibf=np.append(Yib[0:fedge],Yibf)


    ## Processing FES2014
    # pressure gradients from FES2014
    logging.info("Processing FES2014 data.")
    sindex=sdate  #This is important indexing used to create a correct dataset. 
    eindex=sindex+indexing
    arrlen=eindex-sindex
    SD=FD['SD']
    tht=SD['th']; Tft=SD['Tft']
    fdlat=SD['dlola'][1];fdlon=SD['dlola'][0]   
    tht=tht[sindex:eindex,:]
    Tft=gf.num2datetimesecs(2014,3,15,sindex,eindex,Tft*60.)
    Fpgx=[];Fpgy=[]
    for i in range(arrlen):
        dhxt=tht[i,0]-tht[i,1]
        dhyt=tht[i,2]-tht[i,3]
        lat=Yib[i]
        (dy,dx)=gf.latlon2meters(lat,fdlat,fdlon)
        dzdx=dhxt/dx;dzdy=dhyt/dy
        Fpgx=np.append(Fpgx,dzdx)
        Fpgy=np.append(Fpgy,dzdy)

    logging.info("Processing done for FES2014 data. Pressure gradient vectors created.")
    ## Processing GTSM data.
    # obtaining GTSM tidal data in vector form interpolated to buoy locations.
    sindex=sdate-96
    eindex=sindex+indexing
    arrlen=eindex-sindex
    TD=FD['TD']
    Xt=TD['Xt']; Yt=TD['Yt']
    ut=TD['ut']; vt=TD['vt']
    Tt=TD['Tt']; ssht=TD['ssht']
    Tt=gf.num2datetimesecs(2014,3,1,sindex,eindex,Tt)
    ut=ut[sindex:eindex,:];vt=vt[sindex:eindex,:];ssht=ssht[sindex:eindex,:]

    logging.info("Interpolating GTSM data to buoy locations....")
    # processing of gtsm data
    #interpolating the u and v velocities to the buoy locations.
    Buoy_Xi=Xib; Buoy_Yi=Yib#buoy longitude (x position)
    Ut=[];Vt=[];Pgxt=[];Pgyt=[]
    for i in range(arrlen):
        ut_time=ut[i,:];vt_time=vt[i,:];ssht_time=ssht[i,:]
        xk=Buoy_Xi[i];yk=Buoy_Yi[i]
        box_x1=xk-0.45;box_x2=xk+0.45  # box of 0.9 lon length and 0.3 lat height to snap obs in that area for computing interpolation.
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
        #computation of gradient
        # lat=Buoy_Yi[i];dlon=0.02;dlat=0.01 # dlon dlat
        lat=Buoy_Yi[i];dlon=0.2;dlat=0.05
        dzx=sshinterp(Buoy_Xi[i]+dlon,Buoy_Yi[i])-sshinterp(Buoy_Xi[i]-dlon,Buoy_Yi[i])
        dzy=sshinterp(Buoy_Xi[i],Buoy_Yi[i]+dlat)-sshinterp(Buoy_Xi[i],Buoy_Yi[i]-dlat)
        (dy,dx)=gf.latlon2meters(lat,2*dlat,2*dlon)
        dzdx=dzx/dx;dzdy=dzy/dy
        Pgxt=np.append(Pgxt,dzdx);Pgyt=np.append(Pgyt,dzdy)

    logging.info("GTSM Tidal pressure gradients computation done.")        
    logging.info("Interpolation successfully done. Processing complete for GTSM data.")
 
    # Processing of wind data using ERA5  
    sindex=336+int(sdate/4) #converting to 15th march using 336
    eindex=sindex+int((indexing/4))
    arrlen=eindex-sindex
    # obtaining wind data in vector form from dictonary AD
    AD=FD['AD']
    Xa=AD['Xa']; Ya=AD['Ya']
    u10=AD['u10']; v10=AD['v10']
    Ta=AD['Ta']
    Ta=gf.num2datetimesecs(1900,1,1,sindex,eindex,Ta*3600.)
    u10=u10[sindex:eindex,:,:]
    v10=v10[sindex:eindex,:,:]
    logging.info("Interpolating wind data to buoy locations....")

    #interpolating the u and v velocities to the buoy locations. Note that the ERA5 winds are hourly. 
    Buoy_Xi=Xib[::4];     Buoy_Yi=Yib[::4]#buoy longitude in hourly basis (x position)
    Ua=[]; Va=[]
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

    #Processing ocean_currents using CMEMS model of hourly dataset.
    sindex=int(sdate/4)
    eindex=sindex+int(indexing/4)
    arrlen=eindex-sindex
    OD=FD['OD']
    Xo=OD['Xo']; Yo=OD['Yo']
    uo=OD['uo']; vo=OD['vo']
    ssh=OD['ssh']; To=OD['To']
    To=gf.num2datetimesecs(1950,1,1,sindex,eindex,To*3600.)
    uo=uo[sindex:eindex,0,:,:]
    vo=vo[sindex:eindex,0,:,:]
    ssh=ssh[sindex:eindex,0,:,:]
    # processing sea surface heights to calculate pressure gradients
    logging.info("Calculating pressure gradients from sea surface heights.")
    Buoy_Xi=Xib[::4];Buoy_Yi=Yib[::4]  #buoy longitude in hourly basis (x position)    
    Uo=[];Vo=[]
    Pgx=[]; Pgy=[]
    for i in range(arrlen):
    # evaluation for sea surface heights
    # calculating the area of influence for buoy location. similar to tidal calculations,
        xk=Buoy_Xi[i];yk=Buoy_Yi[i]
        box_x1=xk-1.5;box_x2=xk+1.5
        box_y1=yk-0.5;box_y2=yk+0.5
        j1=np.array(np.where(Xo<=box_x1));j2=np.array(np.where(Xo<=box_x2))  #here i n j's are 2D arrays.
        i1=np.array(np.where(Yo<=box_y1));i2=np.array(np.where(Yo<=box_y2))
        Xon=Xo[j1[0,-1]:j2[0,-1]];Yon=Yo[i1[0,-1]:i2[0,-1]]
        sshn=ssh[i,i1[0,-1]:i2[0,-1],j1[0,-1]:j2[0,-1]]
        dzdxmat=np.zeros((len(Yon),len(Xon)))
        for j in range(len(Yon)-1):
            dzdxvec=[]
            dzdxvec=np.append(dzdxvec,0)
            for k in range(len(Xon)-1):
                # dzet=ssh[i,j,k+1]-ssh[i,j,k] #earlier wrong logic
                dzet=sshn[j,k+1]-sshn[j,k]
                lat=Yo[j]
                dlat=0
                # dlon=Xo[k+1]-Xo[k]  ## earlier wrong logic 
                dlon=Xon[k+1]-Xon[k]
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
                # dzet=ssh[i,k+1,j]-ssh[i,k,j]
                dzet=sshn[k+1,j]-sshn[k,j]
                lat=Yo[k]
                # dlat=Yo[k+1]-Yo[k]
                dlat=Yon[k+1]-Yon[k]
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

    # Processing for ice thickness data
    sindex=int(sdate/(4*24))
    eindex=sindex+int((indexing/(4*24)))
    arrlen=eindex-sindex    
    ID=FD['ID']
    Xe=ID['Xe']; Ye=ID['Ye']
    hi=ID['hi']; Te=ID['Te']
    Te=gf.num2datetimesecs(1950,1,1,sindex,eindex,Te*3600.)
    hi=hi[sindex:eindex,:,:]
    logging.info("Interpolating ice thickness to buoy locations....")
    Buoy_Xi=Xib[::4*24] 
    Buoy_Yi=Yib[::4*24]
    Hi=[]
    for i in range(arrlen):
        hi_t=hi[i,:,:]
        hi_t=np.where(hi_t==-32767.,0,hi_t)
        hinterp=interp2d(Xe,Ye,hi_t,kind='linear')
        #calculating ice thickness at buoy locations
        hii=hinterp(Buoy_Xi[i],Buoy_Yi[i])
        if hii==0:
            Hi=np.append(Hi,Hi[i-1])
        else:
            Hi=np.append(Hi,hii)

    logging.info("Interpolation successfully done for ice thickness")
    logging.info("Processing complete for ice thickness")

#  Transformations
## Creating a total vector for ease of representations
    Utvec=np.column_stack((Ut,Vt))
    Pgtvec=np.column_stack((Pgxt,Pgyt))
    Fpgvec=np.column_stack((Fpgx,Fpgy))
    Uibvec=np.column_stack((Uib,Vib))
    Uavec=np.column_stack((Ua,Va))
    Uavec=np.repeat(Uavec,4,axis=0)
    Uovec=np.column_stack((Uo,Vo))
    Uovec=np.repeat(Uovec,4,axis=0)
    Uwvec=Utvec+Uovec
    Pgvec=np.column_stack((Pgx,Pgy))
    Pgvec=np.repeat(Pgvec,4,axis=0)
    Hivec=np.repeat(Hi,4*24,axis=0)

    logging.info("Velocity vectors created successfully.")
# Final Dictionary
    PD={'Utvec':Utvec,'Uibvec':Uibvec,'Uavec':Uavec,'Uovec':Uovec,'Uwvec':Uwvec,
        'Pgvec':Pgvec,'Pgtvec':Pgtvec,'Fpgvec':Fpgvec,'Hivec':Hivec,
        'Xib':Xib,'Yib':Yib,'Tib':Tib,'Tt':Tt,'Ta':Ta,'To':To,'Tft':Tft,'Te':Te}
    return(PD)

def save2excel(PD,path):
    # Create some Pandas dataframes from some data.
    dfd = pd.DataFrame({'Date(GMT)': PD['Tib'][:-1]})
    dfibu = pd.DataFrame({'Uib': PD['Uibvec'][:,0]});dfibv = pd.DataFrame({'Vib': PD['Uibvec'][:,1]})
    dftu = pd.DataFrame({'Ut': PD['Utvec'][:,0]});dftv = pd.DataFrame({'Vt': PD['Utvec'][:,1]})
    dfau = pd.DataFrame({'Ua': PD['Uavec'][:,0]});dfav = pd.DataFrame({'Va': PD['Uavec'][:,1]})
    dfou = pd.DataFrame({'Uo': PD['Uovec'][:,0]});dfov = pd.DataFrame({'Vo': PD['Uovec'][:,1]})
    dfwu = pd.DataFrame({'Uw': PD['Uwvec'][:,0]});dfwv = pd.DataFrame({'Vw': PD['Uwvec'][:,1]})
    dfpgx = pd.DataFrame({'Pgx': PD['Pgvec'][:,0]});dfpgy = pd.DataFrame({'Pgy': PD['Pgvec'][:,1]})
    dfpgxt = pd.DataFrame({'Pgxt': PD['Pgtvec'][:,0]});dfpgyt = pd.DataFrame({'Pgyt': PD['Pgtvec'][:,1]})
    dffpgx = pd.DataFrame({'Fpgx': PD['Fpgvec'][:,0]});dffpgy = pd.DataFrame({'Fpgy': PD['Fpgvec'][:,1]})
    dfibx = pd.DataFrame({'Xib': PD['Xib'][:-1]});dfiby = pd.DataFrame({'Yib': PD['Yib'][:-1]})
    dfta = pd.DataFrame({'Ta': PD['Ta'][:-1]});dftt = pd.DataFrame({'Tt': PD['Tt'][:-1]})
    dfto = pd.DataFrame({'To': PD['To'][:-1]});dftf = pd.DataFrame({'Tft':PD['Tft'][:-1]})
    dfhi=pd.DataFrame({'Hi':PD['Hivec'][:]});dfte = pd.DataFrame({'Te': PD['Te'][:-1]})
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(path+'/Pos_Vel_data.xlsx', engine='xlsxwriter')

    # Position the dataframes in the worksheet.
    dfd.to_excel(writer, startcol=0,startrow=0,index=False)  # Default position, cell A1.
    #buoy
    dfibu.to_excel(writer, startcol=1,startrow=0,index=False);dfibv.to_excel(writer, startcol=2,startrow=0,index=False)
    #GTSM
    dftu.to_excel(writer, startcol=3,startrow=0,index=False);dftv.to_excel(writer, startcol=4,startrow=0,index=False)
    #Ocean
    dfou.to_excel(writer, startcol=5,startrow=0,index=False);dfov.to_excel(writer, startcol=6,startrow=0,index=False)
    #Ocean+Tidal
    dfwu.to_excel(writer, startcol=7,startrow=0,index=False);dfwv.to_excel(writer, startcol=8,startrow=0,index=False)
    #ERA 5 winds
    dfau.to_excel(writer, startcol=9,startrow=0,index=False);dfav.to_excel(writer, startcol=10,startrow=0,index=False)
    # Pressure gradients from ocean currents
    dfpgx.to_excel(writer, startcol=11,startrow=0,index=False);dfpgy.to_excel(writer, startcol=12,startrow=0,index=False)
    # Pressure gradients from tides in GTSM
    dfpgxt.to_excel(writer, startcol=13,startrow=0,index=False);dfpgyt.to_excel(writer, startcol=14,startrow=0,index=False)
    # Pressure gradients from FES2014 
    dffpgx.to_excel(writer,startcol=15,startrow=0,index=False);dffpgy.to_excel(writer,startcol=16,startrow=0,index=False)
    #buoy locs
    dfibx.to_excel(writer, startcol=17,startrow=0,index=False);dfiby.to_excel(writer, startcol=18,startrow=0,index=False)
    # time vecs of wind ocean tidal
    dfta.to_excel(writer, startcol=19,startrow=0,index=False);dftt.to_excel(writer, startcol=20,startrow=0,index=False)
    dfto.to_excel(writer, startcol=21,startrow=0,index=False);dftf.to_excel(writer, startcol=22,startrow=0,index=False)
    # Ice thickness and time for ice thickness
    dfhi.to_excel(writer, startcol=23,startrow=0,index=False);dfte.to_excel(writer, startcol=24,startrow=0,index=False)    
    workbook  = writer.book
    worksheet = writer.sheets['Sheet1']
    merge_format = workbook.add_format({
    'bold': 1,
    'border': 1,
    'align': 'center',
    'valign': 'vcenter',
     })
    workbook.close()

def processingdata(Bnum,indexing,FD,sdate):
    prefix="BUOY_"
    bname=prefix+Bnum
    logging.info("Processing started for buoy: "+ bname)
    # define the name of the directory to be created
    path = "../../generated_data/"+bname
    gf.mkdir_p(path)
    PD=processall(Bnum,indexing,FD,sdate)
    save2excel(PD,path)
    logging.info("All processed data is available in the following loc:"+ path)
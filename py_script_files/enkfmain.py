

#Some issues...WIP
#Enkf_Filter_code
## EnKF simulation for tuning the parameters
import numpy as np
import generalfunc as gf
import genplots as gp
import logging
import icemodel as im
import pandas as pd
import openpyxl
import cartopy.crs as ccrs
import cartopy.feature as cpf 
import matplotlib.pyplot as plt
import settings
rad2deg=180.0/np.pi
deg2rad=np.pi/180.0
ensembleN=10

def enkfsettings():
 s=settings.settings()
 N=s['en']
 logging.info("Number of elements of extended state vector are:"+str(N))
 #mean and cov for state vector addition wk
 cov=0*np.eye(N)
 cov[4,4]=cov[5,5]=1e-8
 cov[0,0]=cov[1,1]=0e01
 mean=np.zeros(N)
 s['cwk']=cov
 s['mwk']=mean
 #mean and cov of obs
 #obs=2 positions
 R=10*np.eye(2)
 meanobs=np.zeros(2)
 s['cvk']=R
 s['mvk']=meanobs
 #H matrix. (obs)
 hrow=np.zeros(N) 
 H=np.row_stack((hrow,hrow))
 H[0,0]=H[1,1]=1
 s['H']=H
 return(s)

def wkcomputations(mean,cov):
 mlen=len(mean)
 wk=np.zeros((mlen))
 for i in range(mlen):
  sigma=np.sqrt(cov[i,i])
  wki=np.random.normal(mean[i],sigma)
  wk[i]=wki
 # sigma=np.sqrt(cov)
 # wk=np.random.normal(mean,sigma)
 return(wk)

def obs2meters(Xib,Yib):
 Xibm=Yibm=[0]
 for i in range(len(Xib)-1):
  Xibc=Xib[i+1]-Xib[i]
  Yibc=Yib[i+1]-Yib[i]
  lat=Yib[i]
  [Yibi,Xibi]=gf.latlon2meters(lat,Yibc,Xibc)
  Xibi=Xibi+Xibm[i];Yibi=Yibi+Yibm[i]
  Xibm=np.append(Xibm,Xibi);Yibm=np.append(Yibm,Yibi) 
 return(Xibm,Yibm)

def Kgain(P_f,H,R):
  # Creation of kalman gain 
 K1=np.matmul(P_f,H.T)
 K2=np.matmul(H,K1)
 K3=K2+R
 K3inv=np.linalg.inv(K3)
 K=np.matmul(K1,K3inv)
 return(K,K3)

#function for the extended enkf model. 
def extmodel(y,wk,consts,s):
  N=s['en']; alpha=1
  My=im.rk2model(y[0:4],consts,s)
  Nk=np.zeros(N);Nk[2]=y[4];Nk[3]=y[5]
  N1k=np.array([alpha*y[4],alpha*y[5]])
  wnk=np.zeros(N-6)
  # wnk[0]=wk[4];wnk[1]=wk[5]
  # for i in range(2,len(wnk)):
  #   wnk[i]=y[i+4]
  Ly=np.concatenate((My,N1k,wnk)) + Nk + wk
  return(Ly)


def simulate(s,Bnum,indexing,Cor):
 [x0,times,PD]=im.initialisation(Bnum,indexing,s)
 Uavec=PD['Uavec']; Utvec=PD['Utvec']; Uovec=PD['Uovec'];Pgvec=PD['Pgvec'];Pgtvec=PD['Pgtvec']
 Hivec=PD['Hivec']
 omega=s['omega']
 tmplier=s['tmplier']
 N=s['en']
 Xis=lon=x0[0]
 Yis=lat=x0[1]
 x=x0 
 x[0]=x[1]=0
 #yo is the extended initial state vector. All the initial values are zero.
 yo=np.zeros(N);yo[0:4]=x 
 #initial values in SI unit assuming the starting latitude as 0m.
 #creation of observation.
 # Xib=PD['Xib'];Yib=PD['Yib']
 Xibf=PD['Xibf'];Yibf=PD['Yibf']
 [Xibm,Yibm]=obs2meters(Xibf,Yibf)
 results=np.zeros((N+1,len(times)))
 zeta_kmat=np.zeros((N,ensembleN))
 for i in range(ensembleN):
  zeta_kmat[:,i]=yo
 zeta_ktot=np.zeros((len(times),N,ensembleN))
 ti=0;y=yo
 # loop
 for t in times[:]:
  zeta_ktot[ti]=zeta_kmat
  #print(t,x)
  if Cor==1:
   f = 2.0*omega*np.sin(lat*deg2rad)
  else:
   f=0
  results[:-1,ti]=y[:];results[-1,ti]=f
  tiv=int(np.floor(ti*tmplier))
  Ua=Uavec[tiv,0];Va=Uavec[tiv,1]
  Ut=Utvec[tiv,0];Vt=Utvec[tiv,1]
  Uo=Uovec[tiv,0];Vo=Uovec[tiv,1]
  Pgx=Pgvec[tiv,0];Pgy=Pgvec[tiv,1]
  Pgxt=Pgtvec[tiv,0];Pgyt=Pgtvec[tiv,1]
  # zeta_kmat[6,:]=Nu[tiv];zeta_kmat[7,:]=Nv[tiv]
  # if Hivec[tiv]>0.2:
  #   h=Hivec[tiv]
  # else:
  #   h=0.2
  h=s['h']
  consts=[f,h,Ua,Va, Ut, Vt, Uo, Vo,Pgx,Pgy,Pgxt,Pgyt]
  obs=np.array([Xibm[tiv],Yibm[tiv]])
  [yn,zeta_hatmat,innovmean,K3]=enkfruns(consts,obs,s,ti,tmplier,zeta_kmat)
  # print(innovmean)
  y  = yn
  # if ((ti*tmplier)%1==0): print (tiv)
  # print(y[6:])
  zeta_kmat=zeta_hatmat
  #coriolis calculator
  dx=yn[0]-results[0,ti]
  dy=yn[1]-results[1,ti]
  #print dx,dy
  (dlon,dlat)=gf.m2latlon(lat,dx,dy)
  lat=lat+dlat
  lon=lon+dlon
  Xis=np.append(Xis,lon)
  Yis=np.append(Yis,lat)
  # print(lon,lat)
  ti+=1
 PD['Xibm']=Xibm; PD['Yibm']=Yibm
 PD['innovmean']=innovmean;PD['K3']=K3

 return(Xis,Yis,results,zeta_ktot,times,PD)

def CFevalsim(s,Bnum,indexing,Cor):
 [x0,times,PD]=im.initialisation(Bnum,indexing,s)
 Uavec=PD['Uavec']; Utvec=PD['Utvec']; Uovec=PD['Uovec']
 Hivec=PD['Hivec'];Pgvec=PD['Pgvec'];Pgtvec=PD['Pgtvec']
 omega=s['omega'];tmplier=s['tmplier']; N=s['en'];h=s['h']
 Xis=lon=x0[0]
 Yis=lat=x0[1]
 x=x0 
 x[0]=x[1]=0
 #yo is the extended initial state vector. All the initial values are zero.
 yo=np.zeros(N);yo[0:4]=x 
 #initial values in SI unit assuming the starting latitude as 0m.
 #creation of observation.
 # Xib=PD['Xib'];Yib=PD['Yib']
 Xibf=PD['Xibf'];Yibf=PD['Yibf']
 [Xibm,Yibm]=obs2meters(Xibf,Yibf)
 resmat=np.zeros((2,len(Xibm)))
 wk1mat=np.zeros((len(times)))
 wk2mat=np.zeros((len(times)))
 zeta_kmat=np.zeros((N,ensembleN))
 Jval=0
 for i in range(ensembleN):
  zeta_kmat[:,i]=yo
 ti=0;y=yo
 # loop
 for t in times[:]:
  #print(t,x)
  if Cor==1:
   f = 2.0*omega*np.sin(lat*deg2rad)
  else:
   f=0
  tiv=int(np.floor(ti*tmplier))
  Ua=Uavec[tiv,0];Va=Uavec[tiv,1]
  Ut=Utvec[tiv,0];Vt=Utvec[tiv,1]
  Uo=Uovec[tiv,0];Vo=Uovec[tiv,1]
  Pgx=Pgvec[tiv,0];Pgy=Pgvec[tiv,1]
  Pgxt=Pgtvec[tiv,0];Pgyt=Pgtvec[tiv,1]
  consts=[f,h,Ua,Va, Ut, Vt, Uo, Vo,Pgx,Pgy,Pgxt,Pgyt]
  obs=np.array([Xibm[tiv],Yibm[tiv]])
  [yn,zeta_hatmat,innovmean,K3]=enkfruns(consts,obs,s,ti,tmplier,zeta_kmat)
  y  = yn
  # if ((ti*tmplier)%1==0): print (tiv)
  zeta_kmat=zeta_hatmat
  [val,res,wk1,wk2]=CFval(s,y,obs,ti)
  if ((ti*tmplier)%1==0):  resmat[:,tiv]=res 
  wk1mat[ti]=wk1;wk2mat[ti]=wk2
  print(resmat[:,tiv])
  Jval=Jval+val
  ti+=1
 # PD['Xibm']=Xibm; PD['Yibm']=Yibm
 return(resmat,wk1mat,wk2mat,Jval,times)

def CFval(s,y,obs,ti):
 H=s['H'];R=s['cvk'];N=s['n'];Q=s['cwk'];h=s['h'];tmplier=s['tmplier']
 Qinv=1./Q[4,4]
 Rinv=R; Rinv[0,0]=1./R[0,0];Rinv[1,1]=1./R[1,1]
 res=np.array([0,0])
 if ((ti*tmplier)%1==0): 
   res=obs-H @ y
   normres=res.T @ Rinv @ res
   wk1=y[-2];wk2=y[-1]
   normwerr=wk1*Qinv*wk1+wk2*Qinv*wk2
   val=normres+normwerr
 else:
   wk1=y[-2];wk2=y[-1]
   normwerr=wk1*Qinv*wk1+wk2*Qinv*wk2
   val=normwerr
 return(val,res,wk1,wk2)

def enkfruns(consts,obs,s,ti,tmplier,zeta_kmat):
 N=s['en'];H=s['H']
 wmean=s['mwk'];Q=s['cwk']
 vmean=s['mvk']; R=s['cvk']
 zeta_hatmat=np.zeros((N,ensembleN))
 innov=np.zeros((2,ensembleN));K3=np.zeros((2,2))
 innovmean=np.mean(innov,axis=1)
 if ((ti*tmplier)%1==0):
  zeta_fkmat=np.zeros((N,ensembleN))
  # for i in range(ensembleN):
  #  yi=zeta_kmat[:,i]
  #  wk=wkcomputations(mean,cov)
  #  zeta_hatmat[:,i] = im.model(xi,consts,s) + wk
  #  zeta_hatmean=np.mean(zeta_hatmat,axis=1)
  #  y=zeta_hatmean
  if (ti==0):
   for i in range(ensembleN):
    yi=zeta_kmat[:,i]
    wk=wkcomputations(wmean,Q)
    zeta_hatmat[:,i] =extmodel(yi,wk,consts,s)
   zeta_hatmean=np.mean(zeta_hatmat,axis=1)
   y=zeta_hatmean
  else:
   for i in range(ensembleN):
    wk=wkcomputations(wmean,Q)
    yi=zeta_kmat[:,i]
    zeta_fkmat[:,i] =extmodel(yi,wk,consts,s)
  #mean and covariance
   yhat_k=np.mean(zeta_fkmat,axis=1).reshape((N,1))
   error_fmat=zeta_fkmat-yhat_k
   error_sum=np.zeros((N,N))
   for i in range(ensembleN):
    error_fimat=np.outer(error_fmat[:,i],error_fmat[:,i].T)
    error_sum=error_sum+error_fimat
   P_f=(1./(ensembleN-1))*error_sum 
   # print("Covariance is:")
   # print(P_f)
   # Kalman gain
   [K,K3]=Kgain(P_f,H,R)
   # print("Kalman gain matrix is:")
   # print(K)
   # print(ti)
   #observation   
   for i in range(ensembleN):
    vk=wkcomputations(vmean,R)
    zk=obs+vk
    zk=zk.reshape((2,1))
    zeta_fk=zeta_fkmat[:,i].reshape((N,1))
    temp=zk-np.matmul(H,zeta_fk)
    zeta_hatk=zeta_fk+np.matmul(K,temp)
    zeta_hatmat[:,i]=zeta_hatk.reshape(N)
    innvtemp=obs.reshape((2,1))- H @ zeta_fk
    innov[:,i]=innvtemp.reshape(2)
   zeta_hatmean=np.mean(zeta_hatmat,axis=1)
   innovmean=np.mean(innov,axis=1)
   print(innovmean)
   innovcov=np.cov(innov)
   # logging.info("Standard deviation of innov from sim is:" )
   # logging.info([np.sqrt(innovcov[0,0]),np.sqrt(innovcov[1,1])])
   # logging.info("Standard deviation of innov from system is:")
   # logging.info([np.sqrt(K3[0,0]),np.sqrt(K3[1,1])])
   y=zeta_hatmean
   # print(x)
   # print(obs)
 else:
  for i in range(ensembleN):
   yi=zeta_kmat[:,i]
   wk=wkcomputations(wmean,Q)
   zeta_hatk=extmodel(yi,wk,consts,s)
   zeta_hatmat[:,i] = zeta_hatk.reshape(N)
  zeta_hatmean=np.mean(zeta_hatmat,axis=1)
  y=zeta_hatmean
 return(y,zeta_hatmat,innovmean,K3)

 
def body(Bnum,indexing,Cor):
# global ensembleN 
# ensembleN=M
 prefix="BUOY_"
 bname=prefix+Bnum
 path = "../generated_data/"+bname
 logging.info("Model simulation started.")
 s=enkfsettings()
 [Xis,Yis,results,zeta_ktot,times,PD]=simulate(s,Bnum,indexing,Cor)
 logging.info("Model Simulations done.")
 xis=results[0,:]; yis=results[1,:]
 Uis=results[2,:]; Vis=results[3,:]
 Nu=results[4,:];Nv=results[5,:]
 Uisvec=np.column_stack((Uis,Vis))    
 Uibvec=PD['Uibvec']
 # Xib=PD['Xib'];Yib=PD['Yib']
 Xibf=PD['Xibf'];Yibf=PD['Yibf']
 tmplierinv=int(1/s['tmplier'])
 logging.info("Plotting started")   
 gp.plticevel(Uisvec,Uibvec,tmplierinv,path)
 gp.plticepos(Xibf,Yibf,Xis,Yis,path)
 gp.pltalpha(Nu,path,'Nu')
 gp.pltalpha(Nv,path,'Nv')
 logging.info("Plotting completed. Files available in:" +path)
 #pltensembspread(zeta_ktot,times,PD,path)
 # simdat2excel(path,results,Xis,Yis,tmplierinv)
 logging.info("Processing completed for buoy: "+ bname )
 # rmspos=poserror(Xib,Yib,Xis,Yis,tmplierinv)
 # logging.info("RMS error in position for "+bname+" is:" + str(rmspos))

def CFbody(Bnum,indexing,Cor):
# global ensembleN 
# ensembleN=M
 prefix="BUOY_"
 bname=prefix+Bnum
 path = "../generated_data/"+bname
 logging.info("Model simulation started.")
 s=enkfsettings()
 [resmat,wk1mat,wk2mat,Jval,times]=CFevalsim(s,Bnum,indexing,Cor) 
 proEnkfCFdata(resmat,wk1mat,wk2mat,Jval,times,path)


def proEnkfCFdata(resmat,wk1mat,wk2mat,Jval,times,path):
 logging.info("Model simulation completed")
 logging.info("Value of Cost function is:"+ str(Jval/2.))
 logging.info("Plotting started.")
 gp.pltalpha(resmat[0,:],path,'res_x')
 gp.pltalpha(resmat[1,:],path,'res_y')
 gp.pltalpha(wk1mat,path,'wk_x')
 gp.pltalpha(wk2mat,path,'wk_y')
 logging.info("Plotting completed. Files available in:" +path)
 # logging.info("Processing completed for buoy: "+ bname )
 return(Jval)




def poserror(Xib,Yib,Xis,Yis,tmplierinv):
 Xis=Xis[::tmplierinv]
 Yis=Yis[::tmplierinv]
 dX=np.sqrt(((Xib - Xis[:-1]) ** 2).mean())
 dY=np.sqrt(((Yib - Yis[:-1]) ** 2).mean())
 rmspos=[dX,dY]
 return(rmspos)

def simdat2excel(path,results,Xis,Yis,PD,tmplierinv):
 tinv=tmplierinv
 Uis=results[2,:]; Vis=results[3,:]
 Uisvec=np.column_stack((Uis,Vis)) 
 Nu=results[4,:];Nv=results[5,:]
 Uis=Uis[::tinv];Vis=Vis[::tinv]
 Xis=Xis[::tinv];Yis=Yis[::tinv]
 Nu=Nu[::tinv];Nv=Nv[::tinv]
 Xibm=PD['Xibm'];Yibm=PD['Yibm']
 Xism=results[0,:];Yism=results[1,:]
 # workbook1 = openpyxl.load_workbook(path+'/Pos_Vel_data.xlsx')
 # writer = pd.ExcelWriter(path+'/Pos_Vel_data.xlsx', engine='openpyxl') 
 # writer.book = workbook1
  # Position the dataframes in the worksheet.
 writer = pd.ExcelWriter(path+'/Sim_data.xlsx', engine='xlsxwriter')
 dfisx = pd.DataFrame({'Xis': Xis})
 dfisy = pd.DataFrame({'Yis': Yis})
 dfusx = pd.DataFrame({'Uis': Uis})
 dfvsy = pd.DataFrame({'Vis': Vis})
 dfnu = pd.DataFrame({'Nu': Nu})  
 dfnv = pd.DataFrame({'Nv': Nv}) 
 dfibmx=pd.DataFrame({'Xibm':Xibm})
 dfibmy=pd.DataFrame({'Yibm':Yibm}) 
 dfismx=pd.DataFrame({'Xism':Xism[::tinv]})
 dfismy=pd.DataFrame({'Yism':Yism[::tinv]})
 dfisx.to_excel(writer, startcol=0,startrow=0,index=False)
 dfisy.to_excel(writer, startcol=1,startrow=0,index=False)
 dfusx.to_excel(writer, startcol=2,startrow=0,index=False,sheet_name='Sheet1')
 dfvsy.to_excel(writer, startcol=3,startrow=0,index=False,sheet_name='Sheet1')
 dfnu.to_excel(writer, startcol=4,startrow=0,index=False,sheet_name='Sheet1')
 dfnv.to_excel(writer, startcol=5,startrow=0,index=False,sheet_name='Sheet1')
 dfibmx.to_excel(writer, startcol=6,startrow=0,index=False,sheet_name='Sheet1')
 dfibmy.to_excel(writer, startcol=7,startrow=0,index=False,sheet_name='Sheet1')
 dfismx.to_excel(writer, startcol=8,startrow=0,index=False,sheet_name='Sheet1')
 dfismy.to_excel(writer, startcol=9,startrow=0,index=False,sheet_name='Sheet1')
 # writer.save()
 # print('Pos_Vel_data file updated.')
 # workbook1.close()
 workbook  = writer.book
 worksheet = writer.sheets['Sheet1']
 merge_format = workbook.add_format({
    'bold': 1,
    'border': 1,
    'align': 'center',
    'valign': 'vcenter',
     })
 workbook.close()


def excelinnov(PD):
  K3=PD['K3']
  innovmean=PD['innovmean']
  K3vec=[np.sqrt(K3)]


def pltensembspread(zeta_ktot,times,PD,path):
 Xib=PD['Xib'];Yib=PD['Yib']
 fig=plt.figure(figsize=(12, 12), frameon=True)
 ax=plt.axes(projection=ccrs.LambertAzimuthalEqualArea(central_longitude=25.0,central_latitude=77.0)) 
 ax.set_extent([15,33,74,81]) 
 plt.plot(Xib,Yib,color='red',transform=ccrs.PlateCarree(),label='buoy drift')
 for i in range(ensembleN):
  ti=1
  Xzi=[];Yzi=[]
  lon=Xzi=Xib[0];lat=Yzi=Yib[0]
  for t in times[:-2]:
   dx=zeta_ktot[ti,0,i]-zeta_ktot[ti-1,0,i]
   dy=zeta_ktot[ti,1,i]-zeta_ktot[ti-1,1,i]
   (dlon,dlat)=gf.m2latlon(lat,dx,dy)
   lat=lat+dlat
   lon=lon+dlon
   Xzi=np.append(Xzi,lon)
   Yzi=np.append(Yzi,lat)
   ti=ti+1
  print(Xzi[0],Yzi[0],Xzi[-1],Yzi[-1])
  plt.plot(Xzi[:],Yzi[:],color='yellow',alpha=0.7,transform=ccrs.PlateCarree())
 # plt.legend()
 service='https://www.gebco.net/data_and_products/gebco_web_services/web_map_service/mapserv?'
 ax.add_wms(service,layers=['GEBCO_LATEST'],wms_kwargs={'width':900*2,'height':600*2}) 
 plt.savefig(path+'/ensemspread.jpg',dpi=500)
 plt.close(fig) 

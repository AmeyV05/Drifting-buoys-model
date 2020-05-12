# optimization_Cws
import numpy as np 
import enkfmain as enkf 
import settings 
import generalfunc as gf
import genplots as gp
import logging
import icemodel as im

s=enkf.enkfsettings() 
def CFval(s,times,PD,results):
 H=s['H'];R=s['cvk'];N=s['n'];Q=s['cwk']
 h=s['h']
 # Uavec=PD['Uavec']; Utvec=PD['Utvec']; Uovec=PD['Uovec'];Pgvec=PD['Pgvec']
 tmplier=s['tmplier']
 logging.info("Value of ice-drag coefficient for computation of cost function is: "+ str(s['iCw']))
 # Qinv=Q; Qinv[4,4]=1./Q[4,4];Qinv[5,5]=1./Q[5,5]
 Qinv=1./Q[4,4]
 Rinv=R; Rinv[0,0]=1./R[0,0];Rinv[1,1]=1./R[1,1]
 Zkmat=np.row_stack((PD['Xibm'],PD['Yibm']))
 Ykmat=results[:-1,:]
 # ti=0
 val=0; simend=len(times)-2*int(1/s['tmplier']+1)
 simstrt=2*int(1/s['tmplier'])
 ti=simstrt
 for t in times[simstrt:simend]:
  tiv=int(np.floor(ti*tmplier))
  # f=results[-1,ti]
  # Ua=Uavec[tiv,0];Va=Uavec[tiv,1]
  # Ut=Utvec[tiv,0];Vt=Utvec[tiv,1]
  # Uo=Uovec[tiv,0];Vo=Uovec[tiv,1]
  # Pgx=Pgvec[tiv,0];Pgy=Pgvec[tiv,1]
  # consts=[f,h,Ua,Va, Ut, Vt, Uo, Vo,Pgx,Pgy]
  Yk=Ykmat[:,ti] 
  # ;Yk1=Ykmat[:,ti+1]
  if ((ti*tmplier)%1==0): 
   Zk=Zkmat[:,tiv]
   res=Zk-H @ Yk
   normres=res.T @ Rinv @ res
   # mYk1=im.expeumodel(Yk,consts,s)
   # err=Yk1-mYk1
   # normerr=err.T @ Qinv @ err
   wk1=Yk[-2];wk2=Yk[-1]
   normwerr=wk1*Qinv*wk1+wk2*Qinv*wk2
   print(normwerr)
   val=val+normres+normwerr
  else:
   # mYk1=im.expeumodel(Yk,consts,s)
   # err=Yk1-mYk1
   # normerr=err.T @ Qinv @ err
   wk1=Yk[-2];wk2=Yk[-1]
   normwerr=wk1*Qinv*wk1+wk2*Qinv*wk2
   print(normwerr)
   val=val+normwerr
  ti+=1
 # tiv=int(np.floor(ti*tmplier))
 # if ((ti*tmplier)%1==0):
 #   Zk=Zkmat[:,tiv]; Yk=Ykmat[:,ti]
 #   res=Zk-H @ Yk
 #   normres=res.T @ Rinv @ res
 # val=val+normres
 Jval=0.5*val
 return(Jval)
  
def proEnkfdata(s,Xis,Yis,results,zeta_ktot,times,PD,path,i):
 logging.info("Model Simulations done for iteration"+ str(i+1))
 xis=results[0,:]; yis=results[1,:]
 Uis=results[2,:]; Vis=results[3,:]
 Nu=results[4,:];Nv=results[5,:]
 Uisvec=np.column_stack((Uis,Vis))    
 Uibvec=PD['Uibvec']
 Xibf=PD['Xibf'];Yibf=PD['Yibf']
 tmplierinv=int(1/s['tmplier'])
 logging.info("Plotting started")   
 gp.plticevel(Uisvec,Uibvec,tmplierinv,path)
 gp.plticepos(Xibf,Yibf,Xis,Yis,path)
 gp.pltalpha(Nu,path,'Nu')
 gp.pltalpha(Nv,path,'Nv')
 logging.info("Plotting completed. Files available in:" +path)
 enkf.simdat2excel(path,results,Xis,Yis,PD,tmplierinv)
 logging.info("Successfully created Excel file for the simulation data.")


    
 
### Ice model 

import numpy as np 
import pandas as pd 
import generalfunc as gf

rad2deg=180.0/np.pi
deg2rad=np.pi/180.0

def readveldata(path):
 D=pd.read_excel(path+'/Pos_Vel_data.xlsx')
 Uib=np.array(D['Uib'])
 Vib=np.array(D['Vib'])
 Ut=np.array(D['Ut'])
 Vt=np.array(D['Vt'])
 Ua=np.array(D['Ua'])
 Va=np.array(D['Va'])
 Uo=np.array(D['Uo'])
 Vo=np.array(D['Vo'])
 Pgx=np.array(D['Pgx']);Pgxt=np.array(D['Pgxt'])
 Pgy=np.array(D['Pgy']);Pgyt=np.array(D['Pgyt'])
 Tib=np.array(D['Date(GMT)'])
 Xib=np.array(D['Xib'])
 Yib=np.array(D['Yib'])
 Xibf=np.array(D['Xibf'])
 Yibf=np.array(D['Yibf'])
 Uibvec=np.column_stack((Uib,Vib))
 Utvec=np.column_stack((Ut,Vt))
 Uavec=np.column_stack((Ua,Va))
 Uovec=np.column_stack((Uo,Vo))
 Pgvec=np.column_stack((Pgx,Pgy))
 Pgtvec=np.column_stack((Pgxt,Pgyt))
 Hivec=np.array(D['Hi'])
 Uwvec=Uovec+Utvec
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
    'Yibf':Yibf}
 return(PD)

def initialisation(Bnum,indexing,s):
 prefix="BUOY_"
 bname=prefix+Bnum
 path = "../generated_data/"+bname
 PD=readveldata(path)
 # s=settings()
 dt=s['dt']
 days=s['days']
 hours=s['hours']
 N=s['n']
 tstart = 0.0 
 [ndays,nhrs]=gf.index2time(indexing)
 tstop=days*ndays+nhrs*hours
 times = np.arange(tstart,tstop,dt)
 Xib=PD['Xib']
 Yib=PD['Yib']
 Uibvec=PD['Uibvec']
 x0=[Xib[0],Yib[0],Uibvec[0,0],Uibvec[0,1]]
 return(x0,times,PD)



def modelFx(x,consts,s):
 dt=s['dt']
 N=s['n']
 g=s['g']
 thetaa=s['thetaa']
 thetaw=s['thetaw']
 rho_water=s['rho_water']
 rho_air=s['rho_air']
 rho_ice=s['rho_ice']
 [f, h, Ua, Va, Ut, Vt, Uo, Vo,Pgx,Pgy,Pgxt,Pgyt] =consts
 [xc,yc,u,v] = x
 Cwi=s['iCw'];Cai=s['iCa']
 #transformation for rotation of Ua
 Ua=Ua*np.cos(thetaa)-Va*np.sin(thetaa)
 Va=Ua*np.sin(thetaa)+Va*np.cos(thetaa)
 #transformation for rotation of Uw
 Uw=Uo+Ut; Vw=Vo+Vt
 Uw=(Uw-u)*np.cos(thetaw)-(Vw-v)*np.sin(thetaw)
 Vw=(Uw-u)*np.sin(thetaw)+(Vw-v)*np.cos(thetaw)
 # Ut=(Ut)*np.cos(thetaw)-(Vt)*np.sin(thetaw)
 # Vt=(Ut)*np.sin(thetaw)+(Vt)*np.cos(thetaw)
 # Uw=Ut+Uo-u; Vw=Vt+Vo-v
 # # new trial with Ro doesn't work.  
 # Ro=s['Ro']/f
 # u=u*(1/(1+Ro*np.sin(thetaw)))-v*(1/(1+Ro*np.cos(thetaw)))
 # v=u*(1/(1+Ro*np.cos(thetaw)))+v*(1/(1+Ro*np.sin(thetaw)))
 # Ut=Ut*(1/(1+Ro*np.sin(thetaw)))-Vt*(1/(1+Ro*np.cos(thetaw)))
 # Vt=Ut*(1/(1+Ro*np.cos(thetaw)))+Vt*(1/(1+Ro*np.sin(thetaw)))
 # Uw=Uw+Uo;Vw=Vw+Vo
 u_mag = np.sqrt((Uw)**2+(Vw)**2)
 Ua_mag=np.sqrt((Ua)**2+(Va)**2)
 Fx=np.zeros(N)
 Fx[0]=u;Fx[1]=v
 Fx[2]=((Cwi*rho_water*Uw*u_mag + 
         Cai*rho_air*Ua*Ua_mag)/(rho_ice*h)+((f*v)))-g*Pgx # -g*Pgxt
 Fx[3]=((Cwi*rho_water*Vw*u_mag+
         Cai*rho_air*Va*Ua_mag)/(rho_ice*h)-((f*u)))-g*Pgy #-g*Pgyt
 #linear shear stress model
 # Fx[2]=((Cwi*rho_water*Uw +
 #         Cai*rho_air*Ua)/(rho_ice*h)+((f*v))) -g*Pgx
 # Fx[3]=((Cwi*rho_water*Vw+
 #         Cai*rho_air*Va)/(rho_ice*h)-((f*u))) -g*Pgy
 # Fx[4] = 0; Fx[5] = 0
 return (Fx)

def expeumodel(x,consts,s):
 dt=s['dt']
 N=s['n']
 Fx=modelFx(x,consts,s) 
 #Nk=np.array(([0,0,x[6],x[7],0,0,0,0])) when Nk forcing is added in EnKF
 # Nk=np.array(([0,0,x[4],x[5],0,0]))
 Mx=x+dt*Fx
 # Mx=x+dt*Fx
 return Mx

def rk2model(x,consts,s):
 dt=s['dt']
 N=s['n']
 Fx2=modelFx(x,consts,s)
 Mx2=x+(dt/2)*Fx2
 Fx=modelFx(Mx2,consts,s)
 # Nk=np.zeros(N);Nk[2]=x[4];Nk[3]=x[5]
 Mx=x+dt*Fx
 # +Nk
 # Mx=x+dt*Fx
 return Mx

def main():
	#doesn't work anymore
 Bnum='02'
 indexing=3152
 s=settings()
 [x0,times,PD]=initialisation(Bnum,indexing,s)
 Uavec=PD['Uavec']
 Utvec=PD['Utvec']
 Uovec=PD['Uovec']
 ti=1
 f=0.0001
 h=0.1
 Ua=Uavec[ti,0];Va=Uavec[ti,1]
 Ut=Utvec[ti,0];Vt=Utvec[ti,1]
 Uo=Uavec[ti,0];Vo=Uovec[ti,1]
 consts=[f,h, Ua, Va, Ut, Vt, Uo, Vo]
 Mnew=model(x0,consts,s)
 print(Mnew)
 
if __name__ == '__main__':
	main()




### the simple ice model script
#This model has 2 parts:
# initialisation, model def. and simulation or running.
import numpy as np 
import pandas as pd 
import generalfunc as gf
import logging
import settings
rad2deg=180.0/np.pi
deg2rad=np.pi/180.0


## Reading the pos velocity data from Pos_Vel_data file
def readposveldata(path):
    D=pd.read_excel(path+'/Pos_Vel_data.xlsx')
    Uib=np.array(D['Uib']);Vib=np.array(D['Vib'])
    Ut=np.array(D['Ut']); Vt=np.array(D['Vt'])
    Ua=np.array(D['Ua']);Va=np.array(D['Va'])
    Uo=np.array(D['Uo']);Vo=np.array(D['Vo'])
    Pgx=np.array(D['Pgx']);Pgxt=np.array(D['Pgxt']);Fpgx=np.array(D['Fpgx'])
    Pgy=np.array(D['Pgy']);Pgyt=np.array(D['Pgyt']);Fpgy=np.array(D['Fpgy'])
    Tib=np.array(D['Date(GMT)']);Hivec=np.array(D['Hi'])
    Xib=np.array(D['Xib']);Yib=np.array(D['Yib'])
    Uibvec=np.column_stack((Uib,Vib))
    Utvec=np.column_stack((Ut,Vt))
    Uavec=np.column_stack((Ua,Va))
    Uovec=np.column_stack((Uo,Vo))
    Pgvec=np.column_stack((Pgx,Pgy))
    Pgtvec=np.column_stack((Pgxt,Pgyt))
    Fpgvec=np.column_stack((Fpgx,Fpgy))
    Uwvec=Uovec+Utvec
    PD={'Utvec':Utvec,'Uibvec':Uibvec,'Uavec':Uavec,'Uovec':Uovec,
        'Uwvec':Uwvec,'Pgvec':Pgvec,'Pgtvec':Pgtvec,'Fpgvec':Fpgvec,
        'Hivec':Hivec,'Tib':Tib,'Xib':Xib,'Yib':Yib}
    return(PD)

## Model initalisation.
#PD is the dictionary which stores all the vector data from pos_vel_data excel file.
def initialisation(s,Bnum,indexing,PD):
    prefix="BUOY_"
    bname=prefix+Bnum
    path = "../../generated_data/"+bname
    dt=s['dt'];days=s['days']
    hours=s['hours'];N=s['n']
    tstart = 0.0 
    [ndays,nhrs]=gf.index2time(indexing)
    tstop=days*ndays+nhrs*hours
    times = np.arange(tstart,tstop,dt)
    x0=[PD['Xib'][0],PD['Yib'][0],PD['Uibvec'][0,0],PD['Uibvec'][0,1]]
    return(x0,times)

## Model is given by x_k+1=Mx_k
# We have two time integration schemes so two models
# explicit euler and another is runge kutta 2.

def expeumodel(x,consts,s):
    dt=s['dt']; N=s['n']
    Fx=modelFx(x,consts,s) 
    Mx=x+dt*Fx
    return Mx

def rk2model(x,consts,s):
    dt=s['dt'];N=s['n']
    Fx2=modelFx(x,consts,s)
    Mx2=x+(dt/2)*Fx2
    Fx=modelFx(Mx2,consts,s)
    Mx=x+dt*Fx
    return Mx

# Model function def. M_x=x+dt*Fx
def modelFx(x,consts,s):
    dt=s['dt']; N=s['n']
    g=s['g'];thetaa=s['thetaa']
    thetaw=s['thetaw'];rho_water=s['rho_water']
    rho_air=s['rho_air']; rho_ice=s['rho_ice']
    [f, h, Ua, Va, Ut, Vt, Uo, Vo,Pgx,Pgy,Pgxt,Pgyt] =consts
    [xc,yc,u,v] = x
    Cwi=s['iCw'];Cai=s['iCa']
    #transformation for rotation of Ua Uaê^(-itheta)
    # Ua=Va=0
    Ua=Ua*np.cos(thetaa)-Va*np.sin(thetaa)
    Va=Ua*np.sin(thetaa)+Va*np.cos(thetaa)
    #transformation for rotation of Uw
    Uw=Ut+Uo; Vw=Vt+Vo
    Uw=(Uw-u)*np.cos(thetaw)-(Vw-v)*np.sin(thetaw)
    Vw=(Uw-u)*np.sin(thetaw)+(Vw-v)*np.cos(thetaw)
    Uw=0 if Ut==0 and Uo==0 else Uw
    Vw=0 if Vt==0 and Vo==0 else Vw
    u_mag = np.sqrt((Uw)**2+(Vw)**2)
    Ua_mag=np.sqrt((Ua)**2+(Va)**2)
    Fx=np.zeros(N)
    Fx[0]=u;Fx[1]=v
    Fx[2]=((Cwi*rho_water*Uw*u_mag + 
         Cai*rho_air*Ua*Ua_mag)/(rho_ice*h))+f*(v)-g*( Pgxt+Pgx)
    Fx[3]=((Cwi*rho_water*Vw*u_mag+
         Cai*rho_air*Va*Ua_mag)/(rho_ice*h))-f*(u)-g*( Pgyt+Pgy)  
    ## linear shear stress model to check and verify
    # Fx[2]=((Cwi*rho_water*Uw +
    #         Cai*rho_air*Ua)/(rho_ice*h)+((f*v)))-g*( Pgxt+Pgx)
    # Fx[3]=((Cwi*rho_water*Vw+
    #         Cai*rho_air*Va)/(rho_ice*h)-((f*u))) -g*( Pgyt+Pgy) 
    return (Fx)


## model simulation
def simulate(s,Bnum,indexing,forcevec,PD):
    [x0,times]=initialisation(s,Bnum,indexing,PD)
    Uavec=PD['Uavec']; Utvec=PD['Utvec']; Uovec=PD['Uovec'];Pgvec=PD['Pgvec'];Pgtvec=PD['Pgtvec']
    ho=s['h'];trate=s['trate'];Yib=PD['Yib'];Fpgvec=PD['Fpgvec'];omega=s['omega'];tmplier=s['tmplier']
    logging.info("Time step is:" +str(s['dt'])+' seconds')
    logging.info("Ice Air drag coefficient is:" +str(s['iCa']))
    logging.info("Ice Water drag coefficient is:" +str(s['iCw']))
    lon=x0[0];lat=x0[1]
    Xis=lon;Yis=lat
    x=x0; x[0]=x[1]=0 #initial values in SI unit assuming the starting latitude as 0m.
    results=np.zeros((len(x),len(times)));hvec=[]
    ti=0
    # loop run to all but last element.
    for t in times[:-1]:
        results[:,ti]=x[:]
        #constants evaluation.
        tiv=int(np.floor(ti*tmplier))
        f = 2.0*omega*np.sin(Yib[tiv]*deg2rad)
        Ua=Uavec[tiv,0];Va=Uavec[tiv,1]
        Ut=Utvec[tiv,0];Vt=Utvec[tiv,1]
        Uo=Uovec[tiv,0];Vo=Uovec[tiv,1]
        Pgx=Pgvec[tiv,0];Pgy=Pgvec[tiv,1]
        Pgxt=Pgtvec[tiv,0];Pgyt=Pgtvec[tiv,1]
        Fpgx=Fpgvec[tiv,0];Fpgy=Fpgvec[tiv,1]
        #  ice thickness with thinning rate (trate); trate=0 for constant thickness 
        if ((ti*tmplier)%1==0):
            h=gf.thinrate(ho,trate)   
        hvec=np.append(hvec,h*forcevec[1]);ho=h
        consts=[f,h,Ua, Va, Ut, Vt, Uo, Vo,Pgx,Pgy,Pgxt,Pgyt]
        consts[10:]=[Pgxt,Pgyt] if s['tidepg']=='GTSM' else [Fpgx,Fpgy]
        consts=np.multiply(forcevec,consts) #To get ice thickness in it.
        if (s['mod']=='ExplicitEuler'):
            xn=expeumodel(x,consts,s)
        else:
            xn=rk2model(x,consts,s) 
        x=xn
        # lonlat caluculations.
        dx=xn[0]-results[0,ti]
        dy=xn[1]-results[1,ti]
        (dlon,dlat)=gf.m2latlon(lat,dx,dy)
        lat=lat+dlat
        lon=lon+dlon
        Xis=np.append(Xis,lon)
        Yis=np.append(Yis,lat)
        ti+=1

    ## storing all the simulated data as updation into PD at buoy time instants.
    tmplierinv=int(1/tmplier)
    Uis=results[2,:];Vis=results[3,:]
    xis=results[0,:];yis=results[1,:]
    Uisvec=np.column_stack((Uis[::tmplierinv],Vis[::tmplierinv]))
    PD['xis']=xis[::tmplierinv];PD['yis']=yis[::tmplierinv]
    PD['Xis']=Xis[::tmplierinv];PD['Yis']=Yis[::tmplierinv]
    PD['Uisvec']=Uisvec; PD['hvec']=hvec[::tmplierinv]
    return(PD)


def main():
    Bnum='16'
    indexing=5952-96*16
    forcevec=[1,1,1,1,1,1,1,1,1,1,1,1]
    prefix="BUOY_"
    bname=prefix+Bnum
    path = "../../generated_data/"+bname
    PD=readposveldata(path)   
    s=settings.settings()
    h=s['h'];trate=s['trate']
    (forcenam,folname)=gf.forcedetail(forcevec,trate,h)
    #creation of the folder for storing the simulated data.
    path=path+'/'+folname
    gf.mkdir_p(path)
    logging.info("Model simulation started.")
    PD=simulate(s,Bnum,indexing,forcevec,PD)
    logging.info("Simulation done")
 
if __name__ == '__main__':
	main()



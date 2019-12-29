#ice_drift_model

import numpy as np 
import generalfunc as gf 


#parameter definitions.

seconds= 1.0 #si unit
minutes= 60.0*seconds
hours  = 60.0*minutes  
days   = 24.0*hours
rad2deg= 180.0/np.pi
deg2rad= np.pi/180.0
thetaa = 0*deg2rad#angle made by geostrophic wind with ice velocity.
thetaw = 0*deg2rad
omega  = (2.0*np.pi)/(24.0*hours)
r_ice_water =5.5e-03
r_air_ice =1.5e-03
rho_ice = 920.0
rho_air = 1.2
rho_water = 1024.0
h =1.5 #ice thickness
tstep  = 15*minutes

## Creation of model for simulation of wind and water velocities and then,

#modification based on geostrophic velocity ( exp(i theta)=sin(theta)+cos(theta)) ( exp-V ;expU)
def model(x,ti,f,Uwvec,Uavec):
 xc = x[0]
 yc = x[1]
 u  = x[2]
 v  = x[3]
 dx=np.zeros(4)
 dx[0] = d_xc=u
 dx[1] = d_yc=v
 u_mag = np.sqrt((u-Uwvec[ti,0])**2+(v-Uwvec[ti,1])**2)
 #transformation for rotation of Ua
 Ua=Uavec[ti,0]*np.cos(thetaa)-Uavec[ti,1]*np.sin(thetaa)
 Va=Uavec[ti,0]*np.sin(thetaa)+Uavec[ti,1]*np.cos(thetaa)
 Uw=(Uwvec[ti,0]-u)*np.cos(thetaw)-(Uwvec[ti,1]-v)*np.sin(thetaw)
 Vw=(Uwvec[ti,0]-u)*np.sin(thetaw)+(Uwvec[ti,1]-v)*np.cos(thetaw)
 dx[2] = d_u  = (r_ice_water*rho_water*Uw*u_mag + 
        r_air_ice*rho_air*Ua*abs(Ua))/(rho_ice*h)-((f*Vw)) 
 dx[3] = d_v  = (r_ice_water*rho_water*Vw*u_mag+
       r_air_ice*rho_air*Va*abs(Va))/(rho_ice*h)+((f*Uw))
 return dx


#run simulation
def simulation(Uavec,Uwvec,slon,slat,su,sv,times,Cor):
 #initialisation
 # starting latitude and longitude
 # computation of latitudes and longitudes from X and Y positions
 Xis=slon
 Yis=slat
 lon=slon
 lat=slat
 x=[0, 0, su, sv] #initial values in SI unit
 results=np.zeros((len(x),len(times)))
 ti=0
 # loop
 for t in times:
  #print(t,x)
  if Cor==1:
   f = 2.0*omega*np.sin(lat*deg2rad)
  else:
   f=0
  results[:,ti]=x[:]
  dx=model(x,ti,f,Uwvec,Uavec)
  xn = x + tstep*dx 
  x  = xn
  #print(x)
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

 #simulated ice positions and velocity
 xis=results[0,:] #in meters
 yis=results[1,:]
 Uis=results[2,:]
 Vis=results[3,:]
 Uisvec=np.column_stack((Uis,Vis))
 return(Xis,Yis,Uisvec)
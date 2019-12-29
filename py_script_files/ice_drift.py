import cartopy.crs as ccrs
import cartopy.feature as cpf 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import scipy.io
import netCDF4 as nc4
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp2d
np.set_printoptions(threshold=np.inf)

#obtaining velocity vectors for water, air and ice from csv files
UDi=pd.read_csv('Uivec.csv')
Uivecstr=np.array(UDi)
Uivec = [float(x) for x in Uivecstr]

UDw=pd.read_csv('Uwvec.csv')
Uwvecstr=np.array(UDw)
Uwvec = [float(x) for x in Uwvecstr]
UDa=pd.read_csv('Uavec.csv')
Uavecstr=np.array(UDa)
Uavec = [float(x) for x in Uavecstr]
#parameter definitions.

seconds= 1.0 #si unit
minutes= 60.0*seconds
hours  = 60.0*minutes  
days   = 24.0*hours
rad2deg= 180.0/np.pi
deg2rad= np.pi/180.0

#simulation parameters
tstart = 0.0 #arbitrary initial time
tstep  = 15*minutes
tstop  = 61.0*days+23.75*hours
lat    = 76.0 
omega  = (2.0*np.pi)/(24.0*hours)
f      = 2.0*omega*np.sin(lat*deg2rad)
r_ice_water = 5.5e-04
r_air_ice = 1.2e-05
rho_ice = 920.0
rho_air = 1.2
rho_water = 1024.0
h = 0.1 #ice thickness

## Creation of model for simulation of wind and water velocities and then,
##compare the resulting ice velocities with buoy velocoities.
def model(x,ti):
    xc = x[0]
    yc = x[1]
    u  = x[2]
    v  = x[3]
    dx=np.zeros(4)
    dx[0] = d_xc=u
    dx[1] = d_yc=v
    u_mag = np.sqrt((u-Uwvec[ti,0])**2+(v-Uwvec[ti,1])**2)
    dx[2] = d_u  = (r_ice_water*rho_water*(Uwvec[ti,0]-u)*u_mag + 
           r_air_ice*rho_air*Uavec[ti,0]**2)/(rho_ice*h)+((f*v)) 
    dx[3] = d_v  = (r_ice_water*rho_water*(Uwvec[ti,1]-v)*u_mag+
           r_air_ice*rho_air*Uavec[ti,1]**2)/(rho_ice*h)-((f*u))
    return dx

#run simulation

times = np.arange(tstart,tstop,tstep)
x=[25.205, 79.1288, -4.44e-06, -1.1e-06] #initial values
results=np.zeros((len(x),len(times)))
ti=0
#global Uwvec
#global Uavec
for t in times:
    global x
    global ti
    print(t,x)
    results[:,ti]=x[:]
    xn = x + tstep*model(x,ti) 
    x  = xn
    ti+=1

#simulated ice positions and velocity
Xis=results[0,:]
Yis=results[1,:]
Uis=results[2,:]
Vis=results[3,:]


#plotting of simulated ice velocity versus buoy data

fig=plt.figure(figsize=(12, 7))
plt.subplot(1,2,1)
plt.plot(Uis[::4],color='r',label='buoy_drift')
plt.plot(Uw[::4],color='b',label='GTSM_vel')
plt.plot(Ua,color='y',label='wind_vel')

plt.legend(loc=1)               
plt.subplot(1,2,2)
plt.plot(Vis[::4],color='r',label='buoy_drift')
plt.plot(Vw[::4],color='b',label='GTSM_vel')
plt.plot(Va,color='y',label='wind_vel')
plt.legend(loc=1)
# general_Settings_for model simulations


import numpy as np
import logging 
def settings():
 s=dict()
 s['rho_ice'] = 920.0
 s['rho_air'] = 1.2
 s['rho_water'] = 1024.0
 s['g']=9.81
 s['seconds']=1.0 #si unit
 s['minutes']= 60.0*s['seconds']
 s['hours']  = 60.0*s['minutes']  
 s['days']   = 24.0*s['hours']
 s['rad2deg']= 180.0/np.pi
 s['deg2rad']= np.pi/180.0
 s['thetaa'] = 0*s['deg2rad']#angle made by geostrophic wind with ice velocity.
 s['beta']=25*s['deg2rad']
 s['thetaw'] = 0*s['deg2rad']                       #np.pi - s['beta']
 s['omega']  = (2.0*np.pi)/(24.0*s['hours'])
 s['tmplier']=1/30
 s['dt'] = 15*s['tmplier']*s['minutes']
 s['iCa']=1.2e-03
 s['iCw']=5.5e-03
 s['Nu0']=0
 s['Nv0']=0
 s['h']=1
 s['trate']= 0     #-1.1e-07
 s['l']=0
 s['n']=4 #number of elements of model state vector.
 s['en']=s['n']+2*(1+s['l']*int(1/s['tmplier'])) #number of elements of extended state vector +2 is for Nu and Nv
 # s['en']=s['n']+2+2
 logging.info("Time step is:" +str(s['dt'])+' seconds')
 logging.info("Ice Air drag coefficient is:" +str(s['iCa']))
 logging.info("Ice Water drag coefficient is:" +str(s['iCw']))
 s['mod']='RungeKutta2'   #'ExplicitEuler'
 s['tidepg']='GTSM'    #'FES2014'
 s['icethck']='modeldef'  #'owndef'
 s['tidedict']= {'M2':12.421,'S2':12.,'MU2':12.871,'O1':25.819,'K1':23.934,'M4':6.2103}
 s['cordict']={'deg1':74.7,'deg2':79}  #latitude for coriolis
 return(s)
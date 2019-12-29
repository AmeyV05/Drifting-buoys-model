# script for creating nc file from a 0006 map file to obtain data only for 0-40 longitude and 73-81 latitudes. 

import pandas as pd
import numpy as np
import netCDF4 as nc4

GTSM_data= nc4.Dataset("P:/1230882-emodnet_hrsm/fromAmey/Modelling_buoy_analysis/gtsm3_5km_Standard/output/gtsm_fine_std_0006_map.nc")
Xw=np.array(GTSM_data.variables["FlowElem_xcc"]) #X is longitdue and Y is latitude
Yw=np.array(GTSM_data.variables["FlowElem_ycc"])

#tw=np.array(GTSM_data.variables["time"])

Xwn=[]
Ywn=[]
for i in range(len(Xw)):
    if (Xw[i]>=0 and Xw[i]<=40 and Yw[i]>=73 and Yw[i]<=81):
        Xwn=np.append(Xwn,Xw[i])
        Ywn=np.append(Ywn,Yw[i])
print(np.shape(Xwn.T))


for j in range(4):
    Uw=np.array(GTSM_data.variables["ucx"][j,:])
    Vw=np.array(GTSM_data.variables["ucy"][j,:])
    Uwn=[]
    Vwn=[]
    for i in range(len(Xw)):
        if (Xw[i]>=0 and Xw[i]<=40 and Yw[i]>=73 and Yw[i]<=81):
            Uwn=np.append(Uwn,Uw[i])
            Vwn=np.append(Vwn,Vw[i])
    if (j==0):
        Uwtotal=Uwn
        Vwtotal=Vwn
    else:
        Uwtotal=np.vstack((Uwtotal,Uwn))
        Vwtotal=np.vstack((Vwtotal,Vwn))
    print(j)



#creation of an Nc file with GTSM dataset only around certain coordinates.
gtsmwrite = nc4.Dataset('gtsm_truncated.nc','w', format='NETCDF4') #'w' stands for write #f is a netcdf dataset object.
velgrp = gtsmwrite.createGroup('GTSM Tidal Velocity truncated data')
velgrp.createDimension('nFlowElem', len(Xwn))
#velgrp.createDimension('lat', len(Ywn))
velgrp.createDimension('time', len(tw))
longitude = velgrp.createVariable('FlowElem_xcc', 'f8', 'nFlowElem')
latitude = velgrp.createVariable('FlowElem_ycc', 'f8', 'nFlowElem')  
ucx = velgrp.createVariable('ucx', 'f8', ('time',  'nFlowElem'))
ucy = velgrp.createVariable('ucy', 'f8', ('time',  'nFlowElem'))
time = velgrp.createVariable('Time', 'f8', 'time')
longitude[:] = Xwn #The "[:]" at the end of the variable instance is necessary
latitude[:] = Ywn
ucx[:,:] = Uwtotal
ucy[:,:] = Vwtotal
time[:]=tw

#Add global attributes
gtsmwrite.description = "Truncated dataset from the GTSM map file"
gtsmwrite.history = "Created on 26/12/2019" # today.strftime("%d/%m/%y")

#Add local attributes to variable instances
longitude.units = 'degrees east'
latitude.units = 'degrees north'
time.units = 'seconds since 2014-03-01 00:00:00'
ucx.units = 'm s-1'
ucy.units= 'm s-1'
#ucx._FillValue=-999
gtsmwrite.close()
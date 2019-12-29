#module to load data from GTSM

import netCDF4 as nc4 
import numpy as np
import generalfunc as gf
np.set_printoptions(threshold=np.inf)

#Note that data only available
#The 1344 gets data from 15th march to end 
def read_GTSM_his(file,sindex,eindex):	
 GTSM_data= nc4.Dataset(file)
 Xt=np.array(GTSM_data.variables["station_x_coordinate"])[sindex:eindex,1] #X is longitdue and Y is latitude
 Yt=np.array(GTSM_data.variables["station_y_coordinate"])[sindex:eindex,1]
 Ut=np.array(GTSM_data.variables["x_velocity"])[sindex:eindex,1]
 Vt=np.array(GTSM_data.variables["y_velocity"])[sindex:eindex,1]
 Tt=np.array(GTSM_data.variables["time"])
 Tt=gf.num2datetimesecs(2014,3,1,sindex,eindex,Tt)
 print("Completed reading GTSM data")
 return (Xt,Yt,Ut,Vt,Tt)

def read_GTSM_map(file):
 GTSM_data= nc4.Dataset(file)
 gtsm_grp=GTSM_data.groups['GTSM Tidal Velocity truncated data']
 Xt=np.array(gtsm_grp.variables["FlowElem_xcc"]) #X is longitdue and Y is latitude
 Yt=np.array(gtsm_grp.variables["FlowElem_ycc"])
 Tt=np.array(gtsm_grp.variables["Time"])[1:]
 Ut=np.array(gtsm_grp.variables["ucx"])[1:,:]
 Vt=np.array(gtsm_grp.variables["ucy"])[1:,:]
 print("Completed reading GTSM data")
 return (Xt,Yt,Ut,Vt,Tt)
 


#reading ERA5 winds data
#obtaining u and v velocity of winds 
# note winds always start from 1900 1 1
def read_wind(file):
 wind_data_u = nc4.Dataset(file)
 Xa=np.array(wind_data_u.variables["longitude"])
 Ya=np.array(wind_data_u.variables["latitude"])
 Ta=np.array(wind_data_u.variables["time"])
 u10=np.array(wind_data_u.variables["u10"])
 v10=np.array(wind_data_u.variables["v10"])
 print("Completed reading wind data")
 return (Xa,Ya,u10,v10,Ta)

#Ocean_currents_CMEMS
#obtaining u and v velocity of ocean
#note ocean starts from 1950 
def read_ocean(file):
 ocean_data = nc4.Dataset(file)
 Xo=np.array(ocean_data.variables["longitude"])
 Yo=np.array(ocean_data.variables["latitude"])
 To=np.array(ocean_data.variables["time"])
 uo=np.array(ocean_data.variables["uo"])
 vo=np.array(ocean_data.variables["vo"])
 print("Completed reading ocean currents data")
 return (Xo,Yo,uo,vo,To)


def main():
 file="../Data_from_models/era5_wind_201403_05.nc"
 [Xw,Yw,Ut,Vt,Tt]=read_wind(file)
 print(len(Xw),len(Yw))

if __name__ == '__main__':
 main()

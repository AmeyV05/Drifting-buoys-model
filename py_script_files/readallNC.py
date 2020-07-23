#module to load data from all the data files.

import netCDF4 as nc4 
import numpy as np
import generalfunc as gf
import logging 
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
 logging.info("Completed reading GTSM data")
 return (Xt,Yt,Ut,Vt,Tt)

def read_GTSM_map(file):
 GTSM_data= nc4.Dataset(file)
 gtsm_grp=GTSM_data.groups['GTSM Tidal truncated data']
 Xt=np.array(gtsm_grp.variables["FlowElem_xcc"]) #X is longitdue and Y is latitude
 Yt=np.array(gtsm_grp.variables["FlowElem_ycc"])
 Tt=np.array(gtsm_grp.variables["Time"]) #[1:]
 Ut=np.array(gtsm_grp.variables["ucx"])# [1:,:]
 Vt=np.array(gtsm_grp.variables["ucy"])# [1:,:]
 ssht=np.array(gtsm_grp.variables["s1"])
 logging.info("Completed reading GTSM data")
 return (Xt,Yt,Ut,Vt,ssht,Tt)
 


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
 logging.info("Completed reading wind data")
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
 ssh=np.array(ocean_data.variables["zos"])
 logging.info("Completed reading ocean currents data")
 return (Xo,Yo,uo,vo,ssh,To)

def read_ice(file):
 ice_data=nc4.Dataset(file)
 Xe=np.array(ice_data.variables["longitude"])
 Ye=np.array(ice_data.variables["latitude"])
 Te=np.array(ice_data.variables["time"])
 hi=np.array(ice_data.variables["sithick"])
 logging.info("Completed reading ice parameter data")
 return(Xe,Ye,hi,Te)

def read_FES(file,Bnum):
	switcher={'02':[1.2,0.24],
				 '03':[2.2,0.24],
				 '09':[1.2,0.24],
				 '13':[1.2,0.24],
				 '14':[1.0,0.2],
				 '16':[1.5,0.24]}
	dlonlat=switcher.get(Bnum,"Invalid Buoy number")
	efile=file+Bnum+"e.nc"
	fesdatae=nc4.Dataset(efile)
	Tft=np.array(fesdatae.variables["time"])
	the=np.array(fesdatae.variables["tide"])
	the=the.diagonal()
	wfile=file+Bnum+"w.nc"
	fesdataw=nc4.Dataset(wfile)
	thw=np.array(fesdataw.variables["tide"])
	thw=thw.diagonal()
	nfile=file+Bnum+"n.nc"
	fesdatan=nc4.Dataset(nfile)
	thn=np.array(fesdatan.variables["tide"])
	thn=thn.diagonal()
	sfile=file+Bnum+"s.nc"
	fesdatas=nc4.Dataset(sfile)
	ths=np.array(fesdatas.variables["tide"])
	ths=ths.diagonal()
	tht=np.vstack((the,thw,thn,ths))
	print(np.shape(tht.T))
	logging.info("Completed reading FES2014 data")
	return(tht.T,Tft,dlonlat)

def readingalldata(fileloc,Bnum):

 ### FES2014 data
 logging.info("Reading FES2014 data")
 file=fileloc+"/FES2014data/BUOY_"
 [tht,Tft,dlonlat]=read_FES(file,Bnum)
 SD={'th':tht,'Tft':Tft,'dlola':dlonlat}
 #####GTSM data
 logging.info("Reading GTSM data")
 #print("Reading GTSM data")
 # file=fileloc+"/gtsm_truncated.nc"
 file=fileloc+"/GTSM_final_truncated.nc"
 # file=r"P:/1230882-emodnet_hrsm/fromAmey/Modelling_buoy_analysis/GTSM_final_truncated.nc"
 [Xt,Yt,ut,vt,ssht,Tt]=read_GTSM_map(file)
 TD={'Xt': Xt,'Yt': Yt,'ut': ut,'vt': vt,'ssht': ssht,'Tt': Tt}
 #### ERA 5 winds
 logging.info("Reading ERA5 wind data.")
 #ERA 5 winds
 file=fileloc+"/era5_wind_201403_05.nc"
 [Xa,Ya,u10,v10,Ta]=read_wind(file)
 AD={'Xa': Xa,'Ya': Ya,'u10': u10,'v10': v10,'Ta': Ta}
 ## Ocean currents
 logging.info("Reading ocean currents data.")
 # file=fileloc+"/Ocean_currents_buoy.nc"
 file=fileloc+"/CMEMS_ocean_data.nc"
 [Xo,Yo,uo,vo,ssh,To]=read_ocean(file)
 OD={'Xo': Xo,'Yo': Yo,'uo': uo,'vo': vo,'ssh':ssh,'To': To}
 logging.info("Reading Ice thickness data.")
 file=fileloc+"/ice-param-mercator.nc"
 [Xe,Ye,hi,Te]=read_ice(file)
 ID={'Xe': Xe,'Ye': Ye,'hi':hi,'Te': Te}
 FD={'TD':TD,'AD':AD,'OD':OD,'ID':ID,'SD':SD}
 return(FD)

def main():
 file="../Data_from_models/era5_wind_201403_05.nc"
 [Xw,Yw,Ut,Vt,Tt]=read_wind(file)
 print(len(Xw),len(Yw))

if __name__ == '__main__':
 main()

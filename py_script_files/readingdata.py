####module to load data from all the .nc files of the various models.
import netCDF4 as nc4 
import numpy as np
import pandas as pd
import generalfunc as gf
import logging 


### wW define different set of functions to read the different nc files. 
# The GTSM map file is truncated to get data between. 0-40 deg lon and 70-80 deg latitude.

def read_GTSM_map(file1,file2):
	GTSM_data= nc4.Dataset(file1)
	gtsm_grp=GTSM_data.groups['GTSM Tidal truncated data']
	Xt=np.array(gtsm_grp.variables["FlowElem_xcc"]) #X is longitdue and Y is latitude
	Yt=np.array(gtsm_grp.variables["FlowElem_ycc"])
	Tt=np.array(gtsm_grp.variables["Time"]) 
	Ut=np.array(gtsm_grp.variables["ucx"])
	Vt=np.array(gtsm_grp.variables["ucy"])
	ssht=np.array(gtsm_grp.variables["s1"]) #sea surface heights.
	GTSM_data= nc4.Dataset(file2)
	gtsm_grp=GTSM_data.groups['GTSM Tidal truncated data']
	wdt=np.array(gtsm_grp.variables["waterdepth"])
	logging.info("Completed reading GTSM data")
	return (Xt,Yt,Ut,Vt,ssht,wdt,Tt)
 

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


# This function reads the ice thickness data from the .nc file 
# provided by the CMEMS ocean model. Same one used for ocean currents but the 
# .nc files are separate. And thickness data only available on daily basis not hourly.
def read_ice(file):
	ice_data=nc4.Dataset(file)
	Xe=np.array(ice_data.variables["longitude"])
	Ye=np.array(ice_data.variables["latitude"])
	Te=np.array(ice_data.variables["time"])
	hi=np.array(ice_data.variables["sithick"])
	logging.info("Completed reading ice parameter data")
	return(Xe,Ye,hi,Te)

# Reading FES2014 data for sea surface heights.
# This data is created based on the trajectories e,w,n,s of each buoy trajectory. 
def read_FES(file,Bnum):
	switcher={'02':[1.2,0.24],
				 '03':[2.2,0.24],
				 '09':[1.2,0.24],
				 '13':[1.2,0.24],
				 '14':[1.0,0.2],
				 '16':[1.2,0.24]}
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
	logging.info("Completed reading FES2014 data")
	#tht corresponds to tidal heights.
	return(tht.T,Tft,dlonlat)

# Reading and (processing) buoy data.
# This is the only function where we process the buoy data in order to obtain the velocities in m/s.
def read_bdata(file,Bnum):
	prefix="BUOY_"
	bname=prefix+Bnum
	logging.info("Reading data for "+bname)
	file=file+bname+'.csv'
	D=pd.read_csv(file)
	Xib=np.array(D['Lon'])
	Yib=np.array(D['Lat'])
	Tib=np.array(D['Date(GMT)'])
	D.head()
	#velocity computation using the buoy data
	dt=15*60 #sampling time step for each buoy observation.
	Uib=[];Vib=[]
	for j in range(len(Xib)-1):
		dlon=Xib[j+1]-Xib[j]
		dlat=Yib[j+1]-Yib[j]
		lat=Yib[j];lon=Xib[j]
		# (dlatm,dlonm)=gf.latlon2meters(lat,lon,dlat,dlon)
		(dlatm,dlonm)=gf.latlon2meters(lat,dlat,dlon)
		Uib=np.append(Uib,(dlonm)/dt)
		Vib=np.append(Vib,(dlatm)/dt)
	return(Xib,Yib,Uib,Vib,Tib)

# main function to read all the data using the above functions. 
def readingdata(fileloc,Bnum):
	## Reading buoy data
	file =fileloc+"/buoy_data/"
	[Xib,Yib,Uib,Vib,Tib]=read_bdata(file,Bnum)
	BD={'Xib': Xib,'Yib': Yib,'Uib': Uib,'Vib': Vib,'Tib': Tib} #dictionary for all the vectors.
	logging.info("Completed reading buoy data.")
	### FES2014 data 
	logging.info("Reading FES2014 data")
	file=fileloc+"/FES2014data/BUOY_"
	[tht,Tft,dlonlat]=read_FES(file,Bnum)
	SD={'th':tht,'Tft':Tft,'dlola':dlonlat} 

	#####GTSM data
	logging.info("Reading GTSM data")
	file1=fileloc+"/GTSM_final_truncated.nc"
	file2=fileloc+"/gtsm_truncated_wd.nc"
	[Xt,Yt,ut,vt,ssht,wdt,Tt]=read_GTSM_map(file1,file2)
	TD={'Xt': Xt,'Yt': Yt,'ut': ut,'vt': vt,'ssht': ssht,'wdt':wdt,'Tt': Tt}

	#### ERA 5 winds
	logging.info("Reading ERA5 wind data.")
	file=fileloc+"/era5_wind_201403_05.nc"
	[Xa,Ya,u10,v10,Ta]=read_wind(file)
	AD={'Xa': Xa,'Ya': Ya,'u10': u10,'v10': v10,'Ta': Ta}

	## Ocean currents
	logging.info("Reading ocean currents data.")
	file=fileloc+"/CMEMS_ocean_data.nc"
	[Xo,Yo,uo,vo,ssh,To]=read_ocean(file)
	OD={'Xo': Xo,'Yo': Yo,'uo': uo,'vo': vo,'ssh':ssh,'To': To}
	## Ice thickness.
	logging.info("Reading Ice thickness data.")
	file=fileloc+"/ice-param-mercator.nc"
	[Xe,Ye,hi,Te]=read_ice(file)
	ID={'Xe': Xe,'Ye': Ye,'hi':hi,'Te': Te}
	# Clubbing together all dictionaries
	FD={'TD':TD,'AD':AD,'OD':OD,'ID':ID,'SD':SD,'BD':BD}
	return(FD)

def main():
	gendataloc='../../generated_data'
	gf.logcreate(gendataloc)
	Bnum='02'
	fileloc="../../container/Model_n_buoy_data"
	print("Running for buoy: "+Bnum)
	readingdata(fileloc,Bnum)

if __name__ == '__main__':
	main()


## obsolete functions.

# #Note that data only available
# #The 1344 gets data from 15th march to end 
# def read_GTSM_his(file,sindex,eindex):	
#  GTSM_data= nc4.Dataset(file)
#  Xt=np.array(GTSM_data.variables["station_x_coordinate"])[sindex:eindex,1] #X is longitdue and Y is latitude
#  Yt=np.array(GTSM_data.variables["station_y_coordinate"])[sindex:eindex,1]
#  Ut=np.array(GTSM_data.variables["x_velocity"])[sindex:eindex,1]
#  Vt=np.array(GTSM_data.variables["y_velocity"])[sindex:eindex,1]
#  Tt=np.array(GTSM_data.variables["time"])
#  Tt=gf.num2datetimesecs(2014,3,1,sindex,eindex,Tt)
#  logging.info("Completed reading GTSM data")
#  return (Xt,Yt,Ut,Vt,Tt)
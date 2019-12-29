## plotting functions for the model 

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf 
import numpy as np


## plotting of x and y ice buoy velocities versus simulated drifts. 

def plticevel(Uisvec,Uibvec,path):
	fig=plt.figure(figsize=(12, 7))
	plt.subplot(1,2,1)
	plt.plot(Uibvec[:,0],color='r',label='buoy_drift')
	plt.plot(Uisvec[:,0],color='g',linestyle=":",label='sim_ice_drift')
	plt.xlabel('Time')
	plt.ylabel('U')
	plt.legend(loc=1)               
	plt.subplot(1,2,2)
	plt.plot(Uibvec[:,1],color='r',label='buoy_drift')
	plt.plot(Uisvec[:,1],color='g',linestyle=":",label='sim_ice_drift')
	plt.xlabel('Time')
	plt.ylabel('V')
	plt.legend(loc=1)
	plt.savefig(path+'/velplot.jpg', format='jpg', dpi=500)
	#plt.show()
	plt.close(fig)


def plticepos(Xib,Yib,Xis,Yis,path):
 fig=plt.figure(figsize=(12, 12))
 ax=plt.axes(projection=ccrs.LambertAzimuthalEqualArea(central_longitude=25.0,central_latitude=77.0)) 
 #ax.set_extent([-15,45,65,85])         
 #feature=cpf.GSHHSFeature(scale='i',levels=[1],facecolor='#AAAAAA')
 #ax.add_feature(feature) 
 #ax.add_feature(cpf.OCEAN, facecolor='blue')
 gl=ax.gridlines()
 plt.plot(Xib,Yib,color='red',transform=ccrs.PlateCarree(),label='buoy drift')
 plt.plot(Xis,Yis,color='yellow',transform=ccrs.PlateCarree(),label='simulated ice drift')
 ax.set_extent([15,33,74,81]) 
 service='https://www.gebco.net/data_and_products/gebco_web_services/web_map_service/mapserv?'
 ax.add_wms(service,layers=['GEBCO_LATEST'],wms_kwargs={'width':900*2,'height':600*2})
 plt.legend()
 plt.savefig(path+'/ice_drift.jpg',dpi=500)
 plt.close(fig) 


import cartopy.crs as ccrs
import cartopy.feature as cpf 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import scipy.io

#def main():
 # ax=figure().gca(projection=cartopy.crs.LambertAzimuthalEqualArea(central_longitude=40.0,central_latitude=77.0, globe=None)) 
 # ax.add_feature(cpf.LAND)
 # ax.add_feature(cpf.OCEAN)
 # ax.add_feature(cpf.COASTLINE)
 # ax.add_feature(cpf.BORDERS, linestyle=':')
 # ax.add_feature(cpf.LAKES, alpha=0.5)
 # ax.add_feature(cpf.RIVERS)
 # ax.coastlines()
 # ax.set_extent([15,30,74,80])
 
if __name__ == '__main__':
 #main()


 D=pd.read_csv('CliSAP_Boje_16.csv')
 Lon=np.array(D['Lon'])
 Lat=np.array(D['Lat'])
 D.head()

 fig=plt.figure()
 ax=plt.axes(projection=ccrs.LambertAzimuthalEqualArea(central_longitude=25.0,central_latitude=77.0)) 
 ax.set_extent([15,35,74,80]) 
 #ax.set_extent([-15,45,65,85])         
 feature=cpf.GSHHSFeature(scale='i',levels=[1],facecolor='#AAAAAA')
 ax.add_feature(feature) 
 ax.add_feature(cpf.OCEAN, facecolor='blue')
 gl=ax.gridlines()
# gl.xlines=True
# gl.xlabels_top=True
# gl.xlocator=mticker.FixedLocator([15,20,25,30,35,40])
# gl.xformatter=LONGITUDE_FORMATTER
 plt.plot(Lon,Lat,color='red',transform=ccrs.PlateCarree())
 mat = scipy.io.loadmat('ForAmey.mat')
 Lon_bath=mat['LONSRC']
 Lat_bath=mat['LATSRC']
 Dep_bath=mat['DPS_SRC']
 levels=np.linspace(0,200,9)
 plt.contourf(Lon_bath, Lat_bath, Dep_bath,levels,
             transform=ccrs.PlateCarree())
 plt.colorbar()
 plt.savefig('plot.png',dpi=500) 
 plt.show()
 # feature=cpf.NaturalEarthFeature(name='land',category='physical',scale='10m',
#                                 facecolor='#AAAAAA')                                                      
 #ax.stock_img()

# ax.add_feature(feature)
# feature=cpf.NaturalEarthFeature(name='ocean',category='physical',scale='10m',
#                                 facecolor='g')
 
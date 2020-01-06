## plotting functions for the model 
from copy import copy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf 
import numpy as np
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom


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
 fig=plt.figure(figsize=(12, 12), frameon=True)
 ax=plt.axes(projection=ccrs.LambertAzimuthalEqualArea(central_longitude=25.0,central_latitude=77.0)) 
 ax.set_extent([15,33,74,81]) 
 # Define gridline locations and draw the lines using cartopy's built-in gridliner:
 # *must* call draw in order to get the axis boundary used to add ticks:
 fig.canvas.draw()
 xticks = [ 0, 5, 10, 15, 20, 25, 30, 35, 40]
 yticks = [72, 74, 76, 78, 80, 82]
 ax.gridlines(xlocs=xticks, ylocs=yticks)
 # Label the end-points of the gridlines using the custom tick makers:
 ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
 ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
 lambert_xticks(ax, xticks)
 lambert_yticks(ax, yticks)

 # Plotting buoy obs and buoy simulated locations
 plt.plot(Xib,Yib,color='red',transform=ccrs.PlateCarree(),label='buoy drift')
 plt.plot(Xis,Yis,color='yellow',transform=ccrs.PlateCarree(),label='simulated ice drift') 
 # gebco wms background 
 service='https://www.gebco.net/data_and_products/gebco_web_services/web_map_service/mapserv?'
 ax.add_wms(service,layers=['GEBCO_LATEST'],wms_kwargs={'width':900*2,'height':600*2})
 plt.legend()
 plt.savefig(path+'/ice_drift.jpg',dpi=500)
 plt.close(fig) 





## functions used to label x and y axis in cartopy for lambert conformal projections. 
# Adopted from cartopy github

def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.
    
    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)],}
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    te = lambda xy: xy[0]
    lc = lambda t, n, b: np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels])
    

def lambert_yticks(ax, ticks):
    """Draw ricks on the left y-axis of a Lamber Conformal projection."""
    te = lambda xy: xy[1]
    lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels])

def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """Get the tick locations and labels for an axis of a Lambert Conformal projection."""
    outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(ccrs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(ccrs.Geodetic(), xy[:, 0], xy[:, 1])
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:    
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels.pop(index)
    return _ticks, ticklabels
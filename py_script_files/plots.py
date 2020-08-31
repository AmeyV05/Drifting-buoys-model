## plotting functions for the model 
import numpy as np
import matplotlib.pyplot as plt
import settings
import cartopy.crs as ccrs
import cartopy.feature as cpf 
import matplotlib.ticker as mticker
import matplotlib.transforms as transforms
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import lineid_plot
from copy import copy
import shapely.geometry as sgeom
## plotting of x and y ice buoy velocities versus simulated drifts. 

def plticevel(Uisvec,Uibvec,Tib,path):
  fig=plt.figure(figsize=(16, 5))
  plt.subplot(1,2,1)
  plt.plot(Tib,Uibvec[:,0],color='r',label='buoy_drift')
  plt.plot(Tib,Uisvec[:,0],color='g',linestyle=":",label='sim_ice_drift')
  plt.xlabel('Time (s)',fontweight='bold')
  plt.ylabel('U (m/s)',fontweight='bold')
  labels=Tib[::800] 
  plt.xticks(labels, labels, rotation='0')
  plt.legend(loc=1)               
  plt.subplot(1,2,2)
  plt.plot(Tib,Uibvec[:,1],color='r',label='buoy_drift')
  plt.plot(Tib,Uisvec[:,1],color='g',linestyle=":",label='sim_ice_drift')
  plt.xlabel('Time (s)',fontweight='bold')
  plt.ylabel('V (m/s)',fontweight='bold')
  labels=Tib[::800] 
  plt.xticks(labels, labels, rotation='0')
  plt.legend(loc=1)
  plt.savefig(path+'/velplot.jpg', format='jpg')
  plt.close(fig)


def plticepos(Xib,Yib,Xis,Yis,path):
  fig=plt.figure(figsize=(12, 12), frameon=True)
  ax=plt.axes(projection=ccrs.LambertAzimuthalEqualArea(central_longitude=25.0,central_latitude=77.0)) 
  # ax.set_extent([15,33,74,81]) 
  ax.set_extent([16,28,74,78]) 
  # Define gridline locations and draw the lines using cartopy's built-in gridliner:
  # *must* call draw in order to get the axis boundary used to add ticks:
  fig.canvas.draw()
  xticks = [  0, 4, 12, 16, 20, 24, 28, 32, 36]
  yticks = [72, 74, 76, 78, 80, 82]
  ax.gridlines(xlocs=xticks, ylocs=yticks)
  # Label the end-points of the gridlines using the custom tick makers:
  ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
  ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
  lambert_xticks(ax, xticks)
  lambert_yticks(ax, yticks)

  # Plotting buoy obs and buoy simulated locations
  plt.plot(Xib,Yib,color='red',transform=ccrs.PlateCarree(),label='buoy drift')
  plt.plot(Xis,Yis,'--',color='yellow',transform=ccrs.PlateCarree(),label='simulated ice drift') 
  # gebco wms background 
  service='https://www.gebco.net/data_and_products/gebco_web_services/web_map_service/mapserv?'
  ax.add_wms(service,layers=['GEBCO_LATEST'],wms_kwargs={'width':900*2,'height':600*2})
  feature=cpf.GSHHSFeature(scale='i',levels=[1],facecolor='#e6e1e1',alpha=1)
  ax.add_feature(feature) 
  plt.legend(prop={"size":16},framealpha=1)
  plt.savefig(path+'/ice_drift.jpg')
  plt.close(fig) 


def pltFT(path,name,xb,xs,units,tvec,argval):
  s=settings.settings()
  fig=plt.figure(figsize=(21,7))
  ax1=plt.subplot2grid((2,1), (1,0), rowspan=1, colspan=1)
#     ax1=fig.add_axes([0.1,0.1, 0.8, 0.8])
  l1=ax1.plot(tvec,xb[1,:],color='r',label="obs")
  l2=ax1.plot(tvec,xs[1,:],color='g',label="sim") 
  handles, labels = ax1.get_legend_handles_labels()
  ax1.legend(handles,labels,loc=1,fontsize='12')
  ax1.set_xlim([5,27])
  ak = lineid_plot.initial_annotate_kwargs()
  tide=np.array(list(s['tidedict'].keys()));cor=np.array(list(s['cordict'].values()))
  tidlab=np.append(tide,cor)
  lineid_plot.plot_line_ids(tvec, xb[1,:], argval, tidlab,ax=ax1,max_iter=30000)
  plt.xlabel('Period [h]',fontsize=12,fontweight='bold')
  plt.ylabel('Phase (deg)',fontsize=12,fontweight='bold')
  ax2=fig.add_axes([0.125, 0.57, 0.775, 0.34],sharex=ax1)
  l1=ax2.plot(tvec,xb[0,:],color='r',label="obs")
  l2=ax2.plot(tvec,xs[0,:],color='g',label="sim")
  handles, labels = ax2.get_legend_handles_labels()
  ax2.legend(handles,labels,loc=1,fontsize='12')
  lineid_plot.plot_line_ids(tvec, xb[0,:], argval, tidlab,ax=ax2,max_iter=30000) #
  plt.ylabel('Amplitude ('+units+')',fontsize=12,fontweight='bold')
  plt.setp(ax2.get_xticklabels(), visible=False)
  plt.xticks(fontsize=12);plt.yticks(fontsize=12)
  plt.savefig(path+'/'+name+'.jpg',format='jpg')
  plt.close(fig)



def pltfilsig(p1_x,p_xfilter,p_xresiduum,path,name):
  fig = plt.figure(figsize=(12, 3))
  plt.subplot(1,2,1)
  p_x_ = [v for i, v in enumerate(p1_x) if i % 10 == 0]
  p_xfilter_ = [v for i, v in enumerate(p_xfilter) if i % 10 == 0]
  p_xresiduum_ = [v for i, v in enumerate(p_xresiduum) if i % 10 == 0]
  plt.plot(p_x_,'r-',label='Longitude')
  plt.plot(p_xfilter_,'k-',label='filter_Signal')
  plt.xlabel('Time')
  plt.ylabel('Longitude')
  plt.legend(loc=1)
  plt.subplot(1,2,2)
  plt.plot(p_xresiduum_,label='Residual')
  plt.xlabel('Time')
  plt.legend(loc=1)
  plt.savefig(path+'/'+name+'.jpg',format='jpg', bbox_inches='tight')
  plt.close(fig)



# functions used to label x and y axis in cartopy for lambert conformal projections. 
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



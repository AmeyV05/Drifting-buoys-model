## plotting functions for the model 
from copy import copy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cpf 
import numpy as np
import matplotlib.ticker as mticker
import matplotlib.transforms as transforms
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom
import lineid_plot
plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": [
         r"\usepackage[utf8x]{inputenc}",
         r"\usepackage[T1]{fontenc}",
         r"\usepackage{cmbright}",
         ]
})
params = {'legend.fontsize': '16',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'16',
         'ytick.labelsize':'16',
         "font.weight":'bold'}
plt.rcParams.update(params)
## plotting of x and y ice buoy velocities versus simulated drifts. 

def plticevel(Uisvec,Uibvec,tmplierinv,path):
 Uibvec=np.repeat(Uibvec,tmplierinv,axis=0)
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
 plt.savefig(path+'/velplot.jpg', format='jpg')
 #plt.show()
 plt.close(fig)


def plticepos(Xib,Yib,Xis,Yis,path):
 fig=plt.figure(figsize=(12, 12), frameon=True)
 ax=plt.axes(projection=ccrs.LambertAzimuthalEqualArea(central_longitude=25.0,central_latitude=77.0)) 
 # ax.set_extent([15,33,74,81]) 
 ax.set_extent([16,28,74,78]) 
 # Define gridline locations and draw the lines using cartopy's built-in gridliner:
 # *must* call draw in order to get the axis boundary used to add ticks:
 fig.canvas.draw()
 # xticks = [ 0, 5, 10, 15, 20, 25, 30, 35, 40]
 # yticks = [72, 74, 76, 78, 80, 82]
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


def pltalpha(alpha,path,name):
 fig=plt.figure(figsize=(12,7),frameon=True)
 plt.plot(alpha,label=name)
 plt.xlabel('time')
 plt.ylabel(name)
 plt.savefig(path+'/'+name+'.jpg',format='jpg')
 plt.close(fig)

def pltFT(path,name,xb,xs,tvec,arg_tide,arg_cor):
 fig=plt.figure(figsize=(24,10),frameon=True)
 ax1=plt.subplot(2,1,1)
 ax1.plot(tvec,xb[0,:],color='r',label=name+"_obs")
 ax1.plot(tvec,xs[0,:],color='g',label=name+"_sim")
 ax1.set_xlim([5,27])
 argdeg1=arg_cor['argdeg1'];argdeg2=arg_cor['argdeg2']
 # plt.axvline(tvec[arg_tide['M2']],color='m', linewidth=1, linestyle='dashed',label='M2')
 # plt.text(tvec[arg_tide['M2']], .5, 'M2', transform=fig.transFigure)
 # plt.axvline(tvec[arg_tide['S2']],color='y',  linewidth=1, linestyle='dashed',label='S2')
 # plt.axvline(tvec[arg_tide['MU2']],color='g',  linewidth=1, linestyle='dashed',label='MU2')
 # plt.axvline(tvec[arg_tide['O1']],color='orange',  linewidth=1,  linestyle='dashed',label='O1')
 # plt.axvline(tvec[arg_tide['K1']],color='blue',  linewidth=1,  linestyle='dashed',label='K1')
 # plt.axvline(tvec[arg_tide['M4']],color='teal',  linewidth=1,  linestyle='dashed',label='M4')
 # plt.axvline(tvec[argdeg1],color='olive',  linewidth=1, linestyle='dashed',label='74.7')
 # plt.axvline(tvec[argdeg2], color='cyan', linewidth=1, linestyle='dashed',label='79')
 ak = lineid_plot.initial_annotate_kwargs()
 ak['arrowprops']['relpos'] = (0.5, 0)
 pk = lineid_plot.initial_plot_kwargs()
 pk['color'] = "magenta"
 argtid=[tvec[arg_tide['M2']],tvec[arg_tide['S2']],tvec[arg_tide['MU2']],
         tvec[arg_tide['O1']],tvec[arg_tide['K1']],tvec[arg_tide['M4']],
         tvec[argdeg1],tvec[argdeg2]]
 tidlab=['M2','S2','MU2','O1','K1','M4','74.7','79']
 lineid_plot.plot_line_ids(tvec, xb[0,:], argtid, tidlab,ax=ax1,plot_kwargs=pk,annotate_kwargs=ak,max_iter=1, box_loc=0)
 plt.xlabel('Period [h]',fontsize=12,fontweight='bold')
 plt.ylabel('Amplitude (deg)',fontsize=14,fontweight='bold')
 plt.xticks(fontsize=12);plt.yticks(fontsize=12)
 # plt.xlim([5,27])

 # plt.legend(loc=1)
 ax2=plt.subplot(2,1,2)
 ax2.plot(tvec,xb[1,:],color='r')
 ax2.plot(tvec,xs[1,:],color='g')
 ax2.set_xlim([5,27])
 argdeg1=arg_cor['argdeg1'];argdeg2=arg_cor['argdeg2']
 # plt.axvline(tvec[arg_tide['M2']],color='m', linewidth=1, linestyle='dashed')
 # plt.axvline(tvec[arg_tide['S2']],color='y',  linewidth=1, linestyle='dashed')
 # plt.axvline(tvec[arg_tide['MU2']],color='g',  linewidth=1, linestyle='dashed')
 # plt.axvline(tvec[arg_tide['O1']],color='orange',  linewidth=1,  linestyle='dashed')
 # plt.axvline(tvec[arg_tide['K1']],color='blue',  linewidth=1,  linestyle='dashed')
 # plt.axvline(tvec[arg_tide['M4']],color='teal',  linewidth=1,  linestyle='dashed')
 # plt.axvline(tvec[argdeg1],color='olive',  linewidth=1, linestyle='dashed')
 # plt.axvline(tvec[argdeg2], color='cyan', linewidth=1, linestyle='dashed')
 lineid_plot.plot_line_ids(tvec, xb[1,:], argtid, tidlab,ax=ax2,plot_kwargs=pk,annotate_kwargs=ak,max_iter=1, box_loc=0)
 plt.xlabel('Period [h]',fontsize=14,fontweight='bold')
 plt.ylabel('Phase (deg)',fontsize=14,fontweight='bold')
 # plt.xlim([5,27])
 plt.xticks(fontsize=12);plt.yticks(fontsize=12)
 # # fig.legend(bbox_to_anchor=(1.0, 0.89),fontsize=12)
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


def convplt(rmsposvec,tmvec,s,path,mod):
 xrms=rmsposvec[::2]
 yrms=rmsposvec[1::2]
 s['tmplier']=tmvec
 tvec = 15*s['tmplier']*s['minutes']
 fig=plt.figure(figsize=(12,12),frameon=True)
 plt.plot(tvec,xrms,'--bo',label='longitude rms')
 plt.plot(tvec,yrms,'--r*',label='latitude rms')
 plt.xlabel('time step size')
 plt.ylabel('rms error')
 plt.legend()
 plt.savefig(path+'/convergence_'+str(s['h'])+'_'+mod+'.jpg')
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
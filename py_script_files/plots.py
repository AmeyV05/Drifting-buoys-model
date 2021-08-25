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
  plt.plot(Tib,Uibvec[:,0],color='r',label='Obs')
  plt.plot(Tib,Uisvec[:,0],color='g',label='Mod')
  plt.xlabel('Time',fontweight='bold')
  plt.ylabel('U (m/s)',fontweight='bold')
  labels=Tib[::800] 
  plt.xticks(labels, labels, rotation='0',fontsize='11')
  plt.yticks(fontsize=11)
  plt.legend(loc=1)               
  plt.subplot(1,2,2)
  plt.plot(Tib,Uibvec[:,1],color='r',label='Obs')
  plt.plot(Tib,Uisvec[:,1],color='g',label='Mod')
  plt.xlabel('Time',fontweight='bold')
  plt.ylabel('V (m/s)',fontweight='bold')
  labels=Tib[::800] 
  plt.xticks(labels, labels, rotation='0',fontsize='11')
  plt.yticks(fontsize=11)
  plt.legend(loc=1)
  plt.savefig(path+'/velplot.jpg', format='jpg')
  plt.close(fig)


def plticepos(Xib,Yib,Xis,Yis,path,Bnum):
  prefix="BUOY_"
  bname=prefix+Bnum
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
  plt.plot(Xib,Yib,color='red',transform=ccrs.PlateCarree(),label='Obs')
  plt.plot(Xis,Yis,color='green',transform=ccrs.PlateCarree(),label='Mod') 
  # gebco wms background 
  service='https://www.gebco.net/data_and_products/gebco_web_services/web_map_service/mapserv?'
  ax.add_wms(service,layers=['GEBCO_LATEST'],wms_kwargs={'width':900*2,'height':600*2})
  feature=cpf.GSHHSFeature(scale='i',levels=[1],facecolor='#e6e1e1',alpha=1)
  ax.add_feature(feature) 
  plt.legend(prop={"size":18},framealpha=1)
  plt.savefig(path+'/sim_ice_drift'+bname+'.jpg',dpi=100)
  plt.close(fig) 


def pltFT(path,name,xs,xb,units,tvec,argval):
  s=settings.settings()
  fig=plt.figure(figsize=(21,7))
  ax1=plt.subplot2grid((2,1), (1,0), rowspan=1, colspan=1)
#     ax1=fig.add_axes([0.1,0.1, 0.8, 0.8])
  l1=ax1.plot(tvec,xb[1,:],color='r',label="obs")
  l2=ax1.plot(tvec,xs[1,:],color='g',label="sim") 
  handles, labels = ax1.get_legend_handles_labels()
  ax1.legend(handles,labels,loc=1,fontsize='12')
  ax1.set_xlim([5,23])
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

def pltFTn(path,name,xs,xb,units,tvec,argval):
  s=settings.settings()
  tidedict=s['tidedict'];cordict=s['cordict']
  Tcordict=[];Ttiddict=[]
  for k in cordict.keys():
      phi=np.deg2rad(cordict[k])
      f=2*s['omega']*np.sin(phi)
      T=2*np.pi/(f*3600)
      Tcordict=np.append(Tcordict,T)
  for k in tidedict.keys():
      Ttiddict=np.append(Ttiddict,tidedict[k])
  freqvec=np.append(Ttiddict,Tcordict)
  fig=plt.figure(figsize=(21,9))
  ax1=plt.subplot2grid((2,1), (1,0), rowspan=1, colspan=1)
  l1=ax1.plot(tvec,xb[1,:],color='r',label="Obs")
  l2=ax1.plot(tvec,xs[1,:],color='g',label="Mod") 
  handles, labels = ax1.get_legend_handles_labels()
  ax1.legend(handles,labels,loc=1,fontsize='18')
  ax1.set_xlim([5,23])
  for xc in freqvec:
      plt.axvline(x=xc,color='k', linestyle='--',alpha=0.5) 
  ak = lineid_plot.initial_annotate_kwargs()
  tide=np.array(list(s['tidedict'].keys()));cor=np.array(list(s['cordict'].values()))
  tidlab=np.append(tide,cor)
  lineid_plot.plot_line_ids(tvec, xb[1,:], freqvec, tidlab,ax=ax1,max_iter=30000,label1_size=18)
  plt.xlabel('Period [h]',fontsize=18,fontweight='bold')
  plt.ylabel('Phase (deg)',fontsize=18,fontweight='bold')
  ax2=fig.add_axes([0.125, 0.57, 0.775, 0.34],sharex=ax1)
  l1=ax2.plot(tvec,xb[0,:],color='r',label="Obs")
  l2=ax2.plot(tvec,xs[0,:],color='g',label="Mod")
  handles, labels = ax2.get_legend_handles_labels()
  ax2.legend(handles,labels,loc=1,fontsize='18')
  for xc in freqvec:
      plt.axvline(x=xc,color='k', linestyle='--',alpha=0.5) 
  lineid_plot.plot_line_ids(tvec, xb[0,:], freqvec, tidlab,ax=ax2,max_iter=30000,label1_size=18) #
  plt.ylabel('Amplitude ('+units+')',fontsize=18,fontweight='bold')
  plt.setp(ax2.get_xticklabels(), visible=False)
  plt.xticks(fontsize=18);plt.yticks(fontsize=18)
  plt.savefig(path+'/'+name+'.jpg',format='jpg')
  plt.show()
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


def pltvelquiver(loc,name,T,spindex,ylim,Usres,Vsres,Ubres,Vbres,Coru,Corv,labelname):
  llim=0;ulim=spindex
  colorswitch={'Obs':'red','Tides':'blue','Mod':'green','stidmodel':'black'}
  Y=np.zeros(len(Usres))
  fig=plt.figure(figsize=(15,4),frameon=True)
  ax1=plt.subplot(2,1,1)
  Tib=T[llim:ulim]
  ax1.quiver(Tib,Y[llim:ulim],Usres[llim:ulim],Vsres[llim:ulim],angles='uv',scale=1,scale_units='y',width=0.0004,color=colorswitch.get(labelname[0],"Invalid label name"),label=labelname[0])
  ax1.quiver(Tib,Y[llim:ulim],Ubres[llim:ulim],Vbres[llim:ulim],angles='uv',scale=1,scale_units='y',width=0.0004,color=colorswitch.get(labelname[1],"Invalid label name"),label=labelname[1])
  ax1.set_xlabel('Time',labelpad=0.05,fontsize=4, fontweight='bold')
  ax1.set_ylabel('Residual vel. vect.',labelpad=0.05,fontsize=4, fontweight='bold')
  ax1.text(0.995, 0.83, 'CorU: '+str(round(Coru[0],2)),verticalalignment='top', horizontalalignment='right',transform=ax1.transAxes,color='black', fontsize=4)
  ax1.text(0.995, 0.78, 'CorV: '+str(round(Corv[0],2)),verticalalignment='top', horizontalalignment='right',transform=ax1.transAxes,color='black', fontsize=4)
  labels=Tib[::96]
  plt.xticks(labels, labels, rotation='0',fontsize=4)
  plt.yticks(fontsize=4)
  plt.ylim([-ylim,ylim])
  ax1.legend(fontsize=4)
  llim=spindex;ulim=-1
  Tib=T[llim:ulim]
  ax2=plt.subplot(2,1,2)
  ax2.quiver(Tib,Y[llim:ulim],Usres[llim:ulim],Vsres[llim:ulim],angles='uv',scale=1,scale_units='y',width=0.0004,color=colorswitch.get(labelname[0],"Invalid label name"),label=labelname[0])
  ax2.quiver(Tib,Y[llim:ulim],Ubres[llim:ulim],Vbres[llim:ulim],angles='uv',scale=1,scale_units='y',width=0.0004,color=colorswitch.get(labelname[1],"Invalid label name"),label=labelname[1])
  ax2.set_xlabel('Time',labelpad=0.05,fontsize=5, fontweight='bold')
  ax2.set_ylabel('Residual velocity',labelpad=0.05,fontsize=5, fontweight='bold')
  ax2.legend(fontsize=4)
  labels=Tib[::96]
  plt.xticks(labels, labels, rotation='0',fontsize=4)
  plt.yticks(fontsize=4)
  plt.ylim([-ylim,ylim])
  plt.savefig(loc+name+'tidal_vel_vect.jpg',format='jpg',dpi=600)
  plt.close(fig)

def pltallbuoytracks(loc,BD):
  Buoys=['02','03','09','07','12','13','14','16']
  fig=plt.figure(figsize=(12, 12), frameon=True)
  ax=plt.axes(projection=ccrs.LambertAzimuthalEqualArea(central_longitude=25.0,central_latitude=77.0)) 
  ax.set_extent([16,28,74.5,80])
  colors=['red','green','magenta','darkblue','lime','orange','yellow','olive']
  fig.canvas.draw()
  xticks = [ 0, 4, 12, 16, 20, 24, 28, 32, 36]
  yticks = [72, 74.5, 77, 79.5, 82]
  ax.gridlines(xlocs=xticks, ylocs=yticks)
  # Label the end-points of the gridlines using the custom tick makers:
  ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
  ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
  lambert_xticks(ax, xticks)
  lambert_yticks(ax, yticks)
  i=0;pt=96*1
  for b in Buoys:
      # Plotting b-uoy obs and buoy simulated locations
      Xib=BD[b+'_x'];Yib=BD[b+'_y']
      plt.scatter(Xib[pt],Yib[pt],color='black',transform=ccrs.PlateCarree())
      plt.text(Xib[pt],Yib[pt], b, transform=ccrs.PlateCarree(),fontsize=15,fontweight='bold')
      plt.plot(Xib[pt:],Yib[pt:],color=colors[i],transform=ccrs.PlateCarree(),label=b)
      i+=1
  # gebco wms background 
  service='https://www.gebco.net/data_and_products/gebco_web_services/web_map_service/mapserv?'
  ax.add_wms(service,layers=['GEBCO_LATEST'],wms_kwargs={'width':900*2,'height':600*2})
  feature=cpf.GSHHSFeature(scale='i',levels=[1],facecolor='#e6e1e1',alpha=1)
  ax.add_feature(feature) 
  plt.savefig(loc+'/allbuoytracks.jpg',format='jpg',dpi=600)
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



## old function
# def pltvelquiver(loc,bname,Tib,Usres,Vsres,Ubres,Vbres):

#     def splt(llim,ulim,ax):
#         Y=np.zeros(len(Usres))
#         Tib1=Tib[llim:ulim]
#         ax.quiver(Tib1,Y[llim:ulim],Usres[llim:ulim],Vsres[llim:ulim],angles='uv',scale=1,scale_units='y',width=0.002,color='blue',label='sim')
#         ax.quiver(Tib1,Y[llim:ulim],Ubres[llim:ulim],Vbres[llim:ulim],angles='uv',scale=1,scale_units='y',width=0.002,color='red',label='obs')
#         # plt.quiver(Tib1,Y[llim:ulim],Ut[llim:ulim],Vt[llim:ulim],scale=1,scale_units='inches',width=0.002,color='green',label='tidal_cur')
#         # plt.quiver(Tib1,-Pgyt[llim:ulim],+Pgxt[llim:ulim],scale=2,scale_units='inches',color='magenta',label='-pgtid')
#         labels=Tib1[::96] 
#         plt.xticks(labels, labels, rotation='7.5',fontsize=4)
#         plt.ylim([-1.3,1.3])

#     m=2;n=7
#     fig=plt.figure(figsize=(15,4),frameon=True) 
#     llim=0
#     for i in range(m*n):
#         j=i+1
#         if j==1:       
#             ulim=llim+288
#             ax0=plt.subplot(m,n,j)
#             splt(llim,ulim,ax0)  
#             ax0.legend(fontsize=4)
#             ax0.set_xlabel('Time',labelpad=0.05,fontsize=4, fontweight='bold')
#             ax0.set_ylabel('Residual vel. vect.',labelpad=0.05,fontsize=4, fontweight='bold')
#             plt.yticks(fontsize=4)
#             llim=ulim
#         else:
#             ulim=llim+288
#             ax=plt.subplot(m,n,j,sharey=ax0)
#             splt(llim,ulim,ax)  
#             plt.setp(ax.get_yticklabels(), visible=False)
#             yticks = ax.yaxis.get_major_ticks()
#             yticks[-1].label1.set_visible(False)
#             llim=ulim      
#     plt.subplots_adjust(wspace=.0)
#     plt.savefig(loc+bname+'tidal_vel_vect.jpg',format='jpg',dpi=600)
#     plt.close(fig)

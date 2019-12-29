%script for interpolating gtsm 5km data to a coarser grid.

% reading nc file 

file='../Data_from_models/gtsm_truncated.nc';
%regular grid for averaging (spherical coordinates)
xmin=0;
xmax=40;
xstep=400;
ymin=70;
ymax=82;
ystep=120;
[X,Y]=meshgrid(linspace(xmin,xmax,xstep),linspace(ymin,ymax,ystep));

Xt=ncread(file,'/GTSM Tidal Velocity truncated data/FlowElem_xcc');
Yt=ncread(file,'/GTSM Tidal Velocity truncated data/FlowElem_ycc');
Ut=ncread(file,'/GTSM Tidal Velocity truncated data/ucx');
Vt=ncread(file,'/GTSM Tidal Velocity truncated data/ucy');
Tt=ncread(file,'/GTSM Tidal Velocity truncated data/Time');
SI=scatteredInterpolant(Xt,Yt,Ut(:,2));
U=SI(X,Y);
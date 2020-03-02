from __future__ import division
import cf, cfplot as cfp
import numpy as np
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap

#Output path
path = '/gws/nopw/j04/klingaman/amulya/scripts/s2s/coupling/output/'

#Ocean mask
mask = Dataset('/gws/nopw/j04/klingaman/amulya/data/surface_s2s_skill/soilm/dubstep_output/UKMO/SA_UKMO_ensemble_sm20_weekly_hindcasts_lead1-5_200001.nc', 'r').variables['week_sm20'][:][0,0,0,:,:]
mask[mask<0] = np.nan
mask[mask>=0] = 1.0

#Land mask
oc = Dataset('/gws/nopw/j04/klingaman/amulya/data/surface_s2s_skill/soilm/dubstep_output/UKMO/SA_UKMO_ensemble_sm20_weekly_hindcasts_lead1-5_200001.nc', 'r').variables['week_sm20'][:][0,0,0,:,:]
oc[oc>=0] = np.nan
oc[oc<0] = 9000

#Load numpy output for each S2S data and all of the coupling strength indices
a = np.load('../output/era5.npz')
era5_gam = a['gam']
era5_tci = a['tci']
era5_tet = a['tet']
era5_pleg = a['pleg']
era5_tleg = a['tleg']

a = np.load('../output/ukmo.npz')
ukmo_gam = a['gam']
ukmo_tci = a['tci']
ukmo_tet = a['tet']
ukmo_pleg = a['pleg']
ukmo_tleg = a['tleg']

a = np.load('../output/ncep.npz')
ncep_gam = a['gam']
ncep_tci = a['tci']
ncep_tet = a['tet']
ncep_pleg = a['pleg']
ncep_tleg = a['tleg']

a = np.load('../output/ecmwf.npz')
ecmwf_gam = a['gam']
ecmwf_tci = a['tci']
ecmwf_tet = a['tet']
ecmwf_pleg = a['pleg']
ecmwf_tleg = a['tleg']

#Stack the coupling strength indices for different S2S data together
gam = np.stack((era5_gam, ukmo_gam, ncep_gam, ecmwf_gam))
tci = np.stack((era5_tci, ukmo_tci, ncep_tci, ecmwf_tci))
tet = np.stack((era5_tet, ukmo_tet, ncep_tet, ecmwf_tet))
pleg = np.stack((era5_pleg, ukmo_pleg, ncep_pleg, ecmwf_pleg))
tleg = np.stack((era5_tleg, ukmo_tleg, ncep_tleg, ecmwf_tleg))

lon = np.arange(lonlim[0],lonlim[1]+0.5,1.5) #longitude array
lat = np.arange(latlim[0],latlim[1]+0.5,1.5) #latitde array
lw = 1 #linewidth
gl = 20 #lat lon skips for drawing on the map
latlim = [-50,15] #latitude limits
lonlim = [-90,-30] #longitude limits
data = ['ERA5','UKMO','NCEP','ECMWF'] #dataset titles

#Zeng's gamma plot
cmap = plt.get_cmap('RdBu'); cmap.set_bad(color = '0.75', alpha = 1.)
cmin = -0.2
cmax = 0.2
cspc = 0.01
cstp = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
label = 'Gamma'
title = 'PR-ET'
units = ''
clabel = label+' '+title+' '+units
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
fig = plt.figure(figsize=(12,9))
ax=plt.subplot(4,5,3)
mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
mymap.drawcoastlines(linewidth=lw)
mymap.drawcountries(linewidth=lw)
mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
x, y = mymap(*np.meshgrid(lon,lat))
uncal = mymap.pcolormesh(x,y,gam[0,0,:,:]*mask,vmin=cmin,vmax=cmax,cmap=cmap,norm=norm)
mymap.pcolormesh(x,y,oc,vmax=cmax,cmap=plt.get_cmap('Greys_r'),norm=norm)
plt.title(str(data[0]), fontsize=10)
for d in range(1,4):
  for w in range(0,5):
    ax=plt.subplot(4,5,(d*5)+w+1)
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    mymap.drawcountries(linewidth=lw)
    if w ==0:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
    else:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    if d == 3:
      mymap.drawmeridians(np.arange(0,360,gl),labels=[1,0,0,1],labelstyle='+/-')
    else:
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(lon,lat))
    uncal = mymap.pcolormesh(x,y,gam[d,w,:,:]*mask,vmin=cmin,vmax=cmax,cmap=cmap,norm=norm)
    mymap.pcolormesh(x,y,oc,vmax=cmax,cmap=plt.get_cmap('Greys_r'),norm=norm)
    plt.title(str(data[d])+' Week-'+str(w+1), fontsize=10)

cbaxes = fig.add_axes([0.92, 0.11, 0.015, 0.77]) 
cb = plt.colorbar(uncal, cax=cbaxes, label=clabel, extend='both', ticks=np.arange(cmin,cmax+cstp,cstp), orientation="vertical")
plt.savefig(path+label+'_'+title+'.png', bbox_inches='tight')
#plt.show()


#TCI plot
cmap = plt.get_cmap('RdBu'); cmap.set_bad(color = '0.75', alpha = 1.)
cmin = -1
cmax = 1
cspc = 0.1
cstp = 0.5
clevs = np.arange(cmin,cmax+cspc,cspc)
label = 'TCI'
title = 'SM-ET'
units = '(mm day$^{-1}$)'
clabel = label+' '+title+' '+units
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
fig = plt.figure(figsize=(12,9))
ax=plt.subplot(4,5,3)
mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
mymap.drawcoastlines(linewidth=lw)
mymap.drawcountries(linewidth=lw)
mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
x, y = mymap(*np.meshgrid(lon,lat))
uncal = mymap.pcolormesh(x,y,tci[0,0,:,:]*mask,vmin=cmin,vmax=cmax,cmap=cmap,norm=norm)
mymap.pcolormesh(x,y,oc,vmax=cmax,cmap=plt.get_cmap('Greys_r'),norm=norm)
plt.title(str(data[0]), fontsize=10)
for d in range(1,4):
  for w in range(0,5):
    ax=plt.subplot(4,5,(d*5)+w+1)
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    mymap.drawcountries(linewidth=lw)
    if w ==0:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
    else:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    if d == 3:
      mymap.drawmeridians(np.arange(0,360,gl),labels=[1,0,0,1],labelstyle='+/-')
    else:
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(lon,lat))
    uncal = mymap.pcolormesh(x,y,tci[d,w,:,:]*mask,vmin=cmin,vmax=cmax,cmap=cmap,norm=norm)
    mymap.pcolormesh(x,y,oc,vmax=cmax,cmap=plt.get_cmap('Greys_r'),norm=norm)
    plt.title(str(data[d])+' Week-'+str(w+1), fontsize=10)

cbaxes = fig.add_axes([0.92, 0.11, 0.015, 0.77]) 
cb = plt.colorbar(uncal, cax=cbaxes, label=clabel, extend='both', ticks=np.arange(cmin,cmax+cstp,cstp), orientation="vertical")
plt.savefig(path+label+'_'+title+'.png', bbox_inches='tight')
#plt.show()


#TET plot
cmap = plt.get_cmap('RdBu'); cmap.set_bad(color = '0.75', alpha = 1.)
cmin = -1
cmax = 1
cspc = 0.1
cstp = 0.5
clevs = np.arange(cmin,cmax+cspc,cspc)
label = 'TET'
title = 'T2-ET'
units = ''
clabel = label+' '+title+' '+units
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
fig = plt.figure(figsize=(12,9))
ax=plt.subplot(4,5,3)
mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
mymap.drawcoastlines(linewidth=lw)
mymap.drawcountries(linewidth=lw)
mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
x, y = mymap(*np.meshgrid(lon,lat))
uncal = mymap.pcolormesh(x,y,tet[0,0,:,:]*mask,vmin=cmin,vmax=cmax,cmap=cmap,norm=norm)
mymap.pcolormesh(x,y,oc,vmax=cmax,cmap=plt.get_cmap('Greys_r'),norm=norm)
plt.title(str(data[0]), fontsize=10)
for d in range(1,4):
  for w in range(0,5):
    ax=plt.subplot(4,5,(d*5)+w+1)
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    mymap.drawcountries(linewidth=lw)
    if w ==0:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
    else:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    if d == 3:
      mymap.drawmeridians(np.arange(0,360,gl),labels=[1,0,0,1],labelstyle='+/-')
    else:
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(lon,lat))
    uncal = mymap.pcolormesh(x,y,tet[d,w,:,:],vmin=cmin,vmax=cmax,cmap=cmap,norm=norm)
    mymap.pcolormesh(x,y,oc,vmax=cmax,cmap=plt.get_cmap('Greys_r'),norm=norm)
    plt.title(str(data[d])+' Week-'+str(w+1), fontsize=10)

cbaxes = fig.add_axes([0.92, 0.11, 0.015, 0.77]) 
cb = plt.colorbar(uncal, cax=cbaxes, label=clabel, extend='both', ticks=np.arange(cmin,cmax+cstp,cstp), orientation="vertical")
plt.savefig(path+label+'_'+title+'.png', bbox_inches='tight')
#plt.show()


#Two-legged metrics

#SM-ET-PRECIP pathway plot
cmap = plt.get_cmap('RdBu'); cmap.set_bad(color = '0.75', alpha = 1.)
cmin = [-1,-4,-2]
cmax = [1,4,2]
cspc = 0.1
cstp = 0.5
label = ['Surface-Leg', 'Atmosphere-Leg', 'Total-Feedback']
title = ['SM-ET', 'ET-PR', 'SM-PR']
units = ['','','']
for i in range(3):
  clabel = label[i]+' '+title[i]+' '+units[i]
  clevs = np.arange(cmin[i],cmax[i]+cspc,cspc)
  norm = BoundaryNorm(boundaries=clevs, ncolors=256)
  fig = plt.figure(figsize=(12,9))
  ax=plt.subplot(4,5,3)
  mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(lon,lat))
  uncal = mymap.pcolormesh(x,y,pleg[0,i,0,:,:]*mask,vmin=cmin[i],vmax=cmax[i],cmap=cmap,norm=norm)
  mymap.pcolormesh(x,y,oc,vmax=cmax[i],cmap=plt.get_cmap('Greys_r'),norm=norm)
  plt.title(str(data[0]), fontsize=10)
  for d in range(1,4):
    for w in range(0,5):
      ax=plt.subplot(4,5,(d*5)+w+1)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      mymap.drawcountries(linewidth=lw)
      if w ==0:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
      else:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      if d == 3:
        mymap.drawmeridians(np.arange(0,360,gl),labels=[1,0,0,1],labelstyle='+/-')
      else:
        mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(lon,lat))
      uncal = mymap.pcolormesh(x,y,pleg[d,i,w,:,:]*mask,vmin=cmin[i],vmax=cmax[i],cmap=cmap,norm=norm)
      mymap.pcolormesh(x,y,oc,vmax=cmax[i],cmap=plt.get_cmap('Greys_r'),norm=norm)
      plt.title(str(data[d])+' Week-'+str(w+1), fontsize=10)
  cbaxes = fig.add_axes([0.92, 0.11, 0.015, 0.77]) 
  cb = plt.colorbar(uncal, cax=cbaxes, label=clabel, extend='both', ticks=np.arange(cmin[i],cmax[i]+cstp,cstp), orientation="vertical")
  plt.savefig(path+'P_'+label[i]+'_'+title[i]+'.png', bbox_inches='tight')
  #plt.show()

#SM-ET-T2M pathway plot
cmap = plt.get_cmap('RdBu'); cmap.set_bad(color = '0.75', alpha = 1.)
cmin = [-1,-2,-1]
cmax = [1,2,1]
cspc = 0.1
cstp = 0.5
label = ['Surface-Leg', 'Atmosphere-Leg', 'Total-Feedback']
title = ['SM-ET', 'ET-T2', 'SM-T2']
units = ['','','']
for i in range(3):
  clabel = label[i]+' '+title[i]+' '+units[i]
  clevs = np.arange(cmin[i],cmax[i]+cspc,cspc)
  norm = BoundaryNorm(boundaries=clevs, ncolors=256)
  fig = plt.figure(figsize=(12,9))
  ax=plt.subplot(4,5,3)
  mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(lon,lat))
  uncal = mymap.pcolormesh(x,y,tleg[0,i,0,:,:]*mask,vmin=cmin[i],vmax=cmax[i],cmap=cmap,norm=norm)
  mymap.pcolormesh(x,y,oc,vmax=cmax[i],cmap=plt.get_cmap('Greys_r'),norm=norm)
  plt.title(str(data[0]), fontsize=10)
  for d in range(1,4):
    for w in range(0,5):
      ax=plt.subplot(4,5,(d*5)+w+1)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
      mymap.drawcoastlines(linewidth=lw)
      mymap.drawcountries(linewidth=lw)
      if w ==0:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
      else:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      if d == 3:
        mymap.drawmeridians(np.arange(0,360,gl),labels=[1,0,0,1],labelstyle='+/-')
      else:
        mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(lon,lat))
      uncal = mymap.pcolormesh(x,y,tleg[d,i,w,:,:]*mask,vmin=cmin[i],vmax=cmax[i],cmap=cmap,norm=norm)
      mymap.pcolormesh(x,y,oc,vmax=cmax[i],cmap=plt.get_cmap('Greys_r'),norm=norm)
      plt.title(str(data[d])+' Week-'+str(w+1), fontsize=10)
  cbaxes = fig.add_axes([0.92, 0.11, 0.015, 0.77]) 
  cb = plt.colorbar(uncal, cax=cbaxes, label=clabel, extend='both', ticks=np.arange(cmin[i],cmax[i]+cstp,cstp), orientation="vertical")
  plt.savefig(path+'T_'+label[i]+'_'+title[i]+'.png', bbox_inches='tight')
  #plt.show()


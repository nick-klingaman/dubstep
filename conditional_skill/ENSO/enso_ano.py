from __future__ import division
#import cf, cfplot as cfp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap

chirps = np.load('enso_chirps.npz')
pel_chirps = chirps['apel_chirps']
pne_chirps = chirps['apne_chirps']
pla_chirps = chirps['apla_chirps']
pno_chirps = chirps['apno_chirps']
chirps = np.stack((np.nanmean(pel_chirps, axis=0),np.nanmean(pne_chirps, axis=0),np.nanmean(pla_chirps, axis=0),np.nanmean(pno_chirps, axis=0)))

bam = np.load('enso_bam.npz')
pel_bam = bam['apel_bam']
pne_bam = bam['apne_bam']
pla_bam = bam['apla_bam']
pno_bam = bam['apno_bam']
bam = np.stack((np.nanmean(pel_bam, axis=0),np.nanmean(pne_bam, axis=0),np.nanmean(pla_bam, axis=0),np.nanmean(pno_bam, axis=0)))

ecmf = np.load('enso_ecmf.npz')
pel_ecmf = ecmf['apel_ecmf']
pne_ecmf = ecmf['apne_ecmf']
pla_ecmf = ecmf['apla_ecmf']
pno_ecmf = ecmf['apno_ecmf']
ecmf = np.stack((np.nanmean(pel_ecmf, axis=0),np.nanmean(pne_ecmf, axis=0),np.nanmean(pla_ecmf, axis=0),np.nanmean(pno_ecmf, axis=0)))

ncep = np.load('enso_ncep.npz')
pel_ncep = ncep['apel_ncep']
pne_ncep = ncep['apne_ncep']
pla_ncep = ncep['apla_ncep']
pno_ncep = ncep['apno_ncep']
ncep = np.stack((np.nanmean(pel_ncep, axis=0),np.nanmean(pne_ncep, axis=0),np.nanmean(pla_ncep, axis=0),np.nanmean(pno_ncep, axis=0)))

ukmo = np.load('enso_ukmo.npz')
pel_ukmo = ukmo['apel_ukmo']
pne_ukmo = ukmo['apne_ukmo']
pla_ukmo = ukmo['apla_ukmo']
pno_ukmo = ukmo['apno_ukmo']
ukmo = np.stack((np.nanmean(pel_ukmo, axis=0),np.nanmean(pne_ukmo, axis=0),np.nanmean(pla_ukmo, axis=0),np.nanmean(pno_ukmo, axis=0)))


chirps_ano = chirps.copy()
bam_ano = (np.nanmean(bam, axis=2))
ecmf_ano = (np.nanmean(ecmf, axis=2))
ncep_ano = (np.nanmean(ncep, axis=2)) 
ukmo_ano = (np.nanmean(ukmo, axis=2))
ano = np.stack((chirps_ano, bam_ano, ecmf_ano, ncep_ano, ukmo_ano))

latlim = [-49.5,19.5]
lonlim = [-90,-20]
lon = np.arange(lonlim[0],lonlim[1]+0.5,1.5)
lat = np.arange(latlim[0],latlim[1]+0.5,1.5)
cmap = matplotlib.cm.BrBG
cmin = -2
cmax = 2
cspc = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel = 'ANO'
dlabel = 'Precipitation anomaly (mm day$^{-1}$)'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
lw = 1
gl = 20
latext = [-50,15]
lonext = [-90,-30]
data = ['CHIRPS','BAM','ECMWF','NCEP','UKMO']
phase = ['ElNino','Neutral','LaNina','ENSO']
WEEK=5

mask = chirps[0,0,:,:].squeeze()
mask[~np.isnan(mask)] = 1.0

for d in range(len(data)):
    print data[d]
    for p in range(len(phase)):
        for w in range(WEEK):
            fig = plt.figure(figsize=(3,3))
            mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latext[0],urcrnrlat=latext[1],llcrnrlon=lonext[0],urcrnrlon=lonext[1])
	    mymap.drawcoastlines(linewidth=lw)
            mymap.drawcountries(linewidth=lw)
            mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
            mymap.drawmeridians(np.arange(0,360,gl),labels=[1,0,0,1],labelstyle='+/-')
            x, y = mymap(*np.meshgrid(lon,lat))
            uncal = mymap.pcolormesh(x,y,ano[d,p,w,:,:],vmin=cmin,vmax=cmax,cmap=cmap,norm=norm)
#            plt.title(str(data[d])+'_'+str(phase[p])+'_'+'Week'+str(w+1)+'_'+str(clabel), fontsize=15)
            plt.savefig(str(data[d])+'_'+str(phase[p])+'_'+'Week'+str(w+1)+'_'+str(clabel)+'.eps')
#            plt.show()
            plt.close()

fig = plt.figure(figsize=(4,1))
ax = fig.add_axes([0.05, 0.5, 0.9, 0.07])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='both',spacing='uniform',orientation='horizontal', label=dlabel)
#cb.set_label(clabel)
plt.savefig(str(clabel)+'_horizontal_colorbar.eps')
#plt.show()

fig = plt.figure(figsize=(1,4))
ax = fig.add_axes([0.15, 0.05, 0.07, 0.9])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='both',spacing='uniform',orientation='vertical',label=dlabel)
#cb.set_label(clabel)
plt.savefig(str(clabel)+'_vertical_colorbar.eps')
#plt.show()


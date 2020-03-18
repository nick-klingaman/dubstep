from __future__ import division
#import cf, cfplot as cfp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap

mjo = np.load('../mjo_bchirps.npz')
bchirps = mjo['mjoa_bchirps']

mjo = np.load('../mjo_chirps.npz')
chirps = mjo['mjoa_chirps']

mjo = np.load('../mjo_ukmo.npz')
ukmo = mjo['mjoa_ukmo']

mjo = np.load('../mjo_ncep.npz')
ncep = mjo['mjoa_ncep']

mjo = np.load('../mjo_ecmf.npz')
ecmf = mjo['mjoa_ecmf']

mjo = np.load('../mjo_bam.npz')
bam = mjo['mjoa_bam']

bchirps_ano = [0]*4;chirps_ano = [0]*4;ukmo_ano = [0]*4;ncep_ano = [0]*4;ecmf_ano = [0]*4;bam_ano = [0]*4
idx1 = [8,2,4,6]
idx2 = [1,3,5,7]
for i in range(len(idx1)):
    bchirps_ano[i] = np.nanmean(np.vstack((bchirps[idx1[i]],bchirps[idx2[i]])), axis=0) 
    chirps_ano[i] = np.nanmean(np.vstack((chirps[idx1[i]],chirps[idx2[i]])), axis=0)
    ukmo_ano[i] = np.nanmean(np.vstack((ukmo[idx1[i]],ukmo[idx2[i]])), axis=(0,2))  
    ncep_ano[i] = np.nanmean(np.vstack((ncep[idx1[i]],ncep[idx2[i]])), axis=(0,2))
    ecmf_ano[i] = np.nanmean(np.vstack((ecmf[idx1[i]],ecmf[idx2[i]])), axis=(0,2))
    bam_ano[i] = np.nanmean(np.vstack((bam[idx1[i]],bam[idx2[i]])), axis=(0,2))

ano = [chirps_ano, bam_ano, ecmf_ano, ncep_ano, ukmo_ano]

latlim = [-49.5,19.5]
lonlim = [-90,-20]
lon = np.arange(lonlim[0],lonlim[1]+0.5,1.5)
lat = np.arange(latlim[0],latlim[1]+0.5,1.5)
cmap = matplotlib.cm.BrBG
cmin = -3
cmax = 3
cspc = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel = 'ANO'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
lw = 1
gl = 20
latext = [-50,15]
lonext = [-90,-30]
data = ['CHIRPS','BAM','ECMWF','NCEP','UKMO']
phase = ['MJO-8+1','MJO-2+3','MJO-4+5','MJO-6+7']
WEEK=5

for d in range(len(data)):
    for p in range(len(phase)):
        for w in range(WEEK):
            fig = plt.figure(figsize=(3,3))
            mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latext[0],urcrnrlat=latext[1],llcrnrlon=lonext[0],urcrnrlon=lonext[1])
	    mymap.drawcoastlines(linewidth=lw)
            mymap.drawcountries(linewidth=lw)
            mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,1],labelstyle='+/-')
            mymap.drawmeridians(np.arange(0,360,gl),labels=[1,0,0,1],labelstyle='+/-')
            x, y = mymap(*np.meshgrid(lon,lat))
            uncal = mymap.pcolormesh(x,y,ano[d][p][w,:,:],vmin=cmin,vmax=cmax,cmap=cmap,norm=norm)
#            plt.title(str(data[d])+'_'+str(phase[p])+'_'+'Week'+str(w+1)+'_'+str(clabel), fontsize=15)
            plt.savefig(str(data[d])+'_'+str(phase[p])+'_'+'Week'+str(w+1)+'_'+str(clabel)+'.png')
#            plt.show()
            plt.close()

fig = plt.figure(figsize=(4,1))
ax = fig.add_axes([0.05, 0.5, 0.9, 0.07])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='both',spacing='uniform',orientation='horizontal')
#cb.set_label(clabel)
plt.savefig(str(clabel)+'_horizontal_colorbar.png',bbox_inches='tight')
#plt.show()

fig = plt.figure(figsize=(1,4))
ax = fig.add_axes([0.4, 0.05, 0.07, 0.9])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='both',spacing='uniform',orientation='vertical')
#cb.set_label(clabel)
plt.savefig(str(clabel)+'_vertical_colorbar.png',bbox_inches='tight')
#plt.show()


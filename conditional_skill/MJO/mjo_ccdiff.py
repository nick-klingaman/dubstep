from __future__ import division
#import cf, cfplot as cfp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
import scipy as sp
from scipy import stats

LAT=47
LON=47
WEEK=5

mjo = np.load('mjo_bchirps.npz')
bchirps = mjo['mjoa_bchirps']

mjo = np.load('mjo_chirps.npz')
chirps = mjo['mjoa_chirps']

mjo = np.load('mjo_ukmo.npz')
ukmo = mjo['mjoa_ukmo']

mjo = np.load('mjo_ncep.npz')
ncep = mjo['mjoa_ncep']

mjo = np.load('mjo_ecmwf.npz')
ecmwf = mjo['mjoa_ecmwf']

mjo = np.load('mjo_bam.npz')
bam = mjo['mjoa_bam']

bchirps_ano = [0]*5;chirps_ano = [0]*5;ukmo_ano = [0]*5;ncep_ano = [0]*5;ecmwf_ano = [0]*5;bam_ano = [0]*5
idx1 = [8,2,4,6]
idx2 = [1,3,5,7]
for i in range(len(idx1)):
    bchirps_ano[i] = np.vstack((bchirps[idx1[i]],bchirps[idx2[i]]))
    chirps_ano[i] = np.vstack((chirps[idx1[i]],chirps[idx2[i]]))
    ukmo_ano[i] = np.nanmean(np.vstack((ukmo[idx1[i]],ukmo[idx2[i]])), axis=2)  
    ncep_ano[i] = np.nanmean(np.vstack((ncep[idx1[i]],ncep[idx2[i]])), axis=2)
    ecmwf_ano[i] = np.nanmean(np.vstack((ecmwf[idx1[i]],ecmwf[idx2[i]])), axis=2)
    bam_ano[i] = np.nanmean(np.vstack((bam[idx1[i]],bam[idx2[i]])), axis=2)

bchirps_ano[4] = bchirps[9].copy()
chirps_ano[4] = chirps[9].copy()
ukmo_ano[4] = np.nanmean(ukmo[9], axis=2)  
ncep_ano[4] = np.nanmean(ncep[9], axis=2)  
ecmwf_ano[4] = np.nanmean(ecmwf[9], axis=2)  
bam_ano[4] = np.nanmean(bam[9], axis=2)  

bam_acc = [0]*5;bam_pval = [0]*5
for n in range(5):
    bam_acc[n] = np.zeros((WEEK,LAT,LON))
    bam_pval[n] = np.zeros((WEEK,LAT,LON))
    for i in range(0,LAT):
        for j in range(0,LON):
            for w in range(0,WEEK):
    	        bam_acc[n][w,i,j] = sp.stats.pearsonr(bam_ano[n][:,w,i,j], bchirps_ano[n][:,w,i,j])[0]
    	        bam_pval[n][w,i,j] = sp.stats.pearsonr(bam_ano[n][:,w,i,j], bchirps_ano[n][:,w,i,j])[1]

ecmwf_acc = [0]*5;ecmwf_pval = [0]*5
for n in range(5):
    ecmwf_acc[n] = np.zeros((WEEK,LAT,LON))
    ecmwf_pval[n] = np.zeros((WEEK,LAT,LON))
    for i in range(0,LAT):
        for j in range(0,LON):
            for w in range(0,WEEK):
    	        ecmwf_acc[n][w,i,j] = sp.stats.pearsonr(ecmwf_ano[n][:,w,i,j], chirps_ano[n][:,w,i,j])[0]
    	        ecmwf_pval[n][w,i,j] = sp.stats.pearsonr(ecmwf_ano[n][:,w,i,j], chirps_ano[n][:,w,i,j])[1]

ncep_acc = [0]*5;ncep_pval = [0]*5
for n in range(5):
    ncep_acc[n] = np.zeros((WEEK,LAT,LON))
    ncep_pval[n] = np.zeros((WEEK,LAT,LON))
    for i in range(0,LAT):
        for j in range(0,LON):
            for w in range(0,WEEK):
    	        ncep_acc[n][w,i,j] = sp.stats.pearsonr(ncep_ano[n][:,w,i,j], chirps_ano[n][:,w,i,j])[0]
    	        ncep_pval[n][w,i,j] = sp.stats.pearsonr(ncep_ano[n][:,w,i,j], chirps_ano[n][:,w,i,j])[1]

ukmo_acc = [0]*5;ukmo_pval = [0]*5
for n in range(5):
    ukmo_acc[n] = np.zeros((WEEK,LAT,LON))
    ukmo_pval[n] = np.zeros((WEEK,LAT,LON))
    for i in range(0,LAT):
        for j in range(0,LON):
            for w in range(0,WEEK):
    	        ukmo_acc[n][w,i,j] = sp.stats.pearsonr(ukmo_ano[n][:,w,i,j], chirps_ano[n][:,w,i,j])[0]
    	        ukmo_pval[n][w,i,j] = sp.stats.pearsonr(ukmo_ano[n][:,w,i,j], chirps_ano[n][:,w,i,j])[1]

acc = [bam_acc, ecmwf_acc, ncep_acc, ukmo_acc]
pval = [bam_pval, ecmwf_pval, ncep_pval, ukmo_pval]
#pval = np.ma.masked_less(pval, 0.05)
#pval[np.isnan(acc)]=np.nan

latlim = [-49.5,19.5]
lonlim = [-90,-20]
lon = np.arange(lonlim[0],lonlim[1]+0.5,1.5)
lat = np.arange(latlim[0],latlim[1]+0.5,1.5)
cmap = matplotlib.cm.RdYlBu
cmin = -0.5
cmax = 0.5
cspc = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel = 'ACC_Difference'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
lw = 1
gl = 20
latext = [-50,15]
lonext = [-90,-30]
data = ['BAM','ECMWF','NCEP','UKMO']
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
            uncal = mymap.pcolormesh(x,y,acc[d][p][w,:,:]-acc[d][4][w,:,:],vmin=cmin,vmax=cmax,cmap=cmap,norm=norm)
#	    mymap.pcolor(x, y, pval[d,p,w,:,:], hatch='.', alpha=0.)
#            plt.title(str(data[d])+'_'+str(phase[p])+'_'+'Week'+str(w+1)+'_'+str(clabel), fontsize=15)
#           fig.tight_layout()
            plt.savefig(str(data[d])+'_'+str(phase[p])+'_minus_Allyrs_'+'Week'+str(w+1)+'_'+str(clabel)+'.ps')
            #plt.show()

fig = plt.figure(figsize=(4,1))
ax = fig.add_axes([0.05, 0.5, 0.9, 0.07])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='both',spacing='uniform',orientation='horizontal')
#cb.set_label(clabel)
#fig.tight_layout()
plt.savefig(str(clabel)+'_horizontal_colorbar.ps')
#plt.show()

fig = plt.figure(figsize=(1,4))
ax = fig.add_axes([0.4, 0.05, 0.07, 0.9])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='both',spacing='uniform',orientation='vertical')
#cb.set_label(clabel)
#fig.tight_layout()
plt.savefig(str(clabel)+'_vertical_colorbar.ps')
#plt.show()

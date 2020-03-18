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
NUM=1000
NEL = 56
NNE = 108
NLA = 56
NNO = 112
NAL = 220
REG = 6

'''
chirps = np.load('enso_chirps.npz')
pel_chirps = chirps['apel_chirps']
pne_chirps = chirps['apne_chirps']
pla_chirps = chirps['apla_chirps']
pno_chirps = chirps['apno_chirps']
pal_chirps = chirps['apal_chirps']

bchirps = np.load('enso_bchirps.npz')
pel_bchirps = bchirps['apel_bchirps']
pne_bchirps = bchirps['apne_bchirps']
pla_bchirps = bchirps['apla_bchirps']
pno_bchirps = bchirps['apno_bchirps']
pal_bchirps = bchirps['apal_bchirps']

bam = np.load('enso_bam.npz')
pel_bam = np.nanmean(bam['apel_bam'], axis=2)
pne_bam = np.nanmean(bam['apne_bam'], axis=2)
pla_bam = np.nanmean(bam['apla_bam'], axis=2)
pno_bam = np.nanmean(bam['apno_bam'], axis=2)
pal_bam = np.nanmean(bam['apal_bam'], axis=2)

ecmf = np.load('enso_ecmf.npz')
pel_ecmf = np.nanmean(ecmf['apel_ecmf'], axis=2)
pne_ecmf = np.nanmean(ecmf['apne_ecmf'], axis=2)
pla_ecmf = np.nanmean(ecmf['apla_ecmf'], axis=2)
pno_ecmf = np.nanmean(ecmf['apno_ecmf'], axis=2)
pal_ecmf = np.nanmean(ecmf['apal_ecmf'], axis=2)

ncep = np.load('enso_ncep.npz')
pel_ncep = np.nanmean(ncep['apel_ncep'], axis=2)
pne_ncep = np.nanmean(ncep['apne_ncep'], axis=2)
pla_ncep = np.nanmean(ncep['apla_ncep'], axis=2)
pno_ncep = np.nanmean(ncep['apno_ncep'], axis=2)
pal_ncep = np.nanmean(ncep['apal_ncep'], axis=2)

ukmo = np.load('enso_ukmo.npz')
pel_ukmo = np.nanmean(ukmo['apel_ukmo'], axis=2)
pne_ukmo = np.nanmean(ukmo['apne_ukmo'], axis=2)
pla_ukmo = np.nanmean(ukmo['apla_ukmo'], axis=2)
pno_ukmo = np.nanmean(ukmo['apno_ukmo'], axis=2)
pal_ukmo = np.nanmean(ukmo['apal_ukmo'], axis=2)

pel_acc_ukmo = np.zeros((WEEK,LAT,LON)) ; pel_pval_ukmo = np.zeros((WEEK,LAT,LON))
pne_acc_ukmo = np.zeros((WEEK,LAT,LON)) ; pne_pval_ukmo = np.zeros((WEEK,LAT,LON))
pla_acc_ukmo = np.zeros((WEEK,LAT,LON)) ; pla_pval_ukmo = np.zeros((WEEK,LAT,LON))
pno_acc_ukmo = np.zeros((WEEK,LAT,LON)) ; pno_pval_ukmo = np.zeros((WEEK,LAT,LON))
pal_acc_ukmo = np.zeros((WEEK,LAT,LON)) ; pal_pval_ukmo = np.zeros((WEEK,LAT,LON))
for i in range(0,LAT):
    for j in range(0,LON):
        for w in range(0,WEEK):
	    pel_acc_ukmo[w,i,j] = sp.stats.pearsonr(pel_ukmo[:,w,i,j], pel_chirps[:,w,i,j])[0]
	    pel_pval_ukmo[w,i,j] = sp.stats.pearsonr(pel_ukmo[:,w,i,j], pel_chirps[:,w,i,j])[1]
	    pne_acc_ukmo[w,i,j] = sp.stats.pearsonr(pne_ukmo[:,w,i,j], pne_chirps[:,w,i,j])[0]
	    pne_pval_ukmo[w,i,j] = sp.stats.pearsonr(pne_ukmo[:,w,i,j], pne_chirps[:,w,i,j])[1]
	    pla_acc_ukmo[w,i,j] = sp.stats.pearsonr(pla_ukmo[:,w,i,j], pla_chirps[:,w,i,j])[0]
	    pla_pval_ukmo[w,i,j] = sp.stats.pearsonr(pla_ukmo[:,w,i,j], pla_chirps[:,w,i,j])[1]
	    pno_acc_ukmo[w,i,j] = sp.stats.pearsonr(pno_ukmo[:,w,i,j], pno_chirps[:,w,i,j])[0]
	    pno_pval_ukmo[w,i,j] = sp.stats.pearsonr(pno_ukmo[:,w,i,j], pno_chirps[:,w,i,j])[1]
	    pal_acc_ukmo[w,i,j] = sp.stats.pearsonr(pal_ukmo[:,w,i,j], pal_chirps[:,w,i,j])[0]
	    pal_pval_ukmo[w,i,j] = sp.stats.pearsonr(pal_ukmo[:,w,i,j], pal_chirps[:,w,i,j])[1]


ukmo_acc = np.stack((pel_acc_ukmo, pne_acc_ukmo, pla_acc_ukmo, pno_acc_ukmo, pal_acc_ukmo))
ukmo_pval = np.stack((pel_pval_ukmo, pne_pval_ukmo, pla_pval_ukmo, pno_pval_ukmo, pal_pval_ukmo))

pel_acc_ncep = np.zeros((WEEK,LAT,LON)) ; pel_pval_ncep = np.zeros((WEEK,LAT,LON))
pne_acc_ncep = np.zeros((WEEK,LAT,LON)) ; pne_pval_ncep = np.zeros((WEEK,LAT,LON))
pla_acc_ncep = np.zeros((WEEK,LAT,LON)) ; pla_pval_ncep = np.zeros((WEEK,LAT,LON))
pno_acc_ncep = np.zeros((WEEK,LAT,LON)) ; pno_pval_ncep = np.zeros((WEEK,LAT,LON))
pal_acc_ncep = np.zeros((WEEK,LAT,LON)) ; pal_pval_ncep = np.zeros((WEEK,LAT,LON))
for i in range(0,LAT):
    for j in range(0,LON):
        for w in range(0,WEEK):
	    pel_acc_ncep[w,i,j] = sp.stats.pearsonr(pel_ncep[:,w,i,j], pel_chirps[:,w,i,j])[0]
	    pel_pval_ncep[w,i,j] = sp.stats.pearsonr(pel_ncep[:,w,i,j], pel_chirps[:,w,i,j])[1]
	    pne_acc_ncep[w,i,j] = sp.stats.pearsonr(pne_ncep[:,w,i,j], pne_chirps[:,w,i,j])[0]
	    pne_pval_ncep[w,i,j] = sp.stats.pearsonr(pne_ncep[:,w,i,j], pne_chirps[:,w,i,j])[1]
	    pla_acc_ncep[w,i,j] = sp.stats.pearsonr(pla_ncep[:,w,i,j], pla_chirps[:,w,i,j])[0]
	    pla_pval_ncep[w,i,j] = sp.stats.pearsonr(pla_ncep[:,w,i,j], pla_chirps[:,w,i,j])[1]
	    pno_acc_ncep[w,i,j] = sp.stats.pearsonr(pno_ncep[:,w,i,j], pno_chirps[:,w,i,j])[0]
	    pno_pval_ncep[w,i,j] = sp.stats.pearsonr(pno_ncep[:,w,i,j], pno_chirps[:,w,i,j])[1]
	    pal_acc_ncep[w,i,j] = sp.stats.pearsonr(pal_ncep[:,w,i,j], pal_chirps[:,w,i,j])[0]
	    pal_pval_ncep[w,i,j] = sp.stats.pearsonr(pal_ncep[:,w,i,j], pal_chirps[:,w,i,j])[1]


ncep_acc = np.stack((pel_acc_ncep, pne_acc_ncep, pla_acc_ncep, pno_acc_ncep, pal_acc_ncep))
ncep_pval = np.stack((pel_pval_ncep, pne_pval_ncep, pla_pval_ncep, pno_pval_ncep, pal_pval_ncep))

pel_acc_ecmf = np.zeros((WEEK,LAT,LON)) ; pel_pval_ecmf = np.zeros((WEEK,LAT,LON))
pne_acc_ecmf = np.zeros((WEEK,LAT,LON)) ; pne_pval_ecmf = np.zeros((WEEK,LAT,LON))
pla_acc_ecmf = np.zeros((WEEK,LAT,LON)) ; pla_pval_ecmf = np.zeros((WEEK,LAT,LON))
pno_acc_ecmf = np.zeros((WEEK,LAT,LON)) ; pno_pval_ecmf = np.zeros((WEEK,LAT,LON))
pal_acc_ecmf = np.zeros((WEEK,LAT,LON)) ; pal_pval_ecmf = np.zeros((WEEK,LAT,LON))
for i in range(0,LAT):
    for j in range(0,LON):
        for w in range(0,WEEK):
	    pel_acc_ecmf[w,i,j] = sp.stats.pearsonr(pel_ecmf[:,w,i,j], pel_chirps[:,w,i,j])[0]
	    pel_pval_ecmf[w,i,j] = sp.stats.pearsonr(pel_ecmf[:,w,i,j], pel_chirps[:,w,i,j])[1]
	    pne_acc_ecmf[w,i,j] = sp.stats.pearsonr(pne_ecmf[:,w,i,j], pne_chirps[:,w,i,j])[0]
	    pne_pval_ecmf[w,i,j] = sp.stats.pearsonr(pne_ecmf[:,w,i,j], pne_chirps[:,w,i,j])[1]
	    pla_acc_ecmf[w,i,j] = sp.stats.pearsonr(pla_ecmf[:,w,i,j], pla_chirps[:,w,i,j])[0]
	    pla_pval_ecmf[w,i,j] = sp.stats.pearsonr(pla_ecmf[:,w,i,j], pla_chirps[:,w,i,j])[1]
	    pno_acc_ecmf[w,i,j] = sp.stats.pearsonr(pno_ecmf[:,w,i,j], pno_chirps[:,w,i,j])[0]
	    pno_pval_ecmf[w,i,j] = sp.stats.pearsonr(pno_ecmf[:,w,i,j], pno_chirps[:,w,i,j])[1]
	    pal_acc_ecmf[w,i,j] = sp.stats.pearsonr(pal_ecmf[:,w,i,j], pal_chirps[:,w,i,j])[0]
	    pal_pval_ecmf[w,i,j] = sp.stats.pearsonr(pal_ecmf[:,w,i,j], pal_chirps[:,w,i,j])[1]


ecmf_acc = np.stack((pel_acc_ecmf, pne_acc_ecmf, pla_acc_ecmf, pno_acc_ecmf, pal_acc_ecmf))
ecmf_pval = np.stack((pel_pval_ecmf, pne_pval_ecmf, pla_pval_ecmf, pno_pval_ecmf, pal_pval_ecmf))

pel_acc_bam = np.zeros((WEEK,LAT,LON)) ; pel_pval_bam = np.zeros((WEEK,LAT,LON))
pne_acc_bam = np.zeros((WEEK,LAT,LON)) ; pne_pval_bam = np.zeros((WEEK,LAT,LON))
pla_acc_bam = np.zeros((WEEK,LAT,LON)) ; pla_pval_bam = np.zeros((WEEK,LAT,LON))
pno_acc_bam = np.zeros((WEEK,LAT,LON)) ; pno_pval_bam = np.zeros((WEEK,LAT,LON))
pal_acc_bam = np.zeros((WEEK,LAT,LON)) ; pal_pval_bam = np.zeros((WEEK,LAT,LON))
for i in range(0,LAT):
    for j in range(0,LON):
        for w in range(0,WEEK):
	    pel_acc_bam[w,i,j] = sp.stats.pearsonr(pel_bam[:,w,i,j], pel_bchirps[:,w,i,j])[0]
	    pel_pval_bam[w,i,j] = sp.stats.pearsonr(pel_bam[:,w,i,j], pel_bchirps[:,w,i,j])[1]
	    pne_acc_bam[w,i,j] = sp.stats.pearsonr(pne_bam[:,w,i,j], pne_bchirps[:,w,i,j])[0]
	    pne_pval_bam[w,i,j] = sp.stats.pearsonr(pne_bam[:,w,i,j], pne_bchirps[:,w,i,j])[1]
	    pla_acc_bam[w,i,j] = sp.stats.pearsonr(pla_bam[:,w,i,j], pla_bchirps[:,w,i,j])[0]
	    pla_pval_bam[w,i,j] = sp.stats.pearsonr(pla_bam[:,w,i,j], pla_bchirps[:,w,i,j])[1]
	    pno_acc_bam[w,i,j] = sp.stats.pearsonr(pno_bam[:,w,i,j], pno_bchirps[:,w,i,j])[0]
	    pno_pval_bam[w,i,j] = sp.stats.pearsonr(pno_bam[:,w,i,j], pno_bchirps[:,w,i,j])[1]
	    pal_acc_bam[w,i,j] = sp.stats.pearsonr(pal_bam[:,w,i,j], pal_bchirps[:,w,i,j])[0]
	    pal_pval_bam[w,i,j] = sp.stats.pearsonr(pal_bam[:,w,i,j], pal_bchirps[:,w,i,j])[1]


bam_acc = np.stack((pel_acc_bam, pne_acc_bam, pla_acc_bam, pno_acc_bam, pal_acc_bam))
bam_pval = np.stack((pel_pval_bam, pne_pval_bam, pla_pval_bam, pno_pval_bam, pal_pval_bam))

acc = np.stack((bam_acc, ecmf_acc, ncep_acc, ukmo_acc))
pval = np.stack((bam_pval, ecmf_pval, ncep_pval, ukmo_pval))
pval = np.ma.masked_less(pval, 0.05)
pval[np.isnan(acc)]=np.nan
'''

bt = np.load('enso_boot.npz')
racc = bt['racc']
rboot = bt['rboot']
#rboot[(rboot>=5)&(rboot<=95)] = np.nan
rboot = (rboot/900)*100

latlim = [-49.5,19.5]
lonlim = [-90,-20]
lon = np.arange(lonlim[0],lonlim[1]+0.5,1.5)
lat = np.arange(latlim[0],latlim[1]+0.5,1.5)
cols = 'RdYlBu'
cmap= matplotlib.cm.RdYlBu
cmin = -0.3
cmax = 0.3
cspc = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel = 'CC-Diff'
dlabel = 'CC difference'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
lw = 1
gl = 20
latext = [-50,15]
lonext = [-90,-30]
data = ['BAM','ECMWF','NCEP','UKMO']
phase = ['ElNino','Neutral','LaNina', 'ENSO']
WEEK=5

rg_names = ['NSA','AMZ','NDE','SESA','AND','PAT']
rg_lon_min = [-80,-67,-47,-60,-75,-75]
rg_lon_max = [-50,-47,-34,-48,-67,-60]
rg_lat_min = [0,-15,-15,-35,-40,-50]
rg_lat_max = [12,-5,-5,-22,-15,-40]
id_lon_min = [7,16,29,20,10,10]
id_lon_max = [27,29,38,28,16,20]
id_lat_min = [33,23,23,10,7,0]
id_lat_max = [41,30,30,19,23,7]

fmt_ls  = ['-s','-P','-X','-o']
ms_ls   = [7,9,9,8]
mw_ls   = [1,1,1,1]
col_ls  = ['black','forestgreen','dodgerblue','firebrick']

week = np.asarray(range(WEEK))+1

#PLOT
for r in range(len(rg_names)):
    fig = plt.figure(figsize=(3.5,3))
    plt.axhline(0, color='0.5', ls='--')
    for d in range(len(data)):
        diff = (racc[d,3,:,r]-racc[d,1,:,r])
	rdiff = rboot[d,3,:,r].copy()
        plt.plot(week,diff,fmt_ls[d],color=col_ls[d],linewidth=2,ms=ms_ls[d],alpha=0.9,mew=mw_ls[d])
        diff[rdiff>70] = np.nan
        plt.plot(week,diff,fmt_ls[d],color=col_ls[d],linewidth=2,ms=ms_ls[d], markeredgecolor=col_ls[d], markerfacecolor='w')
    plt.xticks(week, fontsize=13)
    plt.xlabel('Weeks')
    plt.yticks(clevs, fontsize=13)
    plt.ylabel(dlabel)
#    plt.ylim(ymin=0, ymax=1)
#    plt.legend(loc='lower right', prop={'size': 9})
#    plt.title(rg_names[r], fontsize=15)
    fig.tight_layout()
    plt.savefig(str(rg_names[r])+'_ENSO-Neutral_'+str(clabel)+'.eps')
#    plt.show()
    plt.close()


#PLOT
for r in range(len(rg_names)):
    fig = plt.figure(figsize=(3.5,3))
    plt.axhline(0, color='0.5', ls='--')
    for d in range(len(data)):
        diff = (racc[d,0,:,r]-racc[d,1,:,r])
	rdiff = rboot[d,0,:,r].copy()
        plt.plot(week,diff,fmt_ls[d],color=col_ls[d],linewidth=2,ms=ms_ls[d],alpha=0.9,mew=mw_ls[d])
        diff[rdiff>70] = np.nan
        plt.plot(week,diff,fmt_ls[d],color=col_ls[d],linewidth=2,ms=ms_ls[d], markeredgecolor=col_ls[d], markerfacecolor='w')
    plt.xticks(week, fontsize=13)
    plt.xlabel('Weeks')
    plt.yticks(np.arange(-0.3,0.4,0.1), fontsize=13)
    plt.ylabel(dlabel)
#    plt.ylim(ymin=0, ymax=1)
#    plt.legend(loc='lower right', prop={'size': 9})
#    plt.title(rg_names[r], fontsize=15)
    fig.tight_layout()
    plt.savefig(str(rg_names[r])+'_Elnino-Neutral_'+str(clabel)+'.eps')
#    plt.show()
    plt.close()
    
#PLOT
for r in range(len(rg_names)):
    fig = plt.figure(figsize=(3.5,3))
    plt.axhline(0, color='0.5', ls='--')
    for d in range(len(data)):
        diff = (racc[d,2,:,r]-racc[d,1,:,r])
	rdiff = rboot[d,2,:,r].copy()
        plt.plot(week,diff,fmt_ls[d],color=col_ls[d],linewidth=2,ms=ms_ls[d],alpha=0.9,mew=mw_ls[d])
        diff[rdiff>70] = np.nan
        plt.plot(week,diff,fmt_ls[d],color=col_ls[d],linewidth=2,ms=ms_ls[d], markeredgecolor=col_ls[d], markerfacecolor='w')
    plt.xticks(week, fontsize=13)
    plt.xlabel('Weeks')
    plt.yticks(np.arange(-0.3,0.4,0.1), fontsize=13)
    plt.ylabel(dlabel)
#    plt.ylim(ymin=0, ymax=1)
#    plt.legend(loc='lower right', prop={'size': 9})
#    plt.title(rg_names[r], fontsize=15)
    fig.tight_layout()
    plt.savefig(str(rg_names[r])+'_LaNina-Neutral_'+str(clabel)+'.eps')
#    plt.show()


#PLOT LEGEND
fig = plt.figure(figsize=(3.5,3))
for d in range(len(data)):
    diff = (racc[d,3,:,0]-racc[d,1,:,0])
    plt.plot(week,diff,fmt_ls[d],color=col_ls[d],linewidth=2,ms=ms_ls[d],alpha=0.9,mew=mw_ls[d], label=data[d])
    plt.xticks(week, fontsize=13)
plt.ylim(-5,1)
plt.legend(loc='lower right', prop={'size': 15})
fig.tight_layout()
plt.savefig(str(clabel)+'_legend.eps')
#plt.show()


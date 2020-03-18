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
N0 = 26 ; NB0 = 14
N1 = 37 ; NB1 = 15
N2 = 38 ; NB2 = 22
N3 = 47 ; NB3 = 21
N4 = 72 ; NB4 = 38
N5 = 220; NB5 = 110
REG = 6

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
phase = ['MJO-8+1','MJO-2+3','MJO-4+5','MJO-6+7']
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
ms_ls   = [6,10,10,8]
mw_ls   = [1,2,2,1]
col_ls  = ['black','forestgreen','dodgerblue','firebrick']

week = np.asarray(range(WEEK))+1

bt = np.load('mjo_boot_ph0.npz')
racc = bt['racc']
rboot = bt['rboot']
#rboot[(rboot>=5)&(rboot<=95)] = np.nan
rboot = (rboot/900)*100

#PLOT
for r in range(len(rg_names)):
    print rg_names[r]
    for p in range(len(phase)):
        fig = plt.figure(figsize=(3.5,3))
        plt.axhline(0, color='0.5', ls='--')
        for d in range(len(data)):
            diff = (racc[d,p,:,r]-racc[d,4,:,r])
	    rdiff = rboot[d,p,:,r].copy()
            plt.plot(week,diff,fmt_ls[d],color=col_ls[d],linewidth=2,ms=ms_ls[d],alpha=0.9,mew=mw_ls[d])
            diff[rdiff>70] = np.nan
            plt.plot(week,diff,fmt_ls[d],color=col_ls[d],linewidth=2,ms=ms_ls[d], markeredgecolor=col_ls[d], markerfacecolor='w')
        plt.xticks(week, fontsize=13)
        plt.xlabel('Weeks')
        plt.yticks(clevs, fontsize=13)
        plt.ylabel(dlabel)
#        plt.ylim(ymin=0, ymax=1)
#        plt.legend(loc='lower right', prop={'size': 9})
#        plt.title(rg_names[r], fontsize=15)
        fig.tight_layout()
        plt.savefig(str(rg_names[r])+'_'+str(phase[p])+'_minus_MJO-0_'+str(clabel)+'.eps')
#        plt.show()
        plt.close()

'''
bt = np.load('mjo_boot_all.npz')
racc = bt['racc']
rboot = bt['rboot']
#rboot[(rboot>=5)&(rboot<=95)] = np.nan
rboot = (rboot/900)*100

#PLOT
for r in range(len(rg_names)):
    for p in range(len(phase)):
        fig = plt.figure(figsize=(3.5,3))
        plt.axhline(0, color='0.5', ls='--')
        for d in range(len(data)):
            diff = (racc[d,p,:,r]-racc[d,5,:,r])
	    rdiff = rboot[d,p,:,r].copy()
            plt.plot(week,diff,fmt_ls[d],color=col_ls[d],linewidth=2,ms=ms_ls[d],alpha=0.9,mew=mw_ls[d])
            diff[rdiff>70] = np.nan
            plt.plot(week,diff,fmt_ls[d],color=col_ls[d],linewidth=2,ms=ms_ls[d], markeredgecolor=col_ls[d], markerfacecolor='w')
        plt.xticks(week, fontsize=13)
        plt.xlabel('Weeks')
        plt.yticks(clevs, fontsize=13)
        plt.ylabel(dlabel)
#        plt.ylim(ymin=0, ymax=1)
#        plt.legend(loc='lower right', prop={'size': 9})
#        plt.title(rg_names[r], fontsize=15)
        fig.tight_layout()
        plt.savefig(str(rg_names[r])+'_'+str(phase[p])+'_minus_ALL-YRS_'+str(clabel)+'.eps')
#        plt.show()
        plt.close()
'''

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
plt.close()


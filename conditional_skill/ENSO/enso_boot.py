from __future__ import division
#import cf, cfplot as cfp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
import scipy as sp
from scipy import stats
#import read_enso_index as nino

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


bchirps = np.load('enso_bchirps.npz')
pel_bchirps = bchirps['apel_bchirps']
pne_bchirps = bchirps['apne_bchirps']
pla_bchirps = bchirps['apla_bchirps']
pno_bchirps = bchirps['apno_bchirps']
pal_bchirps = bchirps['apal_bchirps']

chirps = np.load('enso_chirps.npz')
pel_chirps = chirps['apel_chirps']
pne_chirps = chirps['apne_chirps']
pla_chirps = chirps['apla_chirps']
pno_chirps = chirps['apno_chirps']
pal_chirps = chirps['apal_chirps']

bam = np.load('enso_bam.npz')
pel_bam = np.nanmean(bam['apel_bam'], axis=2)
pne_bam = np.nanmean(bam['apne_bam'], axis=2)
pla_bam = np.nanmean(bam['apla_bam'], axis=2)
pno_bam = np.nanmean(bam['apno_bam'], axis=2)
pal_bam = np.nanmean(bam['apal_bam'], axis=2)

ukmo = np.load('../enso_ukmo.npz')
pel_ukmo = np.nanmean(ukmo['apel_ukmo'], axis=2)
pne_ukmo = np.nanmean(ukmo['apne_ukmo'], axis=2)
pla_ukmo = np.nanmean(ukmo['apla_ukmo'], axis=2)
pno_ukmo = np.nanmean(ukmo['apno_ukmo'], axis=2)
pal_ukmo = np.nanmean(ukmo['apal_ukmo'], axis=2)

ecmf = np.load('../enso_ecmf.npz')
pel_ecmf = np.nanmean(ecmf['apel_ecmf'], axis=2)
pne_ecmf = np.nanmean(ecmf['apne_ecmf'], axis=2)
pla_ecmf = np.nanmean(ecmf['apla_ecmf'], axis=2)
pno_ecmf = np.nanmean(ecmf['apno_ecmf'], axis=2)
pal_ecmf = np.nanmean(ecmf['apal_ecmf'], axis=2)

ncep = np.load('../enso_ncep.npz')
pel_ncep = np.nanmean(ncep['apel_ncep'], axis=2)
pne_ncep = np.nanmean(ncep['apne_ncep'], axis=2)
pla_ncep = np.nanmean(ncep['apla_ncep'], axis=2)
pno_ncep = np.nanmean(ncep['apno_ncep'], axis=2)
pal_ncep = np.nanmean(ncep['apal_ncep'], axis=2) 


iel = np.arange(NEL); ine = np.arange(NNE); ila = np.arange(NLA); ino = np.arange(NNO)
biel = np.arange(int(NEL/2)); bine = np.arange(int(NNE/2)); bila = np.arange(int(NLA/2)); bino = np.arange(int(NNO/2))

pel_dist_ukmo = np.zeros((NUM,WEEK,LAT,LON))
pne_dist_ukmo = np.zeros((NUM,WEEK,LAT,LON))
pla_dist_ukmo = np.zeros((NUM,WEEK,LAT,LON))
pno_dist_ukmo = np.zeros((NUM,WEEK,LAT,LON))
pel_dist_ncep = np.zeros((NUM,WEEK,LAT,LON))
pne_dist_ncep = np.zeros((NUM,WEEK,LAT,LON))
pla_dist_ncep = np.zeros((NUM,WEEK,LAT,LON))
pno_dist_ncep = np.zeros((NUM,WEEK,LAT,LON))
pel_dist_ecmf = np.zeros((NUM,WEEK,LAT,LON))
pne_dist_ecmf = np.zeros((NUM,WEEK,LAT,LON))
pla_dist_ecmf = np.zeros((NUM,WEEK,LAT,LON))
pno_dist_ecmf = np.zeros((NUM,WEEK,LAT,LON))
pel_dist_bam = np.zeros((NUM,WEEK,LAT,LON))
pne_dist_bam = np.zeros((NUM,WEEK,LAT,LON))
pla_dist_bam = np.zeros((NUM,WEEK,LAT,LON))
pno_dist_bam = np.zeros((NUM,WEEK,LAT,LON))
for n in range(0,NUM):
    A = np.random.choice(iel, int(NEL), replace=True)
    B = np.random.choice(ine, int(NNE), replace=True)
    C = np.random.choice(ila, int(NLA), replace=True)
    D = np.random.choice(ino, int(NNO), replace=True)
    L = np.random.choice(biel, int(NEL/2), replace=True)
    M = np.random.choice(bine, int(NNE/2), replace=True)
    N = np.random.choice(bila, int(NLA/2), replace=True)
    O = np.random.choice(bino, int(NNO/2), replace=True)
    for i in range(0,LAT):
        for j in range(0,LON):
            for w in range(0,WEEK):
	        pel_dist_ukmo[n,w,i,j] = sp.stats.pearsonr(pel_ukmo[A,w,i,j], pel_chirps[A,w,i,j])[0]
	        pne_dist_ukmo[n,w,i,j] = sp.stats.pearsonr(pne_ukmo[B,w,i,j], pne_chirps[B,w,i,j])[0]
	        pla_dist_ukmo[n,w,i,j] = sp.stats.pearsonr(pla_ukmo[C,w,i,j], pla_chirps[C,w,i,j])[0]
	        pno_dist_ukmo[n,w,i,j] = sp.stats.pearsonr(pno_ukmo[D,w,i,j], pno_chirps[D,w,i,j])[0]
	        pel_dist_ncep[n,w,i,j] = sp.stats.pearsonr(pel_ncep[A,w,i,j], pel_chirps[A,w,i,j])[0]
	        pne_dist_ncep[n,w,i,j] = sp.stats.pearsonr(pne_ncep[B,w,i,j], pne_chirps[B,w,i,j])[0]
	        pla_dist_ncep[n,w,i,j] = sp.stats.pearsonr(pla_ncep[C,w,i,j], pla_chirps[C,w,i,j])[0]
	        pno_dist_ncep[n,w,i,j] = sp.stats.pearsonr(pno_ncep[D,w,i,j], pno_chirps[D,w,i,j])[0]
	        pel_dist_ecmf[n,w,i,j] = sp.stats.pearsonr(pel_ecmf[A,w,i,j], pel_chirps[A,w,i,j])[0]
	        pne_dist_ecmf[n,w,i,j] = sp.stats.pearsonr(pne_ecmf[B,w,i,j], pne_chirps[B,w,i,j])[0]
	        pla_dist_ecmf[n,w,i,j] = sp.stats.pearsonr(pla_ecmf[C,w,i,j], pla_chirps[C,w,i,j])[0]
	        pno_dist_ecmf[n,w,i,j] = sp.stats.pearsonr(pno_ecmf[D,w,i,j], pno_chirps[D,w,i,j])[0]
	        pel_dist_bam[n,w,i,j] = sp.stats.pearsonr(pel_bam[L,w,i,j], pel_bchirps[L,w,i,j])[0]
	        pne_dist_bam[n,w,i,j] = sp.stats.pearsonr(pne_bam[M,w,i,j], pne_bchirps[M,w,i,j])[0]
	        pla_dist_bam[n,w,i,j] = sp.stats.pearsonr(pla_bam[N,w,i,j], pla_bchirps[N,w,i,j])[0]
	        pno_dist_bam[n,w,i,j] = sp.stats.pearsonr(pno_bam[O,w,i,j], pno_bchirps[O,w,i,j])[0]
    print n

pel_acc_ukmo = np.zeros((WEEK,LAT,LON)); pel_boot_ukmo = np.zeros((WEEK,LAT,LON))
pne_acc_ukmo = np.zeros((WEEK,LAT,LON)); pne_boot_ukmo = np.zeros((WEEK,LAT,LON))
pla_acc_ukmo = np.zeros((WEEK,LAT,LON)); pla_boot_ukmo = np.zeros((WEEK,LAT,LON))
pno_acc_ukmo = np.zeros((WEEK,LAT,LON)); pno_boot_ukmo = np.zeros((WEEK,LAT,LON))
pal_acc_ukmo = np.zeros((WEEK,LAT,LON))
pel_acc_bam = np.zeros((WEEK,LAT,LON)); pel_boot_bam = np.zeros((WEEK,LAT,LON))
pne_acc_bam = np.zeros((WEEK,LAT,LON)); pne_boot_bam = np.zeros((WEEK,LAT,LON))
pla_acc_bam = np.zeros((WEEK,LAT,LON)); pla_boot_bam = np.zeros((WEEK,LAT,LON))
pno_acc_bam = np.zeros((WEEK,LAT,LON)); pno_boot_bam = np.zeros((WEEK,LAT,LON))
pal_acc_bam = np.zeros((WEEK,LAT,LON))
pel_acc_ecmf = np.zeros((WEEK,LAT,LON)); pel_boot_ecmf = np.zeros((WEEK,LAT,LON))
pne_acc_ecmf = np.zeros((WEEK,LAT,LON)); pne_boot_ecmf = np.zeros((WEEK,LAT,LON))
pla_acc_ecmf = np.zeros((WEEK,LAT,LON)); pla_boot_ecmf = np.zeros((WEEK,LAT,LON))
pno_acc_ecmf = np.zeros((WEEK,LAT,LON)); pno_boot_ecmf = np.zeros((WEEK,LAT,LON))
pal_acc_ecmf = np.zeros((WEEK,LAT,LON))
pel_acc_ncep = np.zeros((WEEK,LAT,LON)); pel_boot_ncep = np.zeros((WEEK,LAT,LON))
pne_acc_ncep = np.zeros((WEEK,LAT,LON)); pne_boot_ncep = np.zeros((WEEK,LAT,LON))
pla_acc_ncep = np.zeros((WEEK,LAT,LON)); pla_boot_ncep = np.zeros((WEEK,LAT,LON))
pno_acc_ncep = np.zeros((WEEK,LAT,LON)); pno_boot_ncep = np.zeros((WEEK,LAT,LON))
pal_acc_ncep = np.zeros((WEEK,LAT,LON))
for i in range(0,LAT):
    for j in range(0,LON):
        for w in range(0,WEEK):
	    pel_acc_bam[w,i,j] = sp.stats.pearsonr(pel_bam[:,w,i,j], pel_bchirps[:,w,i,j])[0]
	    pne_acc_bam[w,i,j] = sp.stats.pearsonr(pne_bam[:,w,i,j], pne_bchirps[:,w,i,j])[0]
	    pla_acc_bam[w,i,j] = sp.stats.pearsonr(pla_bam[:,w,i,j], pla_bchirps[:,w,i,j])[0]
	    pno_acc_bam[w,i,j] = sp.stats.pearsonr(pno_bam[:,w,i,j], pno_bchirps[:,w,i,j])[0]
	    pal_acc_bam[w,i,j] = sp.stats.pearsonr(pal_bam[:,w,i,j], pal_bchirps[:,w,i,j])[0]
	    pel_acc_ecmf[w,i,j] = sp.stats.pearsonr(pel_ecmf[:,w,i,j], pel_chirps[:,w,i,j])[0]
	    pne_acc_ecmf[w,i,j] = sp.stats.pearsonr(pne_ecmf[:,w,i,j], pne_chirps[:,w,i,j])[0]
	    pla_acc_ecmf[w,i,j] = sp.stats.pearsonr(pla_ecmf[:,w,i,j], pla_chirps[:,w,i,j])[0]
	    pno_acc_ecmf[w,i,j] = sp.stats.pearsonr(pno_ecmf[:,w,i,j], pno_chirps[:,w,i,j])[0]
	    pal_acc_ecmf[w,i,j] = sp.stats.pearsonr(pal_ecmf[:,w,i,j], pal_chirps[:,w,i,j])[0]
	    pel_acc_ncep[w,i,j] = sp.stats.pearsonr(pel_ncep[:,w,i,j], pel_chirps[:,w,i,j])[0]
	    pne_acc_ncep[w,i,j] = sp.stats.pearsonr(pne_ncep[:,w,i,j], pne_chirps[:,w,i,j])[0]
	    pla_acc_ncep[w,i,j] = sp.stats.pearsonr(pla_ncep[:,w,i,j], pla_chirps[:,w,i,j])[0]
	    pno_acc_ncep[w,i,j] = sp.stats.pearsonr(pno_ncep[:,w,i,j], pno_chirps[:,w,i,j])[0]
	    pal_acc_ncep[w,i,j] = sp.stats.pearsonr(pal_ncep[:,w,i,j], pal_chirps[:,w,i,j])[0]
	    pel_acc_ukmo[w,i,j] = sp.stats.pearsonr(pel_ukmo[:,w,i,j], pel_chirps[:,w,i,j])[0]
	    pne_acc_ukmo[w,i,j] = sp.stats.pearsonr(pne_ukmo[:,w,i,j], pne_chirps[:,w,i,j])[0]
	    pla_acc_ukmo[w,i,j] = sp.stats.pearsonr(pla_ukmo[:,w,i,j], pla_chirps[:,w,i,j])[0]
	    pno_acc_ukmo[w,i,j] = sp.stats.pearsonr(pno_ukmo[:,w,i,j], pno_chirps[:,w,i,j])[0]
	    pal_acc_ukmo[w,i,j] = sp.stats.pearsonr(pal_ukmo[:,w,i,j], pal_chirps[:,w,i,j])[0]
            small = pel_dist_ukmo[:,w,i,j][(pel_dist_ukmo[:,w,i,j]<=np.nanpercentile(pel_dist_ukmo[:,w,i,j], 95))&(pel_dist_ukmo[:,w,i,j]>=np.nanpercentile(pel_dist_ukmo[:,w,i,j], 5))]; pel_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ukmo[:,w,i,j], 5))) 
            small = pne_dist_ukmo[:,w,i,j][(pne_dist_ukmo[:,w,i,j]<=np.nanpercentile(pne_dist_ukmo[:,w,i,j], 95))&(pne_dist_ukmo[:,w,i,j]>=np.nanpercentile(pne_dist_ukmo[:,w,i,j], 5))]; pne_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ukmo[:,w,i,j], 5))) 
            small = pla_dist_ukmo[:,w,i,j][(pla_dist_ukmo[:,w,i,j]<=np.nanpercentile(pla_dist_ukmo[:,w,i,j], 95))&(pla_dist_ukmo[:,w,i,j]>=np.nanpercentile(pla_dist_ukmo[:,w,i,j], 5))]; pla_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ukmo[:,w,i,j], 5))) 
            small = pno_dist_ukmo[:,w,i,j][(pno_dist_ukmo[:,w,i,j]<=np.nanpercentile(pno_dist_ukmo[:,w,i,j], 95))&(pno_dist_ukmo[:,w,i,j]>=np.nanpercentile(pno_dist_ukmo[:,w,i,j], 5))]; pno_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ukmo[:,w,i,j], 5))) 
            small = pel_dist_ncep[:,w,i,j][(pel_dist_ncep[:,w,i,j]<=np.nanpercentile(pel_dist_ncep[:,w,i,j], 95))&(pel_dist_ncep[:,w,i,j]>=np.nanpercentile(pel_dist_ncep[:,w,i,j], 5))]; pel_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ncep[:,w,i,j], 5))) 
            small = pne_dist_ncep[:,w,i,j][(pne_dist_ncep[:,w,i,j]<=np.nanpercentile(pne_dist_ncep[:,w,i,j], 95))&(pne_dist_ncep[:,w,i,j]>=np.nanpercentile(pne_dist_ncep[:,w,i,j], 5))]; pne_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ncep[:,w,i,j], 5))) 
            small = pla_dist_ncep[:,w,i,j][(pla_dist_ncep[:,w,i,j]<=np.nanpercentile(pla_dist_ncep[:,w,i,j], 95))&(pla_dist_ncep[:,w,i,j]>=np.nanpercentile(pla_dist_ncep[:,w,i,j], 5))]; pla_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ncep[:,w,i,j], 5))) 
            small = pno_dist_ncep[:,w,i,j][(pno_dist_ncep[:,w,i,j]<=np.nanpercentile(pno_dist_ncep[:,w,i,j], 95))&(pno_dist_ncep[:,w,i,j]>=np.nanpercentile(pno_dist_ncep[:,w,i,j], 5))]; pno_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ncep[:,w,i,j], 5))) 
            small = pel_dist_ecmf[:,w,i,j][(pel_dist_ecmf[:,w,i,j]<=np.nanpercentile(pel_dist_ecmf[:,w,i,j], 95))&(pel_dist_ecmf[:,w,i,j]>=np.nanpercentile(pel_dist_ecmf[:,w,i,j], 5))]; pel_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ecmf[:,w,i,j], 5))) 
            small = pne_dist_ecmf[:,w,i,j][(pne_dist_ecmf[:,w,i,j]<=np.nanpercentile(pne_dist_ecmf[:,w,i,j], 95))&(pne_dist_ecmf[:,w,i,j]>=np.nanpercentile(pne_dist_ecmf[:,w,i,j], 5))]; pne_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ecmf[:,w,i,j], 5))) 
            small = pla_dist_ecmf[:,w,i,j][(pla_dist_ecmf[:,w,i,j]<=np.nanpercentile(pla_dist_ecmf[:,w,i,j], 95))&(pla_dist_ecmf[:,w,i,j]>=np.nanpercentile(pla_dist_ecmf[:,w,i,j], 5))]; pla_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ecmf[:,w,i,j], 5))) 
            small = pno_dist_ecmf[:,w,i,j][(pno_dist_ecmf[:,w,i,j]<=np.nanpercentile(pno_dist_ecmf[:,w,i,j], 95))&(pno_dist_ecmf[:,w,i,j]>=np.nanpercentile(pno_dist_ecmf[:,w,i,j], 5))]; pno_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_ecmf[:,w,i,j], 5))) 
            small = pel_dist_bam[:,w,i,j][(pel_dist_bam[:,w,i,j]<=np.nanpercentile(pel_dist_bam[:,w,i,j], 95))&(pel_dist_bam[:,w,i,j]>=np.nanpercentile(pel_dist_bam[:,w,i,j], 5))]; pel_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_bam[:,w,i,j], 5))) 
            small = pne_dist_bam[:,w,i,j][(pne_dist_bam[:,w,i,j]<=np.nanpercentile(pne_dist_bam[:,w,i,j], 95))&(pne_dist_bam[:,w,i,j]>=np.nanpercentile(pne_dist_bam[:,w,i,j], 5))]; pne_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_bam[:,w,i,j], 5))) 
            small = pla_dist_bam[:,w,i,j][(pla_dist_bam[:,w,i,j]<=np.nanpercentile(pla_dist_bam[:,w,i,j], 95))&(pla_dist_bam[:,w,i,j]>=np.nanpercentile(pla_dist_bam[:,w,i,j], 5))]; pla_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_bam[:,w,i,j], 5))) 
            small = pno_dist_bam[:,w,i,j][(pno_dist_bam[:,w,i,j]<=np.nanpercentile(pno_dist_bam[:,w,i,j], 95))&(pno_dist_bam[:,w,i,j]>=np.nanpercentile(pno_dist_bam[:,w,i,j], 5))]; pno_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(pne_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(pne_dist_bam[:,w,i,j], 5))) 


rg_names = ['NSA','AMZ','NDE','SESA','AND','PAT']
rg_lon_min = [-80,-67,-47,-60,-75,-75]
rg_lon_max = [-50,-47,-34,-48,-67,-60]
rg_lat_min = [0,-15,-15,-35,-40,-50]
rg_lat_max = [12,-5,-5,-22,-15,-40]
id_lon_min = [7,16,29,20,10,10]
id_lon_max = [27,29,38,28,16,20]
id_lat_min = [33,23,23,10,7,0]
id_lat_max = [41,30,30,19,23,7]

pel_rdist_ukmo = np.zeros((NUM,WEEK,REG))
pne_rdist_ukmo = np.zeros((NUM,WEEK,REG))
pla_rdist_ukmo = np.zeros((NUM,WEEK,REG))
pno_rdist_ukmo = np.zeros((NUM,WEEK,REG))
pel_rdist_ncep = np.zeros((NUM,WEEK,REG))
pne_rdist_ncep = np.zeros((NUM,WEEK,REG))
pla_rdist_ncep = np.zeros((NUM,WEEK,REG))
pno_rdist_ncep = np.zeros((NUM,WEEK,REG))
pel_rdist_ecmf = np.zeros((NUM,WEEK,REG))
pne_rdist_ecmf = np.zeros((NUM,WEEK,REG))
pla_rdist_ecmf = np.zeros((NUM,WEEK,REG))
pno_rdist_ecmf = np.zeros((NUM,WEEK,REG))
pel_rdist_bam = np.zeros((NUM,WEEK,REG))
pne_rdist_bam = np.zeros((NUM,WEEK,REG))
pla_rdist_bam = np.zeros((NUM,WEEK,REG))
pno_rdist_bam = np.zeros((NUM,WEEK,REG))
pel_racc_ukmo = np.zeros((WEEK,REG)); pel_rboot_ukmo = np.zeros((WEEK,REG))
pne_racc_ukmo = np.zeros((WEEK,REG)); pne_rboot_ukmo = np.zeros((WEEK,REG))
pla_racc_ukmo = np.zeros((WEEK,REG)); pla_rboot_ukmo = np.zeros((WEEK,REG))
pno_racc_ukmo = np.zeros((WEEK,REG)); pno_rboot_ukmo = np.zeros((WEEK,REG))
pal_racc_ukmo = np.zeros((WEEK,REG))
pel_racc_bam = np.zeros((WEEK,REG)); pel_rboot_bam = np.zeros((WEEK,REG))
pne_racc_bam = np.zeros((WEEK,REG)); pne_rboot_bam = np.zeros((WEEK,REG))
pla_racc_bam = np.zeros((WEEK,REG)); pla_rboot_bam = np.zeros((WEEK,REG))
pno_racc_bam = np.zeros((WEEK,REG)); pno_rboot_bam = np.zeros((WEEK,REG))
pal_racc_bam = np.zeros((WEEK,REG))
pel_racc_ecmf = np.zeros((WEEK,REG)); pel_rboot_ecmf = np.zeros((WEEK,REG))
pne_racc_ecmf = np.zeros((WEEK,REG)); pne_rboot_ecmf = np.zeros((WEEK,REG))
pla_racc_ecmf = np.zeros((WEEK,REG)); pla_rboot_ecmf = np.zeros((WEEK,REG))
pno_racc_ecmf = np.zeros((WEEK,REG)); pno_rboot_ecmf = np.zeros((WEEK,REG))
pal_racc_ecmf = np.zeros((WEEK,REG))
pel_racc_ncep = np.zeros((WEEK,REG)); pel_rboot_ncep = np.zeros((WEEK,REG))
pne_racc_ncep = np.zeros((WEEK,REG)); pne_rboot_ncep = np.zeros((WEEK,REG))
pla_racc_ncep = np.zeros((WEEK,REG)); pla_rboot_ncep = np.zeros((WEEK,REG))
pno_racc_ncep = np.zeros((WEEK,REG)); pno_rboot_ncep = np.zeros((WEEK,REG))
pal_racc_ncep = np.zeros((WEEK,REG))

for r in range(len(rg_names)):
    pel_rdist_ukmo[:,:,r] = np.nanmean(pel_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pne_rdist_ukmo[:,:,r] = np.nanmean(pne_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pla_rdist_ukmo[:,:,r] = np.nanmean(pla_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pno_rdist_ukmo[:,:,r] = np.nanmean(pno_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pel_rdist_ncep[:,:,r] = np.nanmean(pel_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pne_rdist_ncep[:,:,r] = np.nanmean(pne_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pla_rdist_ncep[:,:,r] = np.nanmean(pla_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pno_rdist_ncep[:,:,r] = np.nanmean(pno_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pel_rdist_ecmf[:,:,r] = np.nanmean(pel_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pne_rdist_ecmf[:,:,r] = np.nanmean(pne_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pla_rdist_ecmf[:,:,r] = np.nanmean(pla_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pno_rdist_ecmf[:,:,r] = np.nanmean(pno_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pel_rdist_bam[:,:,r] = np.nanmean(pel_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pne_rdist_bam[:,:,r] = np.nanmean(pne_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pla_rdist_bam[:,:,r] = np.nanmean(pla_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    pno_rdist_bam[:,:,r] = np.nanmean(pno_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    for w in range(0,WEEK):
        pel_racc_ukmo[w,r] = np.nanmean(pel_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pne_racc_ukmo[w,r] = np.nanmean(pne_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pla_racc_ukmo[w,r] = np.nanmean(pla_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pno_racc_ukmo[w,r] = np.nanmean(pno_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pal_racc_ukmo[w,r] = np.nanmean(pal_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pel_racc_ncep[w,r] = np.nanmean(pel_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pne_racc_ncep[w,r] = np.nanmean(pne_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pla_racc_ncep[w,r] = np.nanmean(pla_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pno_racc_ncep[w,r] = np.nanmean(pno_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pal_racc_ncep[w,r] = np.nanmean(pal_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pel_racc_ecmf[w,r] = np.nanmean(pel_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pne_racc_ecmf[w,r] = np.nanmean(pne_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pla_racc_ecmf[w,r] = np.nanmean(pla_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pno_racc_ecmf[w,r] = np.nanmean(pno_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pal_racc_ecmf[w,r] = np.nanmean(pal_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pel_racc_bam[w,r] = np.nanmean(pel_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pne_racc_bam[w,r] = np.nanmean(pne_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pla_racc_bam[w,r] = np.nanmean(pla_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pno_racc_bam[w,r] = np.nanmean(pno_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        pal_racc_bam[w,r] = np.nanmean(pal_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        small = pel_rdist_ukmo[:,w,r][(pel_rdist_ukmo[:,w,r]<=np.nanpercentile(pel_rdist_ukmo[:,w,r], 95))&(pel_rdist_ukmo[:,w,r]>=np.nanpercentile(pel_rdist_ukmo[:,w,r], 5))]; pel_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ukmo[:,w,r], 5))) 
        small = pne_rdist_ukmo[:,w,r][(pne_rdist_ukmo[:,w,r]<=np.nanpercentile(pne_rdist_ukmo[:,w,r], 95))&(pne_rdist_ukmo[:,w,r]>=np.nanpercentile(pne_rdist_ukmo[:,w,r], 5))]; pne_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ukmo[:,w,r], 5))) 
        small = pla_rdist_ukmo[:,w,r][(pla_rdist_ukmo[:,w,r]<=np.nanpercentile(pla_rdist_ukmo[:,w,r], 95))&(pla_rdist_ukmo[:,w,r]>=np.nanpercentile(pla_rdist_ukmo[:,w,r], 5))]; pla_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ukmo[:,w,r], 5))) 
        small = pno_rdist_ukmo[:,w,r][(pno_rdist_ukmo[:,w,r]<=np.nanpercentile(pno_rdist_ukmo[:,w,r], 95))&(pno_rdist_ukmo[:,w,r]>=np.nanpercentile(pno_rdist_ukmo[:,w,r], 5))]; pno_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ukmo[:,w,r], 5))) 
        small = pel_rdist_ncep[:,w,r][(pel_rdist_ncep[:,w,r]<=np.nanpercentile(pel_rdist_ncep[:,w,r], 95))&(pel_rdist_ncep[:,w,r]>=np.nanpercentile(pel_rdist_ncep[:,w,r], 5))]; pel_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ncep[:,w,r], 5))) 
        small = pne_rdist_ncep[:,w,r][(pne_rdist_ncep[:,w,r]<=np.nanpercentile(pne_rdist_ncep[:,w,r], 95))&(pne_rdist_ncep[:,w,r]>=np.nanpercentile(pne_rdist_ncep[:,w,r], 5))]; pne_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ncep[:,w,r], 5))) 
        small = pla_rdist_ncep[:,w,r][(pla_rdist_ncep[:,w,r]<=np.nanpercentile(pla_rdist_ncep[:,w,r], 95))&(pla_rdist_ncep[:,w,r]>=np.nanpercentile(pla_rdist_ncep[:,w,r], 5))]; pla_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ncep[:,w,r], 5))) 
        small = pno_rdist_ncep[:,w,r][(pno_rdist_ncep[:,w,r]<=np.nanpercentile(pno_rdist_ncep[:,w,r], 95))&(pno_rdist_ncep[:,w,r]>=np.nanpercentile(pno_rdist_ncep[:,w,r], 5))]; pno_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ncep[:,w,r], 5))) 
        small = pel_rdist_ecmf[:,w,r][(pel_rdist_ecmf[:,w,r]<=np.nanpercentile(pel_rdist_ecmf[:,w,r], 95))&(pel_rdist_ecmf[:,w,r]>=np.nanpercentile(pel_rdist_ecmf[:,w,r], 5))]; pel_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ecmf[:,w,r], 5))) 
        small = pne_rdist_ecmf[:,w,r][(pne_rdist_ecmf[:,w,r]<=np.nanpercentile(pne_rdist_ecmf[:,w,r], 95))&(pne_rdist_ecmf[:,w,r]>=np.nanpercentile(pne_rdist_ecmf[:,w,r], 5))]; pne_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ecmf[:,w,r], 5))) 
        small = pla_rdist_ecmf[:,w,r][(pla_rdist_ecmf[:,w,r]<=np.nanpercentile(pla_rdist_ecmf[:,w,r], 95))&(pla_rdist_ecmf[:,w,r]>=np.nanpercentile(pla_rdist_ecmf[:,w,r], 5))]; pla_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ecmf[:,w,r], 5))) 
        small = pno_rdist_ecmf[:,w,r][(pno_rdist_ecmf[:,w,r]<=np.nanpercentile(pno_rdist_ecmf[:,w,r], 95))&(pno_rdist_ecmf[:,w,r]>=np.nanpercentile(pno_rdist_ecmf[:,w,r], 5))]; pno_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_ecmf[:,w,r], 5))) 
        small = pel_rdist_bam[:,w,r][(pel_rdist_bam[:,w,r]<=np.nanpercentile(pel_rdist_bam[:,w,r], 95))&(pel_rdist_bam[:,w,r]>=np.nanpercentile(pel_rdist_bam[:,w,r], 5))]; pel_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_bam[:,w,r], 5))) 
        small = pne_rdist_bam[:,w,r][(pne_rdist_bam[:,w,r]<=np.nanpercentile(pne_rdist_bam[:,w,r], 95))&(pne_rdist_bam[:,w,r]>=np.nanpercentile(pne_rdist_bam[:,w,r], 5))]; pne_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_bam[:,w,r], 5))) 
        small = pla_rdist_bam[:,w,r][(pla_rdist_bam[:,w,r]<=np.nanpercentile(pla_rdist_bam[:,w,r], 95))&(pla_rdist_bam[:,w,r]>=np.nanpercentile(pla_rdist_bam[:,w,r], 5))]; pla_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_bam[:,w,r], 5))) 
        small = pno_rdist_bam[:,w,r][(pno_rdist_bam[:,w,r]<=np.nanpercentile(pno_rdist_bam[:,w,r], 95))&(pno_rdist_bam[:,w,r]>=np.nanpercentile(pno_rdist_bam[:,w,r], 5))]; pno_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(pne_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(pne_rdist_bam[:,w,r], 5))) 
    print rg_names[r]

bam_acc = np.stack((pel_acc_bam, pne_acc_bam, pla_acc_bam, pno_acc_bam, pal_acc_bam))
bam_racc = np.stack((pel_racc_bam, pne_racc_bam, pla_racc_bam, pno_racc_bam, pal_racc_bam))
bam_boot = np.stack((pel_boot_bam, pne_boot_bam, pla_boot_bam, pno_boot_bam))
bam_rboot = np.stack((pel_rboot_bam, pne_rboot_bam, pla_rboot_bam, pno_rboot_bam))
ecmf_acc = np.stack((pel_acc_ecmf, pne_acc_ecmf, pla_acc_ecmf, pno_acc_ecmf, pal_acc_ecmf))
ecmf_racc = np.stack((pel_racc_ecmf, pne_racc_ecmf, pla_racc_ecmf, pno_racc_ecmf, pal_racc_ecmf))
ecmf_boot = np.stack((pel_boot_ecmf, pne_boot_ecmf, pla_boot_ecmf, pno_boot_ecmf))
ecmf_rboot = np.stack((pel_rboot_ecmf, pne_rboot_ecmf, pla_rboot_ecmf, pno_rboot_ecmf))
ncep_acc = np.stack((pel_acc_ncep, pne_acc_ncep, pla_acc_ncep, pno_acc_ncep, pal_acc_ncep))
ncep_racc = np.stack((pel_racc_ncep, pne_racc_ncep, pla_racc_ncep, pno_racc_ncep, pal_racc_ncep))
ncep_boot = np.stack((pel_boot_ncep, pne_boot_ncep, pla_boot_ncep, pno_boot_ncep))
ncep_rboot = np.stack((pel_rboot_ncep, pne_rboot_ncep, pla_rboot_ncep, pno_rboot_ncep))
ukmo_acc = np.stack((pel_acc_ukmo, pne_acc_ukmo, pla_acc_ukmo, pno_acc_ukmo, pal_acc_ukmo))
ukmo_racc = np.stack((pel_racc_ukmo, pne_racc_ukmo, pla_racc_ukmo, pno_racc_ukmo, pal_racc_ukmo))
ukmo_boot = np.stack((pel_boot_ukmo, pne_boot_ukmo, pla_boot_ukmo, pno_boot_ukmo))
ukmo_rboot = np.stack((pel_rboot_ukmo, pne_rboot_ukmo, pla_rboot_ukmo, pno_rboot_ukmo))

acc = np.stack((bam_acc, ecmf_acc, ncep_acc, ukmo_acc))
racc = np.stack((bam_racc, ecmf_racc, ncep_racc, ukmo_racc))
boot = np.stack((bam_boot, ecmf_boot, ncep_boot, ukmo_boot))
rboot = np.stack((bam_rboot, ecmf_rboot, ncep_rboot, ukmo_rboot))
np.savez('enso_boot', acc=acc, racc=racc, boot=boot, rboot=rboot)

print "Completed the bootstrapping"


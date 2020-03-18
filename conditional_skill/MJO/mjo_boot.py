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

mjo = np.load('mjo_bchirps.npz')
bchirps = mjo['mjoa_bchirps']

mjo = np.load('mjo_chirps.npz')
chirps = mjo['mjoa_chirps']

mjo = np.load('mjo_ukmo.npz')
ukmo = mjo['mjoa_ukmo']

mjo = np.load('mjo_ncep.npz')
ncep = mjo['mjoa_ncep']

mjo = np.load('mjo_ecmf.npz')
ecmf = mjo['mjoa_ecmf']

mjo = np.load('mjo_bam.npz')
bam = mjo['mjoa_bam']


bchirps_ano = [0]*6;chirps_ano = [0]*6;ukmo_ano = [0]*6;ncep_ano = [0]*6;ecmf_ano = [0]*6;bam_ano = [0]*6
idx1 = [8,2,4,6]
idx2 = [1,3,5,7]
for i in range(len(idx1)):
    bchirps_ano[i] = np.vstack((bchirps[idx1[i]],bchirps[idx2[i]]))
    chirps_ano[i] = np.vstack((chirps[idx1[i]],chirps[idx2[i]]))
    ukmo_ano[i] = np.nanmean(np.vstack((ukmo[idx1[i]],ukmo[idx2[i]])), axis=2)  
    ncep_ano[i] = np.nanmean(np.vstack((ncep[idx1[i]],ncep[idx2[i]])), axis=2)
    ecmf_ano[i] = np.nanmean(np.vstack((ecmf[idx1[i]],ecmf[idx2[i]])), axis=2)
    bam_ano[i] = np.nanmean(np.vstack((bam[idx1[i]],bam[idx2[i]])), axis=2)

bchirps_ano[4] = bchirps[0].copy()
chirps_ano[4] = chirps[0].copy()
ukmo_ano[4] = np.nanmean(ukmo[0], axis=2)  
ncep_ano[4] = np.nanmean(ncep[0], axis=2)  
ecmf_ano[4] = np.nanmean(ecmf[0], axis=2)  
bam_ano[4] = np.nanmean(bam[0], axis=2)  

bchirps_ano[5] = bchirps[9].copy()
chirps_ano[5] = chirps[9].copy()
ukmo_ano[5] = np.nanmean(ukmo[9], axis=2)  
ncep_ano[5] = np.nanmean(ncep[9], axis=2)  
ecmf_ano[5] = np.nanmean(ecmf[9], axis=2)  
bam_ano[5] = np.nanmean(bam[9], axis=2)  

i0 = np.arange(N0); i1 = np.arange(N1); i2 = np.arange(N2); i3 = np.arange(N3) ; i4 = np.arange(N4); i5 = np.arange(N5) 
ib0 = np.arange(NB0); ib1 = np.arange(NB1); ib2 = np.arange(NB2); ib3 = np.arange(NB3) ; ib4 = np.arange(NB4); ib5 = np.arange(NB5) 
p0_dist_ukmo = np.zeros((NUM,WEEK,LAT,LON))
p1_dist_ukmo = np.zeros((NUM,WEEK,LAT,LON))
p2_dist_ukmo = np.zeros((NUM,WEEK,LAT,LON))
p3_dist_ukmo = np.zeros((NUM,WEEK,LAT,LON))
p4_dist_ukmo = np.zeros((NUM,WEEK,LAT,LON))
p5_dist_ukmo = np.zeros((NUM,WEEK,LAT,LON))
p0_dist_ncep = np.zeros((NUM,WEEK,LAT,LON))
p1_dist_ncep = np.zeros((NUM,WEEK,LAT,LON))
p2_dist_ncep = np.zeros((NUM,WEEK,LAT,LON))
p3_dist_ncep = np.zeros((NUM,WEEK,LAT,LON))
p4_dist_ncep = np.zeros((NUM,WEEK,LAT,LON))
p5_dist_ncep = np.zeros((NUM,WEEK,LAT,LON))
p0_dist_ecmf = np.zeros((NUM,WEEK,LAT,LON))
p1_dist_ecmf = np.zeros((NUM,WEEK,LAT,LON))
p2_dist_ecmf = np.zeros((NUM,WEEK,LAT,LON))
p3_dist_ecmf = np.zeros((NUM,WEEK,LAT,LON))
p4_dist_ecmf = np.zeros((NUM,WEEK,LAT,LON))
p5_dist_ecmf = np.zeros((NUM,WEEK,LAT,LON))
p0_dist_bam = np.zeros((NUM,WEEK,LAT,LON))
p1_dist_bam = np.zeros((NUM,WEEK,LAT,LON))
p2_dist_bam = np.zeros((NUM,WEEK,LAT,LON))
p3_dist_bam = np.zeros((NUM,WEEK,LAT,LON))
p4_dist_bam = np.zeros((NUM,WEEK,LAT,LON))
p5_dist_bam = np.zeros((NUM,WEEK,LAT,LON))
for n in range(0,NUM):
    A = np.random.choice(i0, int(N0), replace=True)
    B = np.random.choice(i1, int(N1), replace=True)
    C = np.random.choice(i2, int(N2), replace=True)
    D = np.random.choice(i3, int(N3), replace=True)
    E = np.random.choice(i4, int(N4), replace=True)
    F = np.random.choice(i5, int(N5), replace=True)
    L = np.random.choice(ib0, int(NB0), replace=True)
    M = np.random.choice(ib1, int(NB1), replace=True)
    N = np.random.choice(ib2, int(NB2), replace=True)
    O = np.random.choice(ib3, int(NB3), replace=True)
    P = np.random.choice(ib4, int(NB4), replace=True)
    Q = np.random.choice(ib5, int(NB5), replace=True)
    for i in range(0,LAT):
        for j in range(0,LON):
            for w in range(0,WEEK):
	        p0_dist_ukmo[n,w,i,j] = sp.stats.pearsonr(ukmo_ano[0][A,w,i,j], chirps_ano[0][A,w,i,j])[0]
	        p1_dist_ukmo[n,w,i,j] = sp.stats.pearsonr(ukmo_ano[1][B,w,i,j], chirps_ano[1][B,w,i,j])[0]
	        p2_dist_ukmo[n,w,i,j] = sp.stats.pearsonr(ukmo_ano[2][C,w,i,j], chirps_ano[2][C,w,i,j])[0]
	        p3_dist_ukmo[n,w,i,j] = sp.stats.pearsonr(ukmo_ano[3][D,w,i,j], chirps_ano[3][D,w,i,j])[0]
	        p4_dist_ukmo[n,w,i,j] = sp.stats.pearsonr(ukmo_ano[4][E,w,i,j], chirps_ano[4][E,w,i,j])[0]
	        p5_dist_ukmo[n,w,i,j] = sp.stats.pearsonr(ukmo_ano[5][F,w,i,j], chirps_ano[5][F,w,i,j])[0]
	        p0_dist_ncep[n,w,i,j] = sp.stats.pearsonr(ncep_ano[0][A,w,i,j], chirps_ano[0][A,w,i,j])[0]
	        p1_dist_ncep[n,w,i,j] = sp.stats.pearsonr(ncep_ano[1][B,w,i,j], chirps_ano[1][B,w,i,j])[0]
	        p2_dist_ncep[n,w,i,j] = sp.stats.pearsonr(ncep_ano[2][C,w,i,j], chirps_ano[2][C,w,i,j])[0]
	        p3_dist_ncep[n,w,i,j] = sp.stats.pearsonr(ncep_ano[3][D,w,i,j], chirps_ano[3][D,w,i,j])[0]
	        p4_dist_ncep[n,w,i,j] = sp.stats.pearsonr(ncep_ano[4][E,w,i,j], chirps_ano[4][E,w,i,j])[0]
	        p5_dist_ncep[n,w,i,j] = sp.stats.pearsonr(ncep_ano[5][F,w,i,j], chirps_ano[5][F,w,i,j])[0]
	        p0_dist_ecmf[n,w,i,j] = sp.stats.pearsonr(ecmf_ano[0][A,w,i,j], chirps_ano[0][A,w,i,j])[0]
	        p1_dist_ecmf[n,w,i,j] = sp.stats.pearsonr(ecmf_ano[1][B,w,i,j], chirps_ano[1][B,w,i,j])[0]
	        p2_dist_ecmf[n,w,i,j] = sp.stats.pearsonr(ecmf_ano[2][C,w,i,j], chirps_ano[2][C,w,i,j])[0]
	        p3_dist_ecmf[n,w,i,j] = sp.stats.pearsonr(ecmf_ano[3][D,w,i,j], chirps_ano[3][D,w,i,j])[0]
	        p4_dist_ecmf[n,w,i,j] = sp.stats.pearsonr(ecmf_ano[4][E,w,i,j], chirps_ano[4][E,w,i,j])[0]
	        p5_dist_ecmf[n,w,i,j] = sp.stats.pearsonr(ecmf_ano[5][F,w,i,j], chirps_ano[5][F,w,i,j])[0]
	        p0_dist_bam[n,w,i,j] = sp.stats.pearsonr(bam_ano[0][L,w,i,j], bchirps_ano[0][L,w,i,j])[0]
	        p1_dist_bam[n,w,i,j] = sp.stats.pearsonr(bam_ano[1][M,w,i,j], bchirps_ano[1][M,w,i,j])[0]
	        p2_dist_bam[n,w,i,j] = sp.stats.pearsonr(bam_ano[2][N,w,i,j], bchirps_ano[2][N,w,i,j])[0]
	        p3_dist_bam[n,w,i,j] = sp.stats.pearsonr(bam_ano[3][O,w,i,j], bchirps_ano[3][O,w,i,j])[0]
	        p4_dist_bam[n,w,i,j] = sp.stats.pearsonr(bam_ano[4][P,w,i,j], bchirps_ano[4][P,w,i,j])[0]
	        p5_dist_bam[n,w,i,j] = sp.stats.pearsonr(bam_ano[5][Q,w,i,j], bchirps_ano[5][Q,w,i,j])[0]
    print n

#####MJO PH0#####

p0_acc_ukmo = np.zeros((WEEK,LAT,LON)); p0_boot_ukmo = np.zeros((WEEK,LAT,LON))
p1_acc_ukmo = np.zeros((WEEK,LAT,LON)); p1_boot_ukmo = np.zeros((WEEK,LAT,LON))
p2_acc_ukmo = np.zeros((WEEK,LAT,LON)); p2_boot_ukmo = np.zeros((WEEK,LAT,LON))
p3_acc_ukmo = np.zeros((WEEK,LAT,LON)); p3_boot_ukmo = np.zeros((WEEK,LAT,LON))
p4_acc_ukmo = np.zeros((WEEK,LAT,LON)); p4_boot_ukmo = np.zeros((WEEK,LAT,LON))
p5_acc_ukmo = np.zeros((WEEK,LAT,LON)); p5_boot_ukmo = np.zeros((WEEK,LAT,LON))
p0_acc_ncep = np.zeros((WEEK,LAT,LON)); p0_boot_ncep = np.zeros((WEEK,LAT,LON))
p1_acc_ncep = np.zeros((WEEK,LAT,LON)); p1_boot_ncep = np.zeros((WEEK,LAT,LON))
p2_acc_ncep = np.zeros((WEEK,LAT,LON)); p2_boot_ncep = np.zeros((WEEK,LAT,LON))
p3_acc_ncep = np.zeros((WEEK,LAT,LON)); p3_boot_ncep = np.zeros((WEEK,LAT,LON))
p4_acc_ncep = np.zeros((WEEK,LAT,LON)); p4_boot_ncep = np.zeros((WEEK,LAT,LON))
p5_acc_ncep = np.zeros((WEEK,LAT,LON)); p5_boot_ncep = np.zeros((WEEK,LAT,LON))
p0_acc_ecmf = np.zeros((WEEK,LAT,LON)); p0_boot_ecmf = np.zeros((WEEK,LAT,LON))
p1_acc_ecmf = np.zeros((WEEK,LAT,LON)); p1_boot_ecmf = np.zeros((WEEK,LAT,LON))
p2_acc_ecmf = np.zeros((WEEK,LAT,LON)); p2_boot_ecmf = np.zeros((WEEK,LAT,LON))
p3_acc_ecmf = np.zeros((WEEK,LAT,LON)); p3_boot_ecmf = np.zeros((WEEK,LAT,LON))
p4_acc_ecmf = np.zeros((WEEK,LAT,LON)); p4_boot_ecmf = np.zeros((WEEK,LAT,LON))
p5_acc_ecmf = np.zeros((WEEK,LAT,LON)); p5_boot_ecmf = np.zeros((WEEK,LAT,LON))
p0_acc_bam = np.zeros((WEEK,LAT,LON)); p0_boot_bam = np.zeros((WEEK,LAT,LON))
p1_acc_bam = np.zeros((WEEK,LAT,LON)); p1_boot_bam = np.zeros((WEEK,LAT,LON))
p2_acc_bam = np.zeros((WEEK,LAT,LON)); p2_boot_bam = np.zeros((WEEK,LAT,LON))
p3_acc_bam = np.zeros((WEEK,LAT,LON)); p3_boot_bam = np.zeros((WEEK,LAT,LON))
p4_acc_bam = np.zeros((WEEK,LAT,LON)); p4_boot_bam = np.zeros((WEEK,LAT,LON))
p5_acc_bam = np.zeros((WEEK,LAT,LON)); p5_boot_bam = np.zeros((WEEK,LAT,LON))
for i in range(0,LAT):
    for j in range(0,LON):
        for w in range(0,WEEK):
	    p0_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[0][:,w,i,j], chirps_ano[0][:,w,i,j])[0]
	    p1_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[1][:,w,i,j], chirps_ano[1][:,w,i,j])[0]
	    p2_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[2][:,w,i,j], chirps_ano[2][:,w,i,j])[0]
	    p3_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[3][:,w,i,j], chirps_ano[3][:,w,i,j])[0]
	    p4_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[4][:,w,i,j], chirps_ano[4][:,w,i,j])[0]
	    p5_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[5][:,w,i,j], chirps_ano[5][:,w,i,j])[0]
	    p0_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[0][:,w,i,j], chirps_ano[0][:,w,i,j])[0]
	    p1_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[1][:,w,i,j], chirps_ano[1][:,w,i,j])[0]
	    p2_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[2][:,w,i,j], chirps_ano[2][:,w,i,j])[0]
	    p3_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[3][:,w,i,j], chirps_ano[3][:,w,i,j])[0]
	    p4_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[4][:,w,i,j], chirps_ano[4][:,w,i,j])[0]
	    p5_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[5][:,w,i,j], chirps_ano[5][:,w,i,j])[0]
	    p0_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[0][:,w,i,j], chirps_ano[0][:,w,i,j])[0]
	    p1_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[1][:,w,i,j], chirps_ano[1][:,w,i,j])[0]
	    p2_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[2][:,w,i,j], chirps_ano[2][:,w,i,j])[0]
	    p3_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[3][:,w,i,j], chirps_ano[3][:,w,i,j])[0]
	    p4_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[4][:,w,i,j], chirps_ano[4][:,w,i,j])[0]
	    p5_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[5][:,w,i,j], chirps_ano[5][:,w,i,j])[0]
	    p0_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[0][:,w,i,j], bchirps_ano[0][:,w,i,j])[0]
	    p1_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[1][:,w,i,j], bchirps_ano[1][:,w,i,j])[0]
	    p2_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[2][:,w,i,j], bchirps_ano[2][:,w,i,j])[0]
	    p3_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[3][:,w,i,j], bchirps_ano[3][:,w,i,j])[0]
	    p4_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[4][:,w,i,j], bchirps_ano[4][:,w,i,j])[0]
	    p5_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[5][:,w,i,j], bchirps_ano[5][:,w,i,j])[0]
            small = p0_dist_ukmo[:,w,i,j][(p0_dist_ukmo[:,w,i,j]<=np.nanpercentile(p0_dist_ukmo[:,w,i,j], 95))&(p0_dist_ukmo[:,w,i,j]>=np.nanpercentile(p0_dist_ukmo[:,w,i,j], 5))]; p0_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ukmo[:,w,i,j], 5))) 
            small = p1_dist_ukmo[:,w,i,j][(p1_dist_ukmo[:,w,i,j]<=np.nanpercentile(p1_dist_ukmo[:,w,i,j], 95))&(p1_dist_ukmo[:,w,i,j]>=np.nanpercentile(p1_dist_ukmo[:,w,i,j], 5))]; p1_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ukmo[:,w,i,j], 5))) 
            small = p2_dist_ukmo[:,w,i,j][(p2_dist_ukmo[:,w,i,j]<=np.nanpercentile(p2_dist_ukmo[:,w,i,j], 95))&(p2_dist_ukmo[:,w,i,j]>=np.nanpercentile(p2_dist_ukmo[:,w,i,j], 5))]; p2_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ukmo[:,w,i,j], 5))) 
            small = p3_dist_ukmo[:,w,i,j][(p3_dist_ukmo[:,w,i,j]<=np.nanpercentile(p3_dist_ukmo[:,w,i,j], 95))&(p3_dist_ukmo[:,w,i,j]>=np.nanpercentile(p3_dist_ukmo[:,w,i,j], 5))]; p3_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ukmo[:,w,i,j], 5))) 
            small = p4_dist_ukmo[:,w,i,j][(p4_dist_ukmo[:,w,i,j]<=np.nanpercentile(p4_dist_ukmo[:,w,i,j], 95))&(p4_dist_ukmo[:,w,i,j]>=np.nanpercentile(p4_dist_ukmo[:,w,i,j], 5))]; p4_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ukmo[:,w,i,j], 5))) 
            small = p5_dist_ukmo[:,w,i,j][(p5_dist_ukmo[:,w,i,j]<=np.nanpercentile(p5_dist_ukmo[:,w,i,j], 95))&(p5_dist_ukmo[:,w,i,j]>=np.nanpercentile(p5_dist_ukmo[:,w,i,j], 5))]; p5_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ukmo[:,w,i,j], 5))) 
            small = p0_dist_ncep[:,w,i,j][(p0_dist_ncep[:,w,i,j]<=np.nanpercentile(p0_dist_ncep[:,w,i,j], 95))&(p0_dist_ncep[:,w,i,j]>=np.nanpercentile(p0_dist_ncep[:,w,i,j], 5))]; p0_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ncep[:,w,i,j], 5))) 
            small = p1_dist_ncep[:,w,i,j][(p1_dist_ncep[:,w,i,j]<=np.nanpercentile(p1_dist_ncep[:,w,i,j], 95))&(p1_dist_ncep[:,w,i,j]>=np.nanpercentile(p1_dist_ncep[:,w,i,j], 5))]; p1_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ncep[:,w,i,j], 5))) 
            small = p2_dist_ncep[:,w,i,j][(p2_dist_ncep[:,w,i,j]<=np.nanpercentile(p2_dist_ncep[:,w,i,j], 95))&(p2_dist_ncep[:,w,i,j]>=np.nanpercentile(p2_dist_ncep[:,w,i,j], 5))]; p2_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ncep[:,w,i,j], 5))) 
            small = p3_dist_ncep[:,w,i,j][(p3_dist_ncep[:,w,i,j]<=np.nanpercentile(p3_dist_ncep[:,w,i,j], 95))&(p3_dist_ncep[:,w,i,j]>=np.nanpercentile(p3_dist_ncep[:,w,i,j], 5))]; p3_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ncep[:,w,i,j], 5))) 
            small = p4_dist_ncep[:,w,i,j][(p4_dist_ncep[:,w,i,j]<=np.nanpercentile(p4_dist_ncep[:,w,i,j], 95))&(p4_dist_ncep[:,w,i,j]>=np.nanpercentile(p4_dist_ncep[:,w,i,j], 5))]; p4_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ncep[:,w,i,j], 5))) 
            small = p5_dist_ncep[:,w,i,j][(p5_dist_ncep[:,w,i,j]<=np.nanpercentile(p5_dist_ncep[:,w,i,j], 95))&(p5_dist_ncep[:,w,i,j]>=np.nanpercentile(p5_dist_ncep[:,w,i,j], 5))]; p5_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ncep[:,w,i,j], 5))) 
            small = p0_dist_ecmf[:,w,i,j][(p0_dist_ecmf[:,w,i,j]<=np.nanpercentile(p0_dist_ecmf[:,w,i,j], 95))&(p0_dist_ecmf[:,w,i,j]>=np.nanpercentile(p0_dist_ecmf[:,w,i,j], 5))]; p0_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ecmf[:,w,i,j], 5))) 
            small = p1_dist_ecmf[:,w,i,j][(p1_dist_ecmf[:,w,i,j]<=np.nanpercentile(p1_dist_ecmf[:,w,i,j], 95))&(p1_dist_ecmf[:,w,i,j]>=np.nanpercentile(p1_dist_ecmf[:,w,i,j], 5))]; p1_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ecmf[:,w,i,j], 5))) 
            small = p2_dist_ecmf[:,w,i,j][(p2_dist_ecmf[:,w,i,j]<=np.nanpercentile(p2_dist_ecmf[:,w,i,j], 95))&(p2_dist_ecmf[:,w,i,j]>=np.nanpercentile(p2_dist_ecmf[:,w,i,j], 5))]; p2_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ecmf[:,w,i,j], 5))) 
            small = p3_dist_ecmf[:,w,i,j][(p3_dist_ecmf[:,w,i,j]<=np.nanpercentile(p3_dist_ecmf[:,w,i,j], 95))&(p3_dist_ecmf[:,w,i,j]>=np.nanpercentile(p3_dist_ecmf[:,w,i,j], 5))]; p3_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ecmf[:,w,i,j], 5))) 
            small = p4_dist_ecmf[:,w,i,j][(p4_dist_ecmf[:,w,i,j]<=np.nanpercentile(p4_dist_ecmf[:,w,i,j], 95))&(p4_dist_ecmf[:,w,i,j]>=np.nanpercentile(p4_dist_ecmf[:,w,i,j], 5))]; p4_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ecmf[:,w,i,j], 5))) 
            small = p5_dist_ecmf[:,w,i,j][(p5_dist_ecmf[:,w,i,j]<=np.nanpercentile(p5_dist_ecmf[:,w,i,j], 95))&(p5_dist_ecmf[:,w,i,j]>=np.nanpercentile(p5_dist_ecmf[:,w,i,j], 5))]; p5_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_ecmf[:,w,i,j], 5))) 
            small = p0_dist_bam[:,w,i,j][(p0_dist_bam[:,w,i,j]<=np.nanpercentile(p0_dist_bam[:,w,i,j], 95))&(p0_dist_bam[:,w,i,j]>=np.nanpercentile(p0_dist_bam[:,w,i,j], 5))]; p0_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_bam[:,w,i,j], 5))) 
            small = p1_dist_bam[:,w,i,j][(p1_dist_bam[:,w,i,j]<=np.nanpercentile(p1_dist_bam[:,w,i,j], 95))&(p1_dist_bam[:,w,i,j]>=np.nanpercentile(p1_dist_bam[:,w,i,j], 5))]; p1_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_bam[:,w,i,j], 5))) 
            small = p2_dist_bam[:,w,i,j][(p2_dist_bam[:,w,i,j]<=np.nanpercentile(p2_dist_bam[:,w,i,j], 95))&(p2_dist_bam[:,w,i,j]>=np.nanpercentile(p2_dist_bam[:,w,i,j], 5))]; p2_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_bam[:,w,i,j], 5))) 
            small = p3_dist_bam[:,w,i,j][(p3_dist_bam[:,w,i,j]<=np.nanpercentile(p3_dist_bam[:,w,i,j], 95))&(p3_dist_bam[:,w,i,j]>=np.nanpercentile(p3_dist_bam[:,w,i,j], 5))]; p3_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_bam[:,w,i,j], 5))) 
            small = p4_dist_bam[:,w,i,j][(p4_dist_bam[:,w,i,j]<=np.nanpercentile(p4_dist_bam[:,w,i,j], 95))&(p4_dist_bam[:,w,i,j]>=np.nanpercentile(p4_dist_bam[:,w,i,j], 5))]; p4_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_bam[:,w,i,j], 5))) 
            small = p5_dist_bam[:,w,i,j][(p5_dist_bam[:,w,i,j]<=np.nanpercentile(p5_dist_bam[:,w,i,j], 95))&(p5_dist_bam[:,w,i,j]>=np.nanpercentile(p5_dist_bam[:,w,i,j], 5))]; p5_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p4_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p4_dist_bam[:,w,i,j], 5))) 

rg_names = ['NSA','AMZ','NDE','SESA','AND','PAT']
rg_lon_min = [-80,-67,-47,-60,-75,-75]
rg_lon_max = [-50,-47,-34,-48,-67,-60]
rg_lat_min = [0,-15,-15,-35,-40,-50]
rg_lat_max = [12,-5,-5,-22,-15,-40]
id_lon_min = [7,16,29,20,10,10]
id_lon_max = [27,29,38,28,16,20]
id_lat_min = [33,23,23,10,7,0]
id_lat_max = [41,30,30,19,23,7]

p0_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p1_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p2_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p3_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p4_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p5_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p0_rdist_ncep = np.zeros((NUM,WEEK,REG))
p1_rdist_ncep = np.zeros((NUM,WEEK,REG))
p2_rdist_ncep = np.zeros((NUM,WEEK,REG))
p3_rdist_ncep = np.zeros((NUM,WEEK,REG))
p4_rdist_ncep = np.zeros((NUM,WEEK,REG))
p5_rdist_ncep = np.zeros((NUM,WEEK,REG))
p0_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p1_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p2_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p3_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p4_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p5_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p0_rdist_bam = np.zeros((NUM,WEEK,REG))
p1_rdist_bam = np.zeros((NUM,WEEK,REG))
p2_rdist_bam = np.zeros((NUM,WEEK,REG))
p3_rdist_bam = np.zeros((NUM,WEEK,REG))
p4_rdist_bam = np.zeros((NUM,WEEK,REG))
p5_rdist_bam = np.zeros((NUM,WEEK,REG))
p0_racc_ukmo = np.zeros((WEEK,REG)); p0_rboot_ukmo = np.zeros((WEEK,REG))
p1_racc_ukmo = np.zeros((WEEK,REG)); p1_rboot_ukmo = np.zeros((WEEK,REG))
p2_racc_ukmo = np.zeros((WEEK,REG)); p2_rboot_ukmo = np.zeros((WEEK,REG))
p3_racc_ukmo = np.zeros((WEEK,REG)); p3_rboot_ukmo = np.zeros((WEEK,REG))
p4_racc_ukmo = np.zeros((WEEK,REG)); p4_rboot_ukmo = np.zeros((WEEK,REG))
p5_racc_ukmo = np.zeros((WEEK,REG)); p5_rboot_ukmo = np.zeros((WEEK,REG))
p0_racc_ncep = np.zeros((WEEK,REG)); p0_rboot_ncep = np.zeros((WEEK,REG))
p1_racc_ncep = np.zeros((WEEK,REG)); p1_rboot_ncep = np.zeros((WEEK,REG))
p2_racc_ncep = np.zeros((WEEK,REG)); p2_rboot_ncep = np.zeros((WEEK,REG))
p3_racc_ncep = np.zeros((WEEK,REG)); p3_rboot_ncep = np.zeros((WEEK,REG))
p4_racc_ncep = np.zeros((WEEK,REG)); p4_rboot_ncep = np.zeros((WEEK,REG))
p5_racc_ncep = np.zeros((WEEK,REG)); p5_rboot_ncep = np.zeros((WEEK,REG))
p0_racc_ecmf = np.zeros((WEEK,REG)); p0_rboot_ecmf = np.zeros((WEEK,REG))
p1_racc_ecmf = np.zeros((WEEK,REG)); p1_rboot_ecmf = np.zeros((WEEK,REG))
p2_racc_ecmf = np.zeros((WEEK,REG)); p2_rboot_ecmf = np.zeros((WEEK,REG))
p3_racc_ecmf = np.zeros((WEEK,REG)); p3_rboot_ecmf = np.zeros((WEEK,REG))
p4_racc_ecmf = np.zeros((WEEK,REG)); p4_rboot_ecmf = np.zeros((WEEK,REG))
p5_racc_ecmf = np.zeros((WEEK,REG)); p5_rboot_ecmf = np.zeros((WEEK,REG))
p0_racc_bam = np.zeros((WEEK,REG)); p0_rboot_bam = np.zeros((WEEK,REG))
p1_racc_bam = np.zeros((WEEK,REG)); p1_rboot_bam = np.zeros((WEEK,REG))
p2_racc_bam = np.zeros((WEEK,REG)); p2_rboot_bam = np.zeros((WEEK,REG))
p3_racc_bam = np.zeros((WEEK,REG)); p3_rboot_bam = np.zeros((WEEK,REG))
p4_racc_bam = np.zeros((WEEK,REG)); p4_rboot_bam = np.zeros((WEEK,REG))
p5_racc_bam = np.zeros((WEEK,REG)); p5_rboot_bam = np.zeros((WEEK,REG))
for r in range(len(rg_names)):
    p0_rdist_ukmo[:,:,r] = np.nanmean(p0_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p1_rdist_ukmo[:,:,r] = np.nanmean(p1_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p2_rdist_ukmo[:,:,r] = np.nanmean(p2_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p3_rdist_ukmo[:,:,r] = np.nanmean(p3_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p4_rdist_ukmo[:,:,r] = np.nanmean(p4_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p5_rdist_ukmo[:,:,r] = np.nanmean(p5_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p0_rdist_ncep[:,:,r] = np.nanmean(p0_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p1_rdist_ncep[:,:,r] = np.nanmean(p1_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p2_rdist_ncep[:,:,r] = np.nanmean(p2_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p3_rdist_ncep[:,:,r] = np.nanmean(p3_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p4_rdist_ncep[:,:,r] = np.nanmean(p4_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p5_rdist_ncep[:,:,r] = np.nanmean(p5_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p0_rdist_ecmf[:,:,r] = np.nanmean(p0_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p1_rdist_ecmf[:,:,r] = np.nanmean(p1_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p2_rdist_ecmf[:,:,r] = np.nanmean(p2_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p3_rdist_ecmf[:,:,r] = np.nanmean(p3_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p4_rdist_ecmf[:,:,r] = np.nanmean(p4_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p5_rdist_ecmf[:,:,r] = np.nanmean(p5_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p0_rdist_bam[:,:,r] = np.nanmean(p0_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p1_rdist_bam[:,:,r] = np.nanmean(p1_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p2_rdist_bam[:,:,r] = np.nanmean(p2_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p3_rdist_bam[:,:,r] = np.nanmean(p3_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p4_rdist_bam[:,:,r] = np.nanmean(p4_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p5_rdist_bam[:,:,r] = np.nanmean(p5_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    for w in range(0,WEEK):
        p0_racc_ukmo[w,r] = np.nanmean(p0_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p1_racc_ukmo[w,r] = np.nanmean(p1_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p2_racc_ukmo[w,r] = np.nanmean(p2_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p3_racc_ukmo[w,r] = np.nanmean(p3_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p4_racc_ukmo[w,r] = np.nanmean(p4_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p5_racc_ukmo[w,r] = np.nanmean(p5_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p0_racc_ncep[w,r] = np.nanmean(p0_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p1_racc_ncep[w,r] = np.nanmean(p1_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p2_racc_ncep[w,r] = np.nanmean(p2_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p3_racc_ncep[w,r] = np.nanmean(p3_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p4_racc_ncep[w,r] = np.nanmean(p4_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p5_racc_ncep[w,r] = np.nanmean(p5_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p0_racc_ecmf[w,r] = np.nanmean(p0_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p1_racc_ecmf[w,r] = np.nanmean(p1_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p2_racc_ecmf[w,r] = np.nanmean(p2_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p3_racc_ecmf[w,r] = np.nanmean(p3_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p4_racc_ecmf[w,r] = np.nanmean(p4_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p5_racc_ecmf[w,r] = np.nanmean(p5_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p0_racc_bam[w,r] = np.nanmean(p0_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p1_racc_bam[w,r] = np.nanmean(p1_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p2_racc_bam[w,r] = np.nanmean(p2_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p3_racc_bam[w,r] = np.nanmean(p3_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p4_racc_bam[w,r] = np.nanmean(p4_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p5_racc_bam[w,r] = np.nanmean(p5_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        small = p0_rdist_ukmo[:,w,r][(p0_rdist_ukmo[:,w,r]<=np.nanpercentile(p0_rdist_ukmo[:,w,r], 95))&(p0_rdist_ukmo[:,w,r]>=np.nanpercentile(p0_rdist_ukmo[:,w,r], 5))]; p0_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ukmo[:,w,r], 5))) 
        small = p1_rdist_ukmo[:,w,r][(p1_rdist_ukmo[:,w,r]<=np.nanpercentile(p1_rdist_ukmo[:,w,r], 95))&(p1_rdist_ukmo[:,w,r]>=np.nanpercentile(p1_rdist_ukmo[:,w,r], 5))]; p1_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ukmo[:,w,r], 5))) 
        small = p2_rdist_ukmo[:,w,r][(p2_rdist_ukmo[:,w,r]<=np.nanpercentile(p2_rdist_ukmo[:,w,r], 95))&(p2_rdist_ukmo[:,w,r]>=np.nanpercentile(p2_rdist_ukmo[:,w,r], 5))]; p2_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ukmo[:,w,r], 5))) 
        small = p3_rdist_ukmo[:,w,r][(p3_rdist_ukmo[:,w,r]<=np.nanpercentile(p3_rdist_ukmo[:,w,r], 95))&(p3_rdist_ukmo[:,w,r]>=np.nanpercentile(p3_rdist_ukmo[:,w,r], 5))]; p3_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ukmo[:,w,r], 5))) 
        small = p4_rdist_ukmo[:,w,r][(p4_rdist_ukmo[:,w,r]<=np.nanpercentile(p4_rdist_ukmo[:,w,r], 95))&(p4_rdist_ukmo[:,w,r]>=np.nanpercentile(p4_rdist_ukmo[:,w,r], 5))]; p4_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ukmo[:,w,r], 5))) 
        small = p5_rdist_ukmo[:,w,r][(p5_rdist_ukmo[:,w,r]<=np.nanpercentile(p5_rdist_ukmo[:,w,r], 95))&(p5_rdist_ukmo[:,w,r]>=np.nanpercentile(p5_rdist_ukmo[:,w,r], 5))]; p5_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ukmo[:,w,r], 5))) 
        small = p0_rdist_ncep[:,w,r][(p0_rdist_ncep[:,w,r]<=np.nanpercentile(p0_rdist_ncep[:,w,r], 95))&(p0_rdist_ncep[:,w,r]>=np.nanpercentile(p0_rdist_ncep[:,w,r], 5))]; p0_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ncep[:,w,r], 5))) 
        small = p1_rdist_ncep[:,w,r][(p1_rdist_ncep[:,w,r]<=np.nanpercentile(p1_rdist_ncep[:,w,r], 95))&(p1_rdist_ncep[:,w,r]>=np.nanpercentile(p1_rdist_ncep[:,w,r], 5))]; p1_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ncep[:,w,r], 5))) 
        small = p2_rdist_ncep[:,w,r][(p2_rdist_ncep[:,w,r]<=np.nanpercentile(p2_rdist_ncep[:,w,r], 95))&(p2_rdist_ncep[:,w,r]>=np.nanpercentile(p2_rdist_ncep[:,w,r], 5))]; p2_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ncep[:,w,r], 5))) 
        small = p3_rdist_ncep[:,w,r][(p3_rdist_ncep[:,w,r]<=np.nanpercentile(p3_rdist_ncep[:,w,r], 95))&(p3_rdist_ncep[:,w,r]>=np.nanpercentile(p3_rdist_ncep[:,w,r], 5))]; p3_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ncep[:,w,r], 5))) 
        small = p4_rdist_ncep[:,w,r][(p4_rdist_ncep[:,w,r]<=np.nanpercentile(p4_rdist_ncep[:,w,r], 95))&(p4_rdist_ncep[:,w,r]>=np.nanpercentile(p4_rdist_ncep[:,w,r], 5))]; p4_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ncep[:,w,r], 5))) 
        small = p5_rdist_ncep[:,w,r][(p5_rdist_ncep[:,w,r]<=np.nanpercentile(p5_rdist_ncep[:,w,r], 95))&(p5_rdist_ncep[:,w,r]>=np.nanpercentile(p5_rdist_ncep[:,w,r], 5))]; p5_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ncep[:,w,r], 5))) 
        small = p0_rdist_ecmf[:,w,r][(p0_rdist_ecmf[:,w,r]<=np.nanpercentile(p0_rdist_ecmf[:,w,r], 95))&(p0_rdist_ecmf[:,w,r]>=np.nanpercentile(p0_rdist_ecmf[:,w,r], 5))]; p0_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ecmf[:,w,r], 5))) 
        small = p1_rdist_ecmf[:,w,r][(p1_rdist_ecmf[:,w,r]<=np.nanpercentile(p1_rdist_ecmf[:,w,r], 95))&(p1_rdist_ecmf[:,w,r]>=np.nanpercentile(p1_rdist_ecmf[:,w,r], 5))]; p1_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ecmf[:,w,r], 5))) 
        small = p2_rdist_ecmf[:,w,r][(p2_rdist_ecmf[:,w,r]<=np.nanpercentile(p2_rdist_ecmf[:,w,r], 95))&(p2_rdist_ecmf[:,w,r]>=np.nanpercentile(p2_rdist_ecmf[:,w,r], 5))]; p2_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ecmf[:,w,r], 5))) 
        small = p3_rdist_ecmf[:,w,r][(p3_rdist_ecmf[:,w,r]<=np.nanpercentile(p3_rdist_ecmf[:,w,r], 95))&(p3_rdist_ecmf[:,w,r]>=np.nanpercentile(p3_rdist_ecmf[:,w,r], 5))]; p3_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ecmf[:,w,r], 5))) 
        small = p4_rdist_ecmf[:,w,r][(p4_rdist_ecmf[:,w,r]<=np.nanpercentile(p4_rdist_ecmf[:,w,r], 95))&(p4_rdist_ecmf[:,w,r]>=np.nanpercentile(p4_rdist_ecmf[:,w,r], 5))]; p4_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ecmf[:,w,r], 5))) 
        small = p5_rdist_ecmf[:,w,r][(p5_rdist_ecmf[:,w,r]<=np.nanpercentile(p5_rdist_ecmf[:,w,r], 95))&(p5_rdist_ecmf[:,w,r]>=np.nanpercentile(p5_rdist_ecmf[:,w,r], 5))]; p5_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_ecmf[:,w,r], 5))) 
        small = p0_rdist_bam[:,w,r][(p0_rdist_bam[:,w,r]<=np.nanpercentile(p0_rdist_bam[:,w,r], 95))&(p0_rdist_bam[:,w,r]>=np.nanpercentile(p0_rdist_bam[:,w,r], 5))]; p0_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_bam[:,w,r], 5))) 
        small = p1_rdist_bam[:,w,r][(p1_rdist_bam[:,w,r]<=np.nanpercentile(p1_rdist_bam[:,w,r], 95))&(p1_rdist_bam[:,w,r]>=np.nanpercentile(p1_rdist_bam[:,w,r], 5))]; p1_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_bam[:,w,r], 5))) 
        small = p2_rdist_bam[:,w,r][(p2_rdist_bam[:,w,r]<=np.nanpercentile(p2_rdist_bam[:,w,r], 95))&(p2_rdist_bam[:,w,r]>=np.nanpercentile(p2_rdist_bam[:,w,r], 5))]; p2_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_bam[:,w,r], 5))) 
        small = p3_rdist_bam[:,w,r][(p3_rdist_bam[:,w,r]<=np.nanpercentile(p3_rdist_bam[:,w,r], 95))&(p3_rdist_bam[:,w,r]>=np.nanpercentile(p3_rdist_bam[:,w,r], 5))]; p3_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_bam[:,w,r], 5))) 
        small = p4_rdist_bam[:,w,r][(p4_rdist_bam[:,w,r]<=np.nanpercentile(p4_rdist_bam[:,w,r], 95))&(p4_rdist_bam[:,w,r]>=np.nanpercentile(p4_rdist_bam[:,w,r], 5))]; p4_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_bam[:,w,r], 5))) 
        small = p5_rdist_bam[:,w,r][(p5_rdist_bam[:,w,r]<=np.nanpercentile(p5_rdist_bam[:,w,r], 95))&(p5_rdist_bam[:,w,r]>=np.nanpercentile(p5_rdist_bam[:,w,r], 5))]; p5_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p4_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p4_rdist_bam[:,w,r], 5))) 

ukmo_acc = np.stack((p0_acc_ukmo, p1_acc_ukmo, p2_acc_ukmo, p3_acc_ukmo, p4_acc_ukmo, p5_acc_ukmo))
ukmo_boot = np.stack((p0_boot_ukmo, p1_boot_ukmo, p2_boot_ukmo, p3_boot_ukmo, p4_boot_ukmo, p5_boot_ukmo))
ncep_acc = np.stack((p0_acc_ncep, p1_acc_ncep, p2_acc_ncep, p3_acc_ncep, p4_acc_ncep, p5_acc_ncep))
ncep_boot = np.stack((p0_boot_ncep, p1_boot_ncep, p2_boot_ncep, p3_boot_ncep, p4_boot_ncep, p5_boot_ncep))
ecmf_acc = np.stack((p0_acc_ecmf, p1_acc_ecmf, p2_acc_ecmf, p3_acc_ecmf, p4_acc_ecmf, p5_acc_ecmf))
ecmf_boot = np.stack((p0_boot_ecmf, p1_boot_ecmf, p2_boot_ecmf, p3_boot_ecmf, p4_boot_ecmf, p5_boot_ecmf))
bam_acc = np.stack((p0_acc_bam, p1_acc_bam, p2_acc_bam, p3_acc_bam, p4_acc_bam, p5_acc_bam))
bam_boot = np.stack((p0_boot_bam, p1_boot_bam, p2_boot_bam, p3_boot_bam, p4_boot_bam, p5_boot_bam))
ukmo_racc = np.stack((p0_racc_ukmo, p1_racc_ukmo, p2_racc_ukmo, p3_racc_ukmo, p4_racc_ukmo, p5_racc_ukmo))
ukmo_rboot = np.stack((p0_rboot_ukmo, p1_rboot_ukmo, p2_rboot_ukmo, p3_rboot_ukmo, p4_rboot_ukmo, p5_rboot_ukmo))
ncep_racc = np.stack((p0_racc_ncep, p1_racc_ncep, p2_racc_ncep, p3_racc_ncep, p4_racc_ncep, p5_racc_ncep))
ncep_rboot = np.stack((p0_rboot_ncep, p1_rboot_ncep, p2_rboot_ncep, p3_rboot_ncep, p4_rboot_ncep, p5_rboot_ncep))
ecmf_racc = np.stack((p0_racc_ecmf, p1_racc_ecmf, p2_racc_ecmf, p3_racc_ecmf, p4_racc_ecmf, p5_racc_ecmf))
ecmf_rboot = np.stack((p0_rboot_ecmf, p1_rboot_ecmf, p2_rboot_ecmf, p3_rboot_ecmf, p4_rboot_ecmf, p5_rboot_ecmf))
bam_racc = np.stack((p0_racc_bam, p1_racc_bam, p2_racc_bam, p3_racc_bam, p4_racc_bam, p5_racc_bam))
bam_rboot = np.stack((p0_rboot_bam, p1_rboot_bam, p2_rboot_bam, p3_rboot_bam, p4_rboot_bam, p5_rboot_bam))

acc = np.stack((bam_acc, ecmf_acc, ncep_acc, ukmo_acc))
racc = np.stack((bam_racc, ecmf_racc, ncep_racc, ukmo_racc))
boot = np.stack((bam_boot, ecmf_boot, ncep_boot, ukmo_boot))
rboot = np.stack((bam_rboot, ecmf_rboot, ncep_rboot, ukmo_rboot))
np.savez('mjo_boot_ph0', acc=acc, racc=racc, boot=boot, rboot=rboot)

'''
#####ALL YEARS#####

p0_acc_ukmo = np.zeros((WEEK,LAT,LON)); p0_boot_ukmo = np.zeros((WEEK,LAT,LON))
p1_acc_ukmo = np.zeros((WEEK,LAT,LON)); p1_boot_ukmo = np.zeros((WEEK,LAT,LON))
p2_acc_ukmo = np.zeros((WEEK,LAT,LON)); p2_boot_ukmo = np.zeros((WEEK,LAT,LON))
p3_acc_ukmo = np.zeros((WEEK,LAT,LON)); p3_boot_ukmo = np.zeros((WEEK,LAT,LON))
p4_acc_ukmo = np.zeros((WEEK,LAT,LON)); p4_boot_ukmo = np.zeros((WEEK,LAT,LON))
p5_acc_ukmo = np.zeros((WEEK,LAT,LON)); p5_boot_ukmo = np.zeros((WEEK,LAT,LON))
p0_acc_ncep = np.zeros((WEEK,LAT,LON)); p0_boot_ncep = np.zeros((WEEK,LAT,LON))
p1_acc_ncep = np.zeros((WEEK,LAT,LON)); p1_boot_ncep = np.zeros((WEEK,LAT,LON))
p2_acc_ncep = np.zeros((WEEK,LAT,LON)); p2_boot_ncep = np.zeros((WEEK,LAT,LON))
p3_acc_ncep = np.zeros((WEEK,LAT,LON)); p3_boot_ncep = np.zeros((WEEK,LAT,LON))
p4_acc_ncep = np.zeros((WEEK,LAT,LON)); p4_boot_ncep = np.zeros((WEEK,LAT,LON))
p5_acc_ncep = np.zeros((WEEK,LAT,LON)); p5_boot_ncep = np.zeros((WEEK,LAT,LON))
p0_acc_ecmf = np.zeros((WEEK,LAT,LON)); p0_boot_ecmf = np.zeros((WEEK,LAT,LON))
p1_acc_ecmf = np.zeros((WEEK,LAT,LON)); p1_boot_ecmf = np.zeros((WEEK,LAT,LON))
p2_acc_ecmf = np.zeros((WEEK,LAT,LON)); p2_boot_ecmf = np.zeros((WEEK,LAT,LON))
p3_acc_ecmf = np.zeros((WEEK,LAT,LON)); p3_boot_ecmf = np.zeros((WEEK,LAT,LON))
p4_acc_ecmf = np.zeros((WEEK,LAT,LON)); p4_boot_ecmf = np.zeros((WEEK,LAT,LON))
p5_acc_ecmf = np.zeros((WEEK,LAT,LON)); p5_boot_ecmf = np.zeros((WEEK,LAT,LON))
p0_acc_bam = np.zeros((WEEK,LAT,LON)); p0_boot_bam = np.zeros((WEEK,LAT,LON))
p1_acc_bam = np.zeros((WEEK,LAT,LON)); p1_boot_bam = np.zeros((WEEK,LAT,LON))
p2_acc_bam = np.zeros((WEEK,LAT,LON)); p2_boot_bam = np.zeros((WEEK,LAT,LON))
p3_acc_bam = np.zeros((WEEK,LAT,LON)); p3_boot_bam = np.zeros((WEEK,LAT,LON))
p4_acc_bam = np.zeros((WEEK,LAT,LON)); p4_boot_bam = np.zeros((WEEK,LAT,LON))
p5_acc_bam = np.zeros((WEEK,LAT,LON)); p5_boot_bam = np.zeros((WEEK,LAT,LON))
for i in range(0,LAT):
    for j in range(0,LON):
        for w in range(0,WEEK):
	    p0_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[0][:,w,i,j], chirps_ano[0][:,w,i,j])[0]
	    p1_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[1][:,w,i,j], chirps_ano[1][:,w,i,j])[0]
	    p2_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[2][:,w,i,j], chirps_ano[2][:,w,i,j])[0]
	    p3_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[3][:,w,i,j], chirps_ano[3][:,w,i,j])[0]
	    p4_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[4][:,w,i,j], chirps_ano[4][:,w,i,j])[0]
	    p5_acc_ukmo[w,i,j] = sp.stats.pearsonr(ukmo_ano[5][:,w,i,j], chirps_ano[5][:,w,i,j])[0]
	    p0_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[0][:,w,i,j], chirps_ano[0][:,w,i,j])[0]
	    p1_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[1][:,w,i,j], chirps_ano[1][:,w,i,j])[0]
	    p2_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[2][:,w,i,j], chirps_ano[2][:,w,i,j])[0]
	    p3_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[3][:,w,i,j], chirps_ano[3][:,w,i,j])[0]
	    p4_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[4][:,w,i,j], chirps_ano[4][:,w,i,j])[0]
	    p5_acc_ncep[w,i,j] = sp.stats.pearsonr(ncep_ano[5][:,w,i,j], chirps_ano[5][:,w,i,j])[0]
	    p0_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[0][:,w,i,j], chirps_ano[0][:,w,i,j])[0]
	    p1_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[1][:,w,i,j], chirps_ano[1][:,w,i,j])[0]
	    p2_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[2][:,w,i,j], chirps_ano[2][:,w,i,j])[0]
	    p3_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[3][:,w,i,j], chirps_ano[3][:,w,i,j])[0]
	    p4_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[4][:,w,i,j], chirps_ano[4][:,w,i,j])[0]
	    p5_acc_ecmf[w,i,j] = sp.stats.pearsonr(ecmf_ano[5][:,w,i,j], chirps_ano[5][:,w,i,j])[0]
	    p0_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[0][:,w,i,j], bchirps_ano[0][:,w,i,j])[0]
	    p1_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[1][:,w,i,j], bchirps_ano[1][:,w,i,j])[0]
	    p2_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[2][:,w,i,j], bchirps_ano[2][:,w,i,j])[0]
	    p3_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[3][:,w,i,j], bchirps_ano[3][:,w,i,j])[0]
	    p4_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[4][:,w,i,j], bchirps_ano[4][:,w,i,j])[0]
	    p5_acc_bam[w,i,j] = sp.stats.pearsonr(bam_ano[5][:,w,i,j], bchirps_ano[5][:,w,i,j])[0]
            small = p0_dist_ukmo[:,w,i,j][(p0_dist_ukmo[:,w,i,j]<=np.nanpercentile(p0_dist_ukmo[:,w,i,j], 95))&(p0_dist_ukmo[:,w,i,j]>=np.nanpercentile(p0_dist_ukmo[:,w,i,j], 5))]; p0_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ukmo[:,w,i,j], 5))) 
            small = p1_dist_ukmo[:,w,i,j][(p1_dist_ukmo[:,w,i,j]<=np.nanpercentile(p1_dist_ukmo[:,w,i,j], 95))&(p1_dist_ukmo[:,w,i,j]>=np.nanpercentile(p1_dist_ukmo[:,w,i,j], 5))]; p1_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ukmo[:,w,i,j], 5))) 
            small = p2_dist_ukmo[:,w,i,j][(p2_dist_ukmo[:,w,i,j]<=np.nanpercentile(p2_dist_ukmo[:,w,i,j], 95))&(p2_dist_ukmo[:,w,i,j]>=np.nanpercentile(p2_dist_ukmo[:,w,i,j], 5))]; p2_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ukmo[:,w,i,j], 5))) 
            small = p3_dist_ukmo[:,w,i,j][(p3_dist_ukmo[:,w,i,j]<=np.nanpercentile(p3_dist_ukmo[:,w,i,j], 95))&(p3_dist_ukmo[:,w,i,j]>=np.nanpercentile(p3_dist_ukmo[:,w,i,j], 5))]; p3_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ukmo[:,w,i,j], 5))) 
            small = p4_dist_ukmo[:,w,i,j][(p4_dist_ukmo[:,w,i,j]<=np.nanpercentile(p4_dist_ukmo[:,w,i,j], 95))&(p4_dist_ukmo[:,w,i,j]>=np.nanpercentile(p4_dist_ukmo[:,w,i,j], 5))]; p4_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ukmo[:,w,i,j], 5))) 
            small = p5_dist_ukmo[:,w,i,j][(p5_dist_ukmo[:,w,i,j]<=np.nanpercentile(p5_dist_ukmo[:,w,i,j], 95))&(p5_dist_ukmo[:,w,i,j]>=np.nanpercentile(p5_dist_ukmo[:,w,i,j], 5))]; p5_boot_ukmo[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ukmo[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ukmo[:,w,i,j], 5))) 
            small = p0_dist_ncep[:,w,i,j][(p0_dist_ncep[:,w,i,j]<=np.nanpercentile(p0_dist_ncep[:,w,i,j], 95))&(p0_dist_ncep[:,w,i,j]>=np.nanpercentile(p0_dist_ncep[:,w,i,j], 5))]; p0_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ncep[:,w,i,j], 5))) 
            small = p1_dist_ncep[:,w,i,j][(p1_dist_ncep[:,w,i,j]<=np.nanpercentile(p1_dist_ncep[:,w,i,j], 95))&(p1_dist_ncep[:,w,i,j]>=np.nanpercentile(p1_dist_ncep[:,w,i,j], 5))]; p1_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ncep[:,w,i,j], 5))) 
            small = p2_dist_ncep[:,w,i,j][(p2_dist_ncep[:,w,i,j]<=np.nanpercentile(p2_dist_ncep[:,w,i,j], 95))&(p2_dist_ncep[:,w,i,j]>=np.nanpercentile(p2_dist_ncep[:,w,i,j], 5))]; p2_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ncep[:,w,i,j], 5))) 
            small = p3_dist_ncep[:,w,i,j][(p3_dist_ncep[:,w,i,j]<=np.nanpercentile(p3_dist_ncep[:,w,i,j], 95))&(p3_dist_ncep[:,w,i,j]>=np.nanpercentile(p3_dist_ncep[:,w,i,j], 5))]; p3_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ncep[:,w,i,j], 5))) 
            small = p4_dist_ncep[:,w,i,j][(p4_dist_ncep[:,w,i,j]<=np.nanpercentile(p4_dist_ncep[:,w,i,j], 95))&(p4_dist_ncep[:,w,i,j]>=np.nanpercentile(p4_dist_ncep[:,w,i,j], 5))]; p4_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ncep[:,w,i,j], 5))) 
            small = p5_dist_ncep[:,w,i,j][(p5_dist_ncep[:,w,i,j]<=np.nanpercentile(p5_dist_ncep[:,w,i,j], 95))&(p5_dist_ncep[:,w,i,j]>=np.nanpercentile(p5_dist_ncep[:,w,i,j], 5))]; p5_boot_ncep[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ncep[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ncep[:,w,i,j], 5))) 
            small = p0_dist_ecmf[:,w,i,j][(p0_dist_ecmf[:,w,i,j]<=np.nanpercentile(p0_dist_ecmf[:,w,i,j], 95))&(p0_dist_ecmf[:,w,i,j]>=np.nanpercentile(p0_dist_ecmf[:,w,i,j], 5))]; p0_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ecmf[:,w,i,j], 5))) 
            small = p1_dist_ecmf[:,w,i,j][(p1_dist_ecmf[:,w,i,j]<=np.nanpercentile(p1_dist_ecmf[:,w,i,j], 95))&(p1_dist_ecmf[:,w,i,j]>=np.nanpercentile(p1_dist_ecmf[:,w,i,j], 5))]; p1_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ecmf[:,w,i,j], 5))) 
            small = p2_dist_ecmf[:,w,i,j][(p2_dist_ecmf[:,w,i,j]<=np.nanpercentile(p2_dist_ecmf[:,w,i,j], 95))&(p2_dist_ecmf[:,w,i,j]>=np.nanpercentile(p2_dist_ecmf[:,w,i,j], 5))]; p2_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ecmf[:,w,i,j], 5))) 
            small = p3_dist_ecmf[:,w,i,j][(p3_dist_ecmf[:,w,i,j]<=np.nanpercentile(p3_dist_ecmf[:,w,i,j], 95))&(p3_dist_ecmf[:,w,i,j]>=np.nanpercentile(p3_dist_ecmf[:,w,i,j], 5))]; p3_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ecmf[:,w,i,j], 5))) 
            small = p4_dist_ecmf[:,w,i,j][(p4_dist_ecmf[:,w,i,j]<=np.nanpercentile(p4_dist_ecmf[:,w,i,j], 95))&(p4_dist_ecmf[:,w,i,j]>=np.nanpercentile(p4_dist_ecmf[:,w,i,j], 5))]; p4_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ecmf[:,w,i,j], 5))) 
            small = p5_dist_ecmf[:,w,i,j][(p5_dist_ecmf[:,w,i,j]<=np.nanpercentile(p5_dist_ecmf[:,w,i,j], 95))&(p5_dist_ecmf[:,w,i,j]>=np.nanpercentile(p5_dist_ecmf[:,w,i,j], 5))]; p5_boot_ecmf[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_ecmf[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_ecmf[:,w,i,j], 5))) 
            small = p0_dist_bam[:,w,i,j][(p0_dist_bam[:,w,i,j]<=np.nanpercentile(p0_dist_bam[:,w,i,j], 95))&(p0_dist_bam[:,w,i,j]>=np.nanpercentile(p0_dist_bam[:,w,i,j], 5))]; p0_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_bam[:,w,i,j], 5))) 
            small = p1_dist_bam[:,w,i,j][(p1_dist_bam[:,w,i,j]<=np.nanpercentile(p1_dist_bam[:,w,i,j], 95))&(p1_dist_bam[:,w,i,j]>=np.nanpercentile(p1_dist_bam[:,w,i,j], 5))]; p1_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_bam[:,w,i,j], 5))) 
            small = p2_dist_bam[:,w,i,j][(p2_dist_bam[:,w,i,j]<=np.nanpercentile(p2_dist_bam[:,w,i,j], 95))&(p2_dist_bam[:,w,i,j]>=np.nanpercentile(p2_dist_bam[:,w,i,j], 5))]; p2_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_bam[:,w,i,j], 5))) 
            small = p3_dist_bam[:,w,i,j][(p3_dist_bam[:,w,i,j]<=np.nanpercentile(p3_dist_bam[:,w,i,j], 95))&(p3_dist_bam[:,w,i,j]>=np.nanpercentile(p3_dist_bam[:,w,i,j], 5))]; p3_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_bam[:,w,i,j], 5))) 
            small = p4_dist_bam[:,w,i,j][(p4_dist_bam[:,w,i,j]<=np.nanpercentile(p4_dist_bam[:,w,i,j], 95))&(p4_dist_bam[:,w,i,j]>=np.nanpercentile(p4_dist_bam[:,w,i,j], 5))]; p4_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_bam[:,w,i,j], 5))) 
            small = p5_dist_bam[:,w,i,j][(p5_dist_bam[:,w,i,j]<=np.nanpercentile(p5_dist_bam[:,w,i,j], 95))&(p5_dist_bam[:,w,i,j]>=np.nanpercentile(p5_dist_bam[:,w,i,j], 5))]; p5_boot_bam[w,i,j] = np.count_nonzero(((small)>np.nanpercentile(p5_dist_bam[:,w,i,j], 95)) | ((small)<np.nanpercentile(p5_dist_bam[:,w,i,j], 5))) 

rg_names = ['NSA','AMZ','NDE','SESA','AND','PAT']
rg_lon_min = [-80,-67,-47,-60,-75,-75]
rg_lon_max = [-50,-47,-34,-48,-67,-60]
rg_lat_min = [0,-15,-15,-35,-40,-50]
rg_lat_max = [12,-5,-5,-22,-15,-40]
id_lon_min = [7,16,29,20,10,10]
id_lon_max = [27,29,38,28,16,20]
id_lat_min = [33,23,23,10,7,0]
id_lat_max = [41,30,30,19,23,7]

p0_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p1_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p2_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p3_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p4_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p5_rdist_ukmo = np.zeros((NUM,WEEK,REG))
p0_rdist_ncep = np.zeros((NUM,WEEK,REG))
p1_rdist_ncep = np.zeros((NUM,WEEK,REG))
p2_rdist_ncep = np.zeros((NUM,WEEK,REG))
p3_rdist_ncep = np.zeros((NUM,WEEK,REG))
p4_rdist_ncep = np.zeros((NUM,WEEK,REG))
p5_rdist_ncep = np.zeros((NUM,WEEK,REG))
p0_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p1_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p2_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p3_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p4_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p5_rdist_ecmf = np.zeros((NUM,WEEK,REG))
p0_rdist_bam = np.zeros((NUM,WEEK,REG))
p1_rdist_bam = np.zeros((NUM,WEEK,REG))
p2_rdist_bam = np.zeros((NUM,WEEK,REG))
p3_rdist_bam = np.zeros((NUM,WEEK,REG))
p4_rdist_bam = np.zeros((NUM,WEEK,REG))
p5_rdist_bam = np.zeros((NUM,WEEK,REG))
p0_racc_ukmo = np.zeros((WEEK,REG)); p0_rboot_ukmo = np.zeros((WEEK,REG))
p1_racc_ukmo = np.zeros((WEEK,REG)); p1_rboot_ukmo = np.zeros((WEEK,REG))
p2_racc_ukmo = np.zeros((WEEK,REG)); p2_rboot_ukmo = np.zeros((WEEK,REG))
p3_racc_ukmo = np.zeros((WEEK,REG)); p3_rboot_ukmo = np.zeros((WEEK,REG))
p4_racc_ukmo = np.zeros((WEEK,REG)); p4_rboot_ukmo = np.zeros((WEEK,REG))
p5_racc_ukmo = np.zeros((WEEK,REG)); p5_rboot_ukmo = np.zeros((WEEK,REG))
p0_racc_ncep = np.zeros((WEEK,REG)); p0_rboot_ncep = np.zeros((WEEK,REG))
p1_racc_ncep = np.zeros((WEEK,REG)); p1_rboot_ncep = np.zeros((WEEK,REG))
p2_racc_ncep = np.zeros((WEEK,REG)); p2_rboot_ncep = np.zeros((WEEK,REG))
p3_racc_ncep = np.zeros((WEEK,REG)); p3_rboot_ncep = np.zeros((WEEK,REG))
p4_racc_ncep = np.zeros((WEEK,REG)); p4_rboot_ncep = np.zeros((WEEK,REG))
p5_racc_ncep = np.zeros((WEEK,REG)); p5_rboot_ncep = np.zeros((WEEK,REG))
p0_racc_ecmf = np.zeros((WEEK,REG)); p0_rboot_ecmf = np.zeros((WEEK,REG))
p1_racc_ecmf = np.zeros((WEEK,REG)); p1_rboot_ecmf = np.zeros((WEEK,REG))
p2_racc_ecmf = np.zeros((WEEK,REG)); p2_rboot_ecmf = np.zeros((WEEK,REG))
p3_racc_ecmf = np.zeros((WEEK,REG)); p3_rboot_ecmf = np.zeros((WEEK,REG))
p4_racc_ecmf = np.zeros((WEEK,REG)); p4_rboot_ecmf = np.zeros((WEEK,REG))
p5_racc_ecmf = np.zeros((WEEK,REG)); p5_rboot_ecmf = np.zeros((WEEK,REG))
p0_racc_bam = np.zeros((WEEK,REG)); p0_rboot_bam = np.zeros((WEEK,REG))
p1_racc_bam = np.zeros((WEEK,REG)); p1_rboot_bam = np.zeros((WEEK,REG))
p2_racc_bam = np.zeros((WEEK,REG)); p2_rboot_bam = np.zeros((WEEK,REG))
p3_racc_bam = np.zeros((WEEK,REG)); p3_rboot_bam = np.zeros((WEEK,REG))
p4_racc_bam = np.zeros((WEEK,REG)); p4_rboot_bam = np.zeros((WEEK,REG))
p5_racc_bam = np.zeros((WEEK,REG)); p5_rboot_bam = np.zeros((WEEK,REG))
for r in range(len(rg_names)):
    p0_rdist_ukmo[:,:,r] = np.nanmean(p0_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p1_rdist_ukmo[:,:,r] = np.nanmean(p1_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p2_rdist_ukmo[:,:,r] = np.nanmean(p2_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p3_rdist_ukmo[:,:,r] = np.nanmean(p3_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p4_rdist_ukmo[:,:,r] = np.nanmean(p4_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p5_rdist_ukmo[:,:,r] = np.nanmean(p5_dist_ukmo[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p0_rdist_ncep[:,:,r] = np.nanmean(p0_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p1_rdist_ncep[:,:,r] = np.nanmean(p1_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p2_rdist_ncep[:,:,r] = np.nanmean(p2_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p3_rdist_ncep[:,:,r] = np.nanmean(p3_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p4_rdist_ncep[:,:,r] = np.nanmean(p4_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p5_rdist_ncep[:,:,r] = np.nanmean(p5_dist_ncep[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p0_rdist_ecmf[:,:,r] = np.nanmean(p0_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p1_rdist_ecmf[:,:,r] = np.nanmean(p1_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p2_rdist_ecmf[:,:,r] = np.nanmean(p2_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p3_rdist_ecmf[:,:,r] = np.nanmean(p3_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p4_rdist_ecmf[:,:,r] = np.nanmean(p4_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p5_rdist_ecmf[:,:,r] = np.nanmean(p5_dist_ecmf[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p0_rdist_bam[:,:,r] = np.nanmean(p0_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p1_rdist_bam[:,:,r] = np.nanmean(p1_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p2_rdist_bam[:,:,r] = np.nanmean(p2_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p3_rdist_bam[:,:,r] = np.nanmean(p3_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p4_rdist_bam[:,:,r] = np.nanmean(p4_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    p5_rdist_bam[:,:,r] = np.nanmean(p5_dist_bam[:,:,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]], axis=(2,3))
    for w in range(0,WEEK):
        p0_racc_ukmo[w,r] = np.nanmean(p0_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p1_racc_ukmo[w,r] = np.nanmean(p1_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p2_racc_ukmo[w,r] = np.nanmean(p2_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p3_racc_ukmo[w,r] = np.nanmean(p3_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p4_racc_ukmo[w,r] = np.nanmean(p4_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p5_racc_ukmo[w,r] = np.nanmean(p5_acc_ukmo[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p0_racc_ncep[w,r] = np.nanmean(p0_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p1_racc_ncep[w,r] = np.nanmean(p1_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p2_racc_ncep[w,r] = np.nanmean(p2_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p3_racc_ncep[w,r] = np.nanmean(p3_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p4_racc_ncep[w,r] = np.nanmean(p4_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p5_racc_ncep[w,r] = np.nanmean(p5_acc_ncep[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p0_racc_ecmf[w,r] = np.nanmean(p0_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p1_racc_ecmf[w,r] = np.nanmean(p1_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p2_racc_ecmf[w,r] = np.nanmean(p2_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p3_racc_ecmf[w,r] = np.nanmean(p3_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p4_racc_ecmf[w,r] = np.nanmean(p4_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p5_racc_ecmf[w,r] = np.nanmean(p5_acc_ecmf[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p0_racc_bam[w,r] = np.nanmean(p0_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p1_racc_bam[w,r] = np.nanmean(p1_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p2_racc_bam[w,r] = np.nanmean(p2_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p3_racc_bam[w,r] = np.nanmean(p3_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p4_racc_bam[w,r] = np.nanmean(p4_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        p5_racc_bam[w,r] = np.nanmean(p5_acc_bam[w,id_lat_min[r]:id_lat_max[r],id_lon_min[r]:id_lon_max[r]])
        small = p0_rdist_ukmo[:,w,r][(p0_rdist_ukmo[:,w,r]<=np.nanpercentile(p0_rdist_ukmo[:,w,r], 95))&(p0_rdist_ukmo[:,w,r]>=np.nanpercentile(p0_rdist_ukmo[:,w,r], 5))]; p0_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ukmo[:,w,r], 5))) 
        small = p1_rdist_ukmo[:,w,r][(p1_rdist_ukmo[:,w,r]<=np.nanpercentile(p1_rdist_ukmo[:,w,r], 95))&(p1_rdist_ukmo[:,w,r]>=np.nanpercentile(p1_rdist_ukmo[:,w,r], 5))]; p1_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ukmo[:,w,r], 5))) 
        small = p2_rdist_ukmo[:,w,r][(p2_rdist_ukmo[:,w,r]<=np.nanpercentile(p2_rdist_ukmo[:,w,r], 95))&(p2_rdist_ukmo[:,w,r]>=np.nanpercentile(p2_rdist_ukmo[:,w,r], 5))]; p2_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ukmo[:,w,r], 5))) 
        small = p3_rdist_ukmo[:,w,r][(p3_rdist_ukmo[:,w,r]<=np.nanpercentile(p3_rdist_ukmo[:,w,r], 95))&(p3_rdist_ukmo[:,w,r]>=np.nanpercentile(p3_rdist_ukmo[:,w,r], 5))]; p3_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ukmo[:,w,r], 5))) 
        small = p4_rdist_ukmo[:,w,r][(p4_rdist_ukmo[:,w,r]<=np.nanpercentile(p4_rdist_ukmo[:,w,r], 95))&(p4_rdist_ukmo[:,w,r]>=np.nanpercentile(p4_rdist_ukmo[:,w,r], 5))]; p4_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ukmo[:,w,r], 5))) 
        small = p5_rdist_ukmo[:,w,r][(p5_rdist_ukmo[:,w,r]<=np.nanpercentile(p5_rdist_ukmo[:,w,r], 95))&(p5_rdist_ukmo[:,w,r]>=np.nanpercentile(p5_rdist_ukmo[:,w,r], 5))]; p5_rboot_ukmo[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ukmo[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ukmo[:,w,r], 5))) 
        small = p0_rdist_ncep[:,w,r][(p0_rdist_ncep[:,w,r]<=np.nanpercentile(p0_rdist_ncep[:,w,r], 95))&(p0_rdist_ncep[:,w,r]>=np.nanpercentile(p0_rdist_ncep[:,w,r], 5))]; p0_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ncep[:,w,r], 5))) 
        small = p1_rdist_ncep[:,w,r][(p1_rdist_ncep[:,w,r]<=np.nanpercentile(p1_rdist_ncep[:,w,r], 95))&(p1_rdist_ncep[:,w,r]>=np.nanpercentile(p1_rdist_ncep[:,w,r], 5))]; p1_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ncep[:,w,r], 5))) 
        small = p2_rdist_ncep[:,w,r][(p2_rdist_ncep[:,w,r]<=np.nanpercentile(p2_rdist_ncep[:,w,r], 95))&(p2_rdist_ncep[:,w,r]>=np.nanpercentile(p2_rdist_ncep[:,w,r], 5))]; p2_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ncep[:,w,r], 5))) 
        small = p3_rdist_ncep[:,w,r][(p3_rdist_ncep[:,w,r]<=np.nanpercentile(p3_rdist_ncep[:,w,r], 95))&(p3_rdist_ncep[:,w,r]>=np.nanpercentile(p3_rdist_ncep[:,w,r], 5))]; p3_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ncep[:,w,r], 5))) 
        small = p4_rdist_ncep[:,w,r][(p4_rdist_ncep[:,w,r]<=np.nanpercentile(p4_rdist_ncep[:,w,r], 95))&(p4_rdist_ncep[:,w,r]>=np.nanpercentile(p4_rdist_ncep[:,w,r], 5))]; p4_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ncep[:,w,r], 5))) 
        small = p5_rdist_ncep[:,w,r][(p5_rdist_ncep[:,w,r]<=np.nanpercentile(p5_rdist_ncep[:,w,r], 95))&(p5_rdist_ncep[:,w,r]>=np.nanpercentile(p5_rdist_ncep[:,w,r], 5))]; p5_rboot_ncep[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ncep[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ncep[:,w,r], 5))) 
        small = p0_rdist_ecmf[:,w,r][(p0_rdist_ecmf[:,w,r]<=np.nanpercentile(p0_rdist_ecmf[:,w,r], 95))&(p0_rdist_ecmf[:,w,r]>=np.nanpercentile(p0_rdist_ecmf[:,w,r], 5))]; p0_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ecmf[:,w,r], 5))) 
        small = p1_rdist_ecmf[:,w,r][(p1_rdist_ecmf[:,w,r]<=np.nanpercentile(p1_rdist_ecmf[:,w,r], 95))&(p1_rdist_ecmf[:,w,r]>=np.nanpercentile(p1_rdist_ecmf[:,w,r], 5))]; p1_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ecmf[:,w,r], 5))) 
        small = p2_rdist_ecmf[:,w,r][(p2_rdist_ecmf[:,w,r]<=np.nanpercentile(p2_rdist_ecmf[:,w,r], 95))&(p2_rdist_ecmf[:,w,r]>=np.nanpercentile(p2_rdist_ecmf[:,w,r], 5))]; p2_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ecmf[:,w,r], 5))) 
        small = p3_rdist_ecmf[:,w,r][(p3_rdist_ecmf[:,w,r]<=np.nanpercentile(p3_rdist_ecmf[:,w,r], 95))&(p3_rdist_ecmf[:,w,r]>=np.nanpercentile(p3_rdist_ecmf[:,w,r], 5))]; p3_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ecmf[:,w,r], 5))) 
        small = p4_rdist_ecmf[:,w,r][(p4_rdist_ecmf[:,w,r]<=np.nanpercentile(p4_rdist_ecmf[:,w,r], 95))&(p4_rdist_ecmf[:,w,r]>=np.nanpercentile(p4_rdist_ecmf[:,w,r], 5))]; p4_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ecmf[:,w,r], 5))) 
        small = p5_rdist_ecmf[:,w,r][(p5_rdist_ecmf[:,w,r]<=np.nanpercentile(p5_rdist_ecmf[:,w,r], 95))&(p5_rdist_ecmf[:,w,r]>=np.nanpercentile(p5_rdist_ecmf[:,w,r], 5))]; p5_rboot_ecmf[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_ecmf[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_ecmf[:,w,r], 5))) 
        small = p0_rdist_bam[:,w,r][(p0_rdist_bam[:,w,r]<=np.nanpercentile(p0_rdist_bam[:,w,r], 95))&(p0_rdist_bam[:,w,r]>=np.nanpercentile(p0_rdist_bam[:,w,r], 5))]; p0_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_bam[:,w,r], 5))) 
        small = p1_rdist_bam[:,w,r][(p1_rdist_bam[:,w,r]<=np.nanpercentile(p1_rdist_bam[:,w,r], 95))&(p1_rdist_bam[:,w,r]>=np.nanpercentile(p1_rdist_bam[:,w,r], 5))]; p1_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_bam[:,w,r], 5))) 
        small = p2_rdist_bam[:,w,r][(p2_rdist_bam[:,w,r]<=np.nanpercentile(p2_rdist_bam[:,w,r], 95))&(p2_rdist_bam[:,w,r]>=np.nanpercentile(p2_rdist_bam[:,w,r], 5))]; p2_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_bam[:,w,r], 5))) 
        small = p3_rdist_bam[:,w,r][(p3_rdist_bam[:,w,r]<=np.nanpercentile(p3_rdist_bam[:,w,r], 95))&(p3_rdist_bam[:,w,r]>=np.nanpercentile(p3_rdist_bam[:,w,r], 5))]; p3_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_bam[:,w,r], 5))) 
        small = p4_rdist_bam[:,w,r][(p4_rdist_bam[:,w,r]<=np.nanpercentile(p4_rdist_bam[:,w,r], 95))&(p4_rdist_bam[:,w,r]>=np.nanpercentile(p4_rdist_bam[:,w,r], 5))]; p4_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_bam[:,w,r], 5))) 
        small = p5_rdist_bam[:,w,r][(p5_rdist_bam[:,w,r]<=np.nanpercentile(p5_rdist_bam[:,w,r], 95))&(p5_rdist_bam[:,w,r]>=np.nanpercentile(p5_rdist_bam[:,w,r], 5))]; p5_rboot_bam[w,r] = np.count_nonzero(((small)>np.nanpercentile(p5_rdist_bam[:,w,r], 95)) | ((small)<np.nanpercentile(p5_rdist_bam[:,w,r], 5))) 

ukmo_acc = np.stack((p0_acc_ukmo, p1_acc_ukmo, p2_acc_ukmo, p3_acc_ukmo, p4_acc_ukmo, p5_acc_ukmo))
ukmo_boot = np.stack((p0_boot_ukmo, p1_boot_ukmo, p2_boot_ukmo, p3_boot_ukmo, p4_boot_ukmo, p5_boot_ukmo))
ncep_acc = np.stack((p0_acc_ncep, p1_acc_ncep, p2_acc_ncep, p3_acc_ncep, p4_acc_ncep, p5_acc_ncep))
ncep_boot = np.stack((p0_boot_ncep, p1_boot_ncep, p2_boot_ncep, p3_boot_ncep, p4_boot_ncep, p5_boot_ncep))
ecmf_acc = np.stack((p0_acc_ecmf, p1_acc_ecmf, p2_acc_ecmf, p3_acc_ecmf, p4_acc_ecmf, p5_acc_ecmf))
ecmf_boot = np.stack((p0_boot_ecmf, p1_boot_ecmf, p2_boot_ecmf, p3_boot_ecmf, p4_boot_ecmf, p5_boot_ecmf))
bam_acc = np.stack((p0_acc_bam, p1_acc_bam, p2_acc_bam, p3_acc_bam, p4_acc_bam, p5_acc_bam))
bam_boot = np.stack((p0_boot_bam, p1_boot_bam, p2_boot_bam, p3_boot_bam, p4_boot_bam, p5_boot_bam))
ukmo_racc = np.stack((p0_racc_ukmo, p1_racc_ukmo, p2_racc_ukmo, p3_racc_ukmo, p4_racc_ukmo, p5_racc_ukmo))
ukmo_rboot = np.stack((p0_rboot_ukmo, p1_rboot_ukmo, p2_rboot_ukmo, p3_rboot_ukmo, p4_rboot_ukmo, p5_rboot_ukmo))
ncep_racc = np.stack((p0_racc_ncep, p1_racc_ncep, p2_racc_ncep, p3_racc_ncep, p4_racc_ncep, p5_racc_ncep))
ncep_rboot = np.stack((p0_rboot_ncep, p1_rboot_ncep, p2_rboot_ncep, p3_rboot_ncep, p4_rboot_ncep, p5_rboot_ncep))
ecmf_racc = np.stack((p0_racc_ecmf, p1_racc_ecmf, p2_racc_ecmf, p3_racc_ecmf, p4_racc_ecmf, p5_racc_ecmf))
ecmf_rboot = np.stack((p0_rboot_ecmf, p1_rboot_ecmf, p2_rboot_ecmf, p3_rboot_ecmf, p4_rboot_ecmf, p5_rboot_ecmf))
bam_racc = np.stack((p0_racc_bam, p1_racc_bam, p2_racc_bam, p3_racc_bam, p4_racc_bam, p5_racc_bam))
bam_rboot = np.stack((p0_rboot_bam, p1_rboot_bam, p2_rboot_bam, p3_rboot_bam, p4_rboot_bam, p5_rboot_bam))

acc = np.stack((bam_acc, ecmf_acc, ncep_acc, ukmo_acc))
racc = np.stack((bam_racc, ecmf_racc, ncep_racc, ukmo_racc))
boot = np.stack((bam_boot, ecmf_boot, ncep_boot, ukmo_boot))
rboot = np.stack((bam_rboot, ecmf_rboot, ncep_rboot, ukmo_rboot))
np.savez('mjo_boot_all', acc=acc, racc=racc, boot=boot, rboot=rboot)

'''
print "Completed the bootstrapping"


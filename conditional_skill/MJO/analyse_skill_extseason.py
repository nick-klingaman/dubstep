'''
Analyse skill of weekly S2S precipitation forecasts
from UKMO, NCEP and ECMWF over a given region

1) Reads in data saved over SAmerica
'save_weekly_forecasts_ukmo_ncep.py'
'save_weekly_forecasts_ecmf.py'
'save_weekly_forecasts_bam.py'

2) Computes mean precipitation, bias, anomaly correlation coefficient
and brier skill scores

3) Plots results as spatial maps

M. Young 29/06/2018
'''
from __future__ import division
import glob
import numpy as np
from netCDF4 import Dataset
from netCDF4 import MFDataset
import time as tt
from datetime import datetime, timedelta
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
import calendar
import os.path
import sys
import read_mjo_index as rmm
from scipy.stats import pearsonr
execfile('date_str.py')
execfile('grab_data.py')
execfile('read_hindcasts_for_analysis.py')
def mask_percentiles(data,percentile,arg):
  '''
  Function: 'mask_percentiles':
  Function takes in a 2d spatial forecast array and 2d spatial array of corresponding percentiles
  and masks the forecast array according to input arguments, 'arg':
  'between': mask the forecast between 2 percentiles
  'below': mask the forecast below a percentile
  'above': mask the forecast above a percentile
  '''
  data_mask= np.copy(data)
  if arg == 'between': # find between 2 percentiles
    # the method of masking below doesn't work
    #data_mask[(data_mask >= percentile[0,:,:]) & (data_mask <= percentile[1,:,:])] = np.nan
    idx = np.where((data_mask >= percentile[0,:,:]) & (data_mask <= percentile[1,:,:]))
    data_mask[idx] = np.nan
    data_mask[np.isnan(data_mask)==False] = 0
    data_mask[np.isnan(data_mask)==True] = 1
  elif arg == 'below': # find data below percentile
    data_mask[data_mask < percentile] = np.nan
    data_mask[data_mask >= percentile] = 0
    data_mask[np.isnan(data_mask)==True] = 1
  elif arg == 'above': # find data above percentile
    data_mask[data_mask > percentile] = np.nan
    data_mask[data_mask <= percentile] = 0
    data_mask[np.isnan(data_mask)==True] = 1
  return data_mask


# dir_in	= '/group_workspaces/jasmin2/ncas_climate/users/myoung02/datasets/S2S_forecasts/weekly_dubstep_style/'
dir_in = '/gws/nopw/j04/ncas_climate_vol1/users/myoung02/datasets/DUBSTEP_data/'
# dir_out	= '/home/users/myoung02/DUBSTEP_paper_results/'
dir_out	= '/home/users/myoung02/DUBSTEP_paper_results_revisions_March2020/'

dir_out2 = dir_out+'individual_panels/'
dir_prob	= dir_in+'hindcast_probabilities/'

model_ls= ['BAM','ECMWF','NCEP','UKMO']
# number of ensemble members for each forecast
ukmo_nmembers = 7
ncep_nmembers = 4
ecmf_nmembers = 11
bam_nmembers = 11
# number of lagged ensembles
ukmo_lags = 1
ncep_lags = 7
ecmf_lags = 3
bam_lags = 1
ukmo_memlag = ukmo_lags*ukmo_nmembers
ncep_memlag = ncep_lags*ncep_nmembers
ecmf_memlag = ecmf_lags*ecmf_nmembers
bam_memlag = bam_lags*bam_nmembers

# Define region for analysis over Brazil
# region	= 'Brazil'
# latlim	= [-40,20]
# lonlim	= [-90,-20]
region = 'SAmerica'
latlim = [-49.5,19.5] # [-40,20]
lonlim = [-90,-20]
lonlim2 =[-90,-30] # for plotting

# open single file to get lats and lons
nc_fid = Dataset(dir_in+region+'_CHIRPS_weekly_200102.nc', 'r')
all_lat = np.array(nc_fid.variables['latitude'][:])
all_lon = np.array(nc_fid.variables['longitude'][:])
nc_fid.close()

years	= np.arange(2000,2010+1,1)# december will always correspond to year-1
nleads	= 5	# number of lead times (in weeks) in the data
#season	= 'JJA'
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
season_ls = ['NDJFM','MJJAS']

season = season_ls[0]
if season == 'NDJFM':
  keep_mon = [11,12,1,2,3]
elif season == 'MJJAS':
  keep_mon = [5,6,7,8,9]

#for season in season_ls:
(all_chirps,all_ukmo_ens,all_ncep_ens,all_ecmf_ens,all_chirps_bam,all_bam_ens,week_mid_dates,week_mid_dates_bam) = read_hindcasts(season)
nweeks = all_chirps.shape[2]

all_ukmo_ens[all_ukmo_ens < 0] = 0
all_ncep_ens[all_ncep_ens < 0] = 0
all_ecmf_ens[all_ecmf_ens < 0] = 0
all_bam_ens[all_bam_ens < 0] = 0


# Forecast ensemble and lagged means
# compute the arithmetic mean along the specific axis, ignoring NaNs
all_ukmo = np.nanmean(all_ukmo_ens,axis=3)# (11, 5, 20, 40, 47)
all_ncep = np.nanmean(all_ncep_ens,axis=(3,4))
all_ecmf = np.nanmean(all_ecmf_ens,axis=(3,4))
all_bam = np.nanmean(all_bam_ens,axis=4)

# find climatological mean of each ensemble member
# all_ukmo_ens_clim = np.nanmean(all_ukmo_ens,axis=0)
# all_ncep_ens_clim = np.nanmean(all_ncep_ens,axis=0)
# all_ecmf_ens_clim = np.nanmean(all_ecmf_ens,axis=0)
# all_bam_ens_clim = np.nanmean(all_bam_ens,axis=0)

# average forecasts over lags only (leaving ensembles)
lag_ncep = np.nanmean(all_ncep_ens,axis=3)# (11, 5, 20, 4, 40, 47)
lag_ecmf = np.nanmean(all_ecmf_ens,axis=3)# (11, 5, 20, 11, 40, 47)

# Forecast means over years and ensembles (& lags for ncep/ecmf)
# week_ukmo_mean = np.nanmean(all_ukmo_ens,axis=(0,3))# (5, 20, 40, 47)
# week_ncep_mean = np.nanmean(all_ncep_ens,axis=(0,3,4))
# week_ecmf_mean = np.nanmean(all_ecmf_ens,axis=(0,3,4))
# week_bam_mean = np.nanmean(all_bam_ens,axis=(0,4))
# week_mean = [week_bam_mean,week_ecmf_mean,week_ncep_mean,week_ukmo_mean]

# chirps mean for each week
week_chirps_mean = np.nanmean(all_chirps,axis=0)# (5, 20, 40, 47)
week_chirps_bam_mean = np.nanmean(all_chirps_bam,axis=0)# (5, 20, 40, 47)

land_sea = np.nanmean(week_chirps_mean,axis=(0,1))
land_sea[land_sea>=0] = 1
land_sea[land_sea<0] = np.nan
# Create a dry mask (chirps seasonal clim < 1 mmd)
# this is used to mask Brier Skill Score later
# as BSS results are dodgy over dry regions
dry_mask = np.nanmean(np.copy(week_chirps_mean),axis=(0,1))
dry_mask[dry_mask <  1] = np.nan
dry_mask[dry_mask >= 1] = 1

dry_mask_bam = np.nanmean(np.copy(week_chirps_bam_mean),axis=(0,1))
dry_mask_bam[dry_mask_bam <  1] = np.nan
dry_mask_bam[dry_mask_bam >= 1] = 1

'''
Reading the RMM indices and identifying MJO amplitude and phase
'''
SY = 2000-1
EY = 2010
mweeki = (np.array([1,9,17,25]))-1 #weeks in the month
bweeki = (np.array([1,16]))-1 #weeks in a month for the bam model
TIME = 55
if season == 'NDJFM':
  MON = np.array([1,2,3,11,12])
elif season == 'MJJAS':
  MON = np.array([5,6,7,8,9])

am,ph = rmm.rmm_index(SY,EY,MON)

ma=np.zeros(((EY-SY)+1,len(MON),len(mweeki)))
mp=np.zeros(((EY-SY)+1,len(MON),len(mweeki)))
ba=np.zeros(((EY-SY)+1,len(MON),len(bweeki))) 
bp=np.zeros(((EY-SY)+1,len(MON),len(bweeki)))
for y in range((EY-SY)+1):
    for m in range(len(MON)):
        ma[y,m,:]=am[y][m][mweeki]
        mp[y,m,:]=ph[y][m][mweeki]
        ba[y,m,:]=am[y][m][bweeki]
        bp[y,m,:]=ph[y][m][bweeki]

ma = ma.reshape((((EY-SY)+1)*len(MON),len(mweeki)))[3:-2,:].reshape((TIME*len(mweeki)))
mp = mp.reshape((((EY-SY)+1)*len(MON),len(mweeki)))[3:-2,:].reshape((TIME*len(mweeki)))
ba = ba.reshape((((EY-SY)+1)*len(MON),len(bweeki)))[3:-2,:].reshape((TIME*len(bweeki)))
bp = bp.reshape((((EY-SY)+1)*len(MON),len(bweeki)))[3:-2,:].reshape((TIME*len(bweeki)))

'''
Compute the Anomaly correlation coefficient (ACC), bias and RMSE
'''
# set up empty arrays
acc = np.zeros((4,nleads,len(all_lat),len(all_lon)))
acc_pre = np.zeros((4,3,nleads,len(all_lat),len(all_lon)))
bias = np.zeros((4,nleads,len(all_lat),len(all_lon)))
week_bias = np.zeros((4,nleads,nweeks,len(all_lat),len(all_lon)))
rmse = np.zeros((4,nleads,len(all_lat),len(all_lon)))
rmse_mean = np.zeros((4,nleads,len(all_lat),len(all_lon)))

rmse_anom = np.zeros((4,nleads,len(all_lat),len(all_lon)))
rmse_anom_mean = np.zeros((4,nleads,len(all_lat),len(all_lon)))

nsamp = np.zeros((4,nleads))# count number of samples at each lead time
nsamp_week = np.zeros((4,nleads,nweeks))
nsamp_week_only = np.zeros((nleads))

nsamp_bam = np.zeros((nleads))
nsamp_week_bam = np.zeros((nleads,5,2))
nsamp_week_only_bam = np.zeros((nleads))

anom_ukmo = np.nan*np.zeros(all_ukmo.shape)
anom_ncep = np.nan*np.zeros(all_ncep.shape)
anom_ecmf = np.nan*np.zeros(all_ecmf.shape)
anom_bam = np.nan*np.zeros(all_bam.shape)
anom_chirps_bam = np.nan*np.zeros(all_chirps_bam.shape)
anom_chirps = np.nan*np.zeros(all_chirps.shape)

anom_ukmo_ens = np.nan*np.zeros(all_ukmo_ens.shape)
anom_ncep_ens = np.nan*np.zeros(all_ncep_ens.shape)
anom_ecmf_ens = np.nan*np.zeros(all_ecmf_ens.shape)
anom_bam_ens = np.nan*np.zeros(all_bam_ens.shape)


for y in np.arange(0,len(years)):
  print years[y]

  # compute independent climatologies for ACC and RMSE
  week_chirps_bam_mean = []
  week_chirps_bam_mean = np.nanmean(np.delete(all_chirps_bam,y,axis=0),axis=0)
  week_bam_mean = []
  week_bam_mean = np.nanmean(np.delete(all_bam_ens,y,axis=0),axis=(0,4))
  all_bam_ens_clim = []
  all_bam_ens_clim = np.nanmean(np.delete(all_bam_ens,y,axis=0),axis=0)

  week_chirps_mean = []
  week_chirps_mean = np.nanmean(np.delete(all_chirps,y,axis=0),axis=0)

  # climatology of ensemble mean with year of interest removed
  week_ecmf_mean = []
  week_ecmf_mean = np.nanmean(np.delete(all_ecmf_ens,y,axis=0),axis=(0,3,4))
  week_ncep_mean = []
  week_ncep_mean = np.nanmean(np.delete(all_ncep_ens,y,axis=0),axis=(0,3,4))
  week_ukmo_mean = []
  week_ukmo_mean = np.nanmean(np.delete(all_ukmo_ens,y,axis=0),axis=(0,3))

  # climatology of individual ensemble members with year of interest removed
  all_ukmo_ens_clim = []
  all_ukmo_ens_clim = np.nanmean(np.delete(all_ukmo_ens,y,axis=0),axis=0)
  all_ncep_ens_clim = []
  all_ncep_ens_clim = np.nanmean(np.delete(all_ncep_ens,y,axis=0),axis=0)
  all_ecmf_ens_clim = []
  all_ecmf_ens_clim = np.nanmean(np.delete(all_ecmf_ens,y,axis=0),axis=0)

  curr_n = np.zeros((nleads))# count the number of samples for each lead time in each year
  for l in np.arange(0,nleads):
    if season == 'NDJFM':
      for wbam in np.arange(0,5):
        for n in [0,1]:
          #if week_mid_dates_bam[y,l,wbam,n,1] in keep_mon:
          nsamp_week_only_bam[l] = nsamp_week_only_bam[l]+1	# count total number of weeks in season at each lead

          curr_chirps_anom = []
          curr_chirps_anom = all_chirps_bam[y,l,wbam,n,:,:] - week_chirps_bam_mean[l,wbam,n,:,:]
          anom_chirps_bam[y,l,wbam,n,:,:] = np.copy(curr_chirps_anom)

          curr_bam_anom = []
          curr_bam_anom = all_bam[y,l,wbam,n,:,:] - week_bam_mean[l,wbam,n,:,:]
          anom_bam[y,l,wbam,n,:,:] = np.copy(curr_bam_anom)

          nnan1 = np.sum(np.isnan(curr_chirps_anom))
          nnan2 = np.sum(np.isnan(curr_bam_anom))
          if (nnan1 < (len(all_lon)*len(all_lat))) & (nnan2 < (len(all_lon)*len(all_lat))):
            acc_pre[0,0,l,:,:] = acc_pre[0,0,l,:,:] + (curr_bam_anom*curr_chirps_anom)
            acc_pre[0,1,l,:,:] = acc_pre[0,1,l,:,:] + curr_bam_anom**2
            acc_pre[0,2,l,:,:] = acc_pre[0,2,l,:,:] + curr_chirps_anom**2
            rmse_mean[0,l,:,:] = rmse_mean[0,l,:,:] + (all_bam[y,l,wbam,n,:,:]-all_chirps_bam[y,l,wbam,n,:,:])**2

            # Compute current anomalies for chirps
            # overall bias throughout whole season, weekly bias and overall rmse
            for e in np.arange(0,bam_nmembers):
              nsamp_bam[l] = nsamp_bam[l] + 1# counts the number of samples at each week lead time
              nsamp_week_bam[l,wbam,n] = nsamp_week_bam[l,wbam,n] + 1
              bias[0,l,:,:] = bias[0,l,:,:] + (all_bam_ens[y,l,wbam,n,e,:,:] - all_chirps_bam[y,l,wbam,n,:,:])
              rmse[0,l,:,:] = rmse[0,l,:,:] + (all_bam_ens[y,l,wbam,n,e,:,:] - all_chirps_bam[y,l,wbam,n,:,:])**2
              curr_bam_anom = []
              curr_bam_anom = all_bam_ens[y,l,wbam,n,e,:,:] - all_bam_ens_clim[l,wbam,n,e,:,:]
              rmse_anom[0,l,:,:] = rmse[0,l,:,:] + (curr_bam_anom - curr_chirps_anom)**2
	      anom_bam_ens[y,l,wbam,n,e,:,:] = np.copy(curr_bam_anom)

          #else:
          #  nsamp_week_bam[l,wbam,n] = np.nan
    
    else:
      acc_pre[0,:,:,:,:] = np.nan
      bias[0,:,:,:] = np.nan
      rmse[0,:,:,:] = np.nan
      rmse_anom[0,:,:,:] = np.nan

    for w in np.arange(0,nweeks):# loop through starts
      # Only use forecasts in DJF since data is masked with NaN's for different months:
      if week_mid_dates[y,l,w,1] in [11,12,1,2,3,4,5]:
        #print 'lead time week '+str(l+1)+' mid-week date: '+str(int(week_mid_dates[y,l,w,0])).zfill(2)+'/'+str(int(week_mid_dates[y,l,w,1])).zfill(2)+'/'+str(int(week_mid_dates[y,l,w,2]))
        nsamp_week_only[l] = nsamp_week_only[l] + 1	# count total number of weeks in season at each lead
        curr_n[l] = curr_n[l] + 1
        # Compute current anomalies for chirps
        curr_chirps_anom = []
        curr_chirps_anom = all_chirps[y,l,w,:,:] - week_chirps_mean[l,w,:,:]
        anom_chirps[y,l,w,:,:] = np.copy(curr_chirps_anom)

        # overall bias throughout whole season, weekly bias and overall rmse
        for e in np.arange(0,ecmf_nmembers):
          if e < ukmo_nmembers:
            nsamp[3,l] = nsamp[3,l] + 1	# counts the number of samples at each week lead time
            nsamp_week[3,l,w] = nsamp_week[3,l,w] + 1
            bias[3,l,:,:] = bias[3,l,:,:] + (all_ukmo_ens[y,l,w,e,:,:] - all_chirps[y,l,w,:,:])
            week_bias[3,l,w,:,:] = week_bias[3,l,w,:,:] + (all_ukmo_ens[y,l,w,e,:,:] - all_chirps[y,l,w,:,:])
            rmse[3,l,:,:] = rmse[3,l,:,:] + (all_ukmo_ens[y,l,w,e,:,:] - all_chirps[y,l,w,:,:])**2
            curr_ukmo_anom = []
            curr_ukmo_anom = all_ukmo_ens[y,l,w,e,:,:] - all_ukmo_ens_clim[l,w,e,:,:]
            rmse_anom[3,l,:,:] = rmse_anom[3,l,:,:] + (curr_ukmo_anom - curr_chirps_anom)**2
	    anom_ukmo_ens[y,l,w,e,:,:] = np.copy(curr_ukmo_anom)


          if e < ncep_nmembers:
            for g in np.arange(0,ncep_lags):
              nsamp[2,l] = nsamp[2,l] + 1
              nsamp_week[2,l,w] = nsamp_week[2,l,w] + 1
              bias[2,l,:,:] = bias[2,l,:,:] + (all_ncep_ens[y,l,w,g,e,:,:] - all_chirps[y,l,w,:,:])
              week_bias[2,l,w,:,:] = week_bias[2,l,w,:,:] + (all_ncep_ens[y,l,w,g,e,:,:] - all_chirps[y,l,w,:,:])
              rmse[2,l,:,:] = rmse[2,l,:,:] + (all_ncep_ens[y,l,w,g,e,:,:] - all_chirps[y,l,w,:,:])**2
              curr_ncep_anom = []
              curr_ncep_anom = all_ncep_ens[y,l,w,g,e,:,:] - all_ncep_ens_clim[l,w,g,e,:,:]
              rmse_anom[2,l,:,:] = rmse_anom[2,l,:,:] + (curr_ncep_anom - curr_chirps_anom)**2
	      anom_ncep_ens[y,l,w,g,e,:,:] = np.copy(curr_ncep_anom)

          for g in np.arange(0,ecmf_lags):
            nsamp[1,l] = nsamp[1,l]+1
            nsamp_week[1,l,w] = nsamp_week[1,l,w]+1
            bias[1,l,:,:] = bias[1,l,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_chirps[y,l,w,:,:])
            week_bias[1,l,w,:,:] = week_bias[1,l,w,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_chirps[y,l,w,:,:])
            rmse[1,l,:,:] = rmse[1,l,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_chirps[y,l,w,:,:])**2
            curr_ecmf_anom = []
            curr_ecmf_anom = all_ecmf_ens[y,l,w,g,e,:,:] - all_ecmf_ens_clim[l,w,g,e,:,:]
            rmse_anom[1,l,:,:] = rmse_anom[1,l,:,:] + (curr_ecmf_anom - curr_chirps_anom)**2
	    anom_ecmf_ens[y,l,w,g,e,:,:] = np.copy(curr_ecmf_anom)

      else: # if not DJF
        week_bias[0:4,l,w,:,:] = np.nan
        nsamp_week[:,l,w] = np.nan


chirps_ano = np.swapaxes(anom_chirps,1,2).reshape((len(years),len(keep_mon),len(mweeki),nleads,len(all_lat),len(all_lon))).reshape((len(years)*len(keep_mon)*len(mweeki),nleads,len(all_lat),len(all_lon))) 
chirps_all =  np.swapaxes(all_chirps,1,2).reshape((len(years),len(keep_mon),len(mweeki),nleads,len(all_lat),len(all_lon))).reshape((len(years)*len(keep_mon)*len(mweeki),nleads,len(all_lat),len(all_lon)))
mjo_chirps = [0]*10; mjoa_chirps = [0]*10
idx = np.where((ma<1))
mjo_chirps[0] = chirps_all[idx,:,:,:].squeeze()
mjoa_chirps[0] = chirps_ano[idx,:,:,:].squeeze()
mjo_chirps[9] = chirps_all.copy()
mjoa_chirps[9] = chirps_ano.copy()
for i in range(1,9):
    idx = np.where((ma>=1) & (mp==(i)))
    mjo_chirps[i] = chirps_all[idx,:,:,:].squeeze()
    mjoa_chirps[i] = chirps_ano[idx,:,:,:].squeeze()

np.savez('mjo_chirps', mjo_chirps=mjo_chirps, mjoa_chirps=mjoa_chirps)

#bchrips
bchirps_ano = np.swapaxes(np.swapaxes(anom_chirps_bam,1,2),2,3).reshape((len(years)*len(keep_mon)*len(bweeki),nleads,len(all_lat),len(all_lon))) 
bchirps_all = np.swapaxes(np.swapaxes(all_chirps_bam,1,2),2,3).reshape((len(years)*len(keep_mon)*len(bweeki),nleads,len(all_lat),len(all_lon))) 
mjo_bchirps = [0]*10; mjoa_bchirps = [0]*10
idx = np.where((ba<1))
mjo_bchirps[0] = bchirps_all[idx,:,:,:].squeeze()
mjoa_bchirps[0] = bchirps_ano[idx,:,:,:].squeeze()
mjo_bchirps[9] = bchirps_all.copy()
mjoa_bchirps[9] = bchirps_ano.copy()
for i in range(1,9):
    idx = np.where((ba>=1) & (bp==(i)))
    mjo_bchirps[i] = bchirps_all[idx,:,:,:].squeeze()
    mjoa_bchirps[i] = bchirps_ano[idx,:,:,:].squeeze()

np.savez('mjo_bchirps', mjo_bchirps=mjo_bchirps, mjoa_bchirps=mjoa_bchirps)

print 'chirps'

ecmf_ano = np.swapaxes((anom_ecmf_ens.reshape((len(years),nleads,nweeks,ecmf_memlag,len(all_lat),len(all_lon)))),1,2).reshape((len(years),len(keep_mon),len(mweeki),nleads,ecmf_memlag,len(all_lat),len(all_lon))).reshape((len(years)*len(keep_mon)*len(mweeki),nleads,ecmf_memlag,len(all_lat),len(all_lon))) 
ecmf_all = np.swapaxes((all_ecmf_ens.reshape((len(years),nleads,nweeks,ecmf_memlag,len(all_lat),len(all_lon)))),1,2).reshape((len(years),len(keep_mon),len(mweeki),nleads,ecmf_memlag,len(all_lat),len(all_lon))).reshape((len(years)*len(keep_mon)*len(mweeki),nleads,ecmf_memlag,len(all_lat),len(all_lon)))
mjo_ecmf = [0]*10; mjoa_ecmf = [0]*10
idx = np.where((ma<1))
mjo_ecmf[0] = ecmf_all[idx,:,:,:,:].squeeze()
mjoa_ecmf[0] = ecmf_ano[idx,:,:,:,:].squeeze()
mjo_ecmf[9] = ecmf_all.copy()
mjoa_ecmf[9] = ecmf_ano.copy()
for i in range(1,9):
    idx = np.where((ma>=1) & (mp==(i)))
    mjo_ecmf[i] = ecmf_all[idx,:,:,:,:].squeeze()
    mjoa_ecmf[i] = ecmf_ano[idx,:,:,:,:].squeeze()

np.savez('mjo_ecmf', mjo_ecmf=mjo_ecmf, mjoa_ecmf=mjoa_ecmf)
print 'ecmf'

ncep_ano = np.swapaxes((anom_ncep_ens.reshape((len(years),nleads,nweeks,ncep_memlag,len(all_lat),len(all_lon)))),1,2).reshape((len(years),len(keep_mon),len(mweeki),nleads,ncep_memlag,len(all_lat),len(all_lon))).reshape((len(years)*len(keep_mon)*len(mweeki),nleads,ncep_memlag,len(all_lat),len(all_lon))) 
ncep_all = np.swapaxes((all_ncep_ens.reshape((len(years),nleads,nweeks,ncep_memlag,len(all_lat),len(all_lon)))),1,2).reshape((len(years),len(keep_mon),len(mweeki),nleads,ncep_memlag,len(all_lat),len(all_lon))).reshape((len(years)*len(keep_mon)*len(mweeki),nleads,ncep_memlag,len(all_lat),len(all_lon)))
mjo_ncep = [0]*10; mjoa_ncep = [0]*10
idx = np.where((ma<1))
mjo_ncep[0] = ncep_all[idx,:,:,:,:].squeeze()
mjoa_ncep[0] = ncep_ano[idx,:,:,:,:].squeeze()
mjo_ncep[9] = ncep_all.copy()
mjoa_ncep[9] = ncep_ano.copy()
for i in range(1,9):
    idx = np.where((ma>=1) & (mp==(i)))
    mjo_ncep[i] = ncep_all[idx,:,:,:,:].squeeze()
    mjoa_ncep[i] = ncep_ano[idx,:,:,:,:].squeeze()

np.savez('mjo_ncep', mjo_ncep=mjo_ncep, mjoa_ncep=mjoa_ncep)
print 'ncep'

ukmo_ano = np.swapaxes((anom_ukmo_ens),1,2).reshape((len(years),len(keep_mon),len(mweeki),nleads,ukmo_memlag,len(all_lat),len(all_lon))).reshape((len(years)*len(keep_mon)*len(mweeki),nleads,ukmo_memlag,len(all_lat),len(all_lon))) 
ukmo_all = np.swapaxes((all_ukmo_ens),1,2).reshape((len(years),len(keep_mon),len(mweeki),nleads,ukmo_memlag,len(all_lat),len(all_lon))).reshape((len(years)*len(keep_mon)*len(mweeki),nleads,ukmo_memlag,len(all_lat),len(all_lon)))
mjo_ukmo = [0]*10; mjoa_ukmo = [0]*10
idx = np.where((ma<1))
mjo_ukmo[0] = ukmo_all[idx,:,:,:,:].squeeze()
mjoa_ukmo[0] = ukmo_ano[idx,:,:,:,:].squeeze()
mjo_ukmo[9] = ukmo_all.copy()
mjoa_ukmo[9] = ukmo_ano.copy()
for i in range(1,9):
    idx = np.where((ma>=1) & (mp==(i)))
    mjo_ukmo[i] = ukmo_all[idx,:,:,:,:].squeeze()
    mjoa_ukmo[i] = ukmo_ano[idx,:,:,:,:].squeeze()

np.savez('mjo_ukmo', mjo_ukmo=mjo_ukmo, mjoa_ukmo=mjoa_ukmo)
print 'ukmo'

#bam
bam_ano = np.swapaxes(np.swapaxes(anom_bam_ens,1,2),2,3).reshape((len(years)*len(keep_mon)*len(bweeki),nleads,bam_memlag,len(all_lat),len(all_lon))) 
bam_all = np.swapaxes(np.swapaxes(all_bam_ens,1,2),2,3).reshape((len(years)*len(keep_mon)*len(bweeki),nleads,bam_memlag,len(all_lat),len(all_lon))) 
mjo_bam = [0]*10; mjoa_bam = [0]*10
idx = np.where((ba<1))
mjo_bam[0] = bam_all[idx,:,:,:,:].squeeze()
mjoa_bam[0] = bam_ano[idx,:,:,:,:].squeeze()
mjo_bam[9] = bam_all.copy()
mjoa_bam[9] = bam_ano.copy()
for i in range(1,9):
    idx = np.where((ba>=1) & (bp==(i)))
    mjo_bam[i] = bam_all[idx,:,:,:,:].squeeze()
    mjoa_bam[i] = bam_ano[idx,:,:,:,:].squeeze()

np.savez('mjo_bam', mjo_bam=mjo_bam, mjoa_bam=mjoa_bam)
print 'bam'


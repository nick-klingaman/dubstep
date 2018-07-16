'''
Analyse skill of weekly S2S precipitation forecasts
from UKMO, NCEP and ECMWF over a given region

1) Reads in data saved over the given region from
'save_weekly_forecasts_ukmo_ncep.py'
'save_weekly_forecasts_ecmf.py'

2) Computes mean precipitation, bias, anomaly correlation coefficient
and brier skill scores

3) Plots results as spatial maps

M. Young 29/06/2018
'''
from __future__ import division
import glob
import numpy as np
from netCDF4 import Dataset
import time as tt
from datetime import datetime, timedelta
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap
import calendar
execfile('date_str.py')
execfile('grab_data.py')

dir_in ='/group_workspaces/jasmin2/ncas_climate/users/myoung02/datasets/DUBSTEP/' # forecast/obs dir
dir_out = '/home/users/myoung02/DUBSTEP/' # outputdir to save figs

# number of ensemble members for each forecast
ukmo_nmembers = 7
ncep_nmembers = 4
ecmf_nmembers = 11
# number of lagged ensembles
ncep_lags = 7
ecmf_lags = 3

# Define region for analysis over Brazil
region = 'Brazil'
latlim = [-40,20]
lonlim = [-90,-20]
years = np.arange(2000,2010+1,1) # december will always correspond to year-1
nleads = 5 # number of lead times (in weeks) in the data

starts_ukmo = ['1025','1101','1109','1117','1125','1201',
               '1209','1217','1225','0101','0109','0117',
               '0125','0201','0209','0217','0225']


'''
Function: 'mask_percentiles':
Function takes in a 2d spatial forecast array and 2d spatial array of corresponding percentiles
and masks the forecast array according to input arguments, 'arg':
 'between': mask the forecast between 2 percentiles
'below': mask the forecast below a percentile
'above': mask the forecast above a percentile

'''
def mask_percentiles(data,percentile,arg):
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

'''
1. Read in all the data for all years
'''
print 'reading years...'
# Array which stores the start-date of each week in the data dims[year,leadtime]
week_mid_dates = np.zeros((len(years),nleads,len(starts_ukmo),3))
yc=0
for y in years:
  print y
  yc = yc + 1
  # create a mask to mask out all dates not within DJF
  for i in np.arange(0,len(starts_ukmo)):
    curr_start = starts_ukmo[i] # get current start date
    curr_mo = int(curr_start[0:2]) # month from start date
    curr_dd = int(curr_start[2:4]) # day from start date
    # use previous year if month between July and December
    if curr_mo in [7,8,9,10,11,12]:
      y_ec = y-1
    else:
      y_ec = np.copy(y)

    # Get start date as date object
    curr_ec_datetime =datetime.strptime(str(y_ec)+curr_start,'%Y%m%d')

    for l in np.arange(0,nleads):
      # get dates of the 4th day of each week (this keeps samples consistent at each lead time)
      # Note:
      #  1) if you use the date of the first day in the week this excludes forecasts
      #     with start-dates e.g. 30th November, causing different samples at each lead time.
      #  2) if you use the date of the last day in the week this excludes forecasts where
      #     the end of the week stretches over into the next month.
      week_mid_dates[yc-1,l,i,0] = (curr_ec_datetime + timedelta(days=((l)*7)+3)).day
      week_mid_dates[yc-1,l,i,1] = (curr_ec_datetime + timedelta(days=((l)*7)+3)).month
      week_mid_dates[yc-1,l,i,2] = (curr_ec_datetime + timedelta(days=((l)*7)+3)).year
      # curr_mid = str(week_mid_dates[yc-1,l,i,0])+'/'+str(week_mid_dates[yc-1,l,i,1])+'/'+str(week_mid_dates[yc-1,l,i,2])
      # print 'lead time week '+str(l+1)+' '+curr_mid

  # UKMO Data, needs converted to mm/d, since units are 1 kg/m2/s = 86400 mm/day.
  # UKMO has 7 ensemble members
  ukmo_file = dir_in+region+'_UKMO_ensemble_weekly_hindcasts_lead1-5_DJF_'+str(y)+'.nc'
  nc_fid = Dataset(ukmo_file, 'r')
  if y == years[0]:
    ukmo_lat = np.array(nc_fid.variables['latitude'][:])
    ukmo_lon = np.array(nc_fid.variables['longitude'][:])
    all_ukmo_ens = np.zeros((len(years),nleads,len(starts_ukmo),ukmo_nmembers,len(ukmo_lat),len(ukmo_lon)))
  all_ukmo_ens[yc-1,:,:,:,:,:] = np.array(nc_fid.variables['week_precip'][:])*86400# convert kg/m2/s to mm/d
  nc_fid.close()

  # NCEP units don't need converted. NCEP has 7 lags and 4 ensemble members
  ncep_file = dir_in+region+'_NCEP_ensemble_weekly_hindcasts_lead1-5_DJF_'+str(y)+'.nc'
  nc_fid = Dataset(ncep_file, 'r')
  if y == years[0]:
    ncep_lat = np.array(nc_fid.variables['latitude'][:])
    ncep_lon = np.array(nc_fid.variables['longitude'][:])
    all_ncep_ens = np.zeros((len(years),nleads,len(starts_ukmo),ncep_lags,ncep_nmembers,len(ncep_lat),len(ncep_lon)))
  all_ncep_ens[yc-1,:,:,:,:,:] = np.array(nc_fid.variables['week_precip'][:])
  nc_fid.close()

  # ECMWF units don't need converted. There are 3 lags and 11 ensemble members
  ecmf_file = dir_in+region+'_ECMWF_ensemble_weekly_hindcasts_lead1-5_DJF_'+str(y)+'.nc'
  nc_fid = Dataset(ecmf_file, 'r')
  if y == years[0]:
    ecmf_lat = np.array(nc_fid.variables['latitude'][:])
    ecmf_lon = np.array(nc_fid.variables['longitude'][:])
    all_ecmf_ens = np.zeros((len(years),nleads,len(starts_ukmo),ecmf_lags,ecmf_nmembers,len(ncep_lat),len(ncep_lon)))
  all_ecmf_ens[yc-1,:,:,:,:,:] = np.array(nc_fid.variables['week_precip'][:])
  nc_fid.close()

  # GPCP at resolution of UKMO (1 degree)
  gpcp_file = dir_in+region+'_GPCP_weekly_DJF_UKMO_times_'+str(y)+'.nc'
  nc_fid = Dataset(gpcp_file, 'r')
  if y == years[0]:
    gpcp_lat = np.array(nc_fid.variables['latitude'][:])
    gpcp_lon = np.array(nc_fid.variables['longitude'][:])
    all_gpcp4ukmo = np.zeros((len(years),nleads,len(starts_ukmo),len(gpcp_lat),len(gpcp_lon)))
  all_gpcp4ukmo[yc-1,:,:,:,:] = np.array(nc_fid.variables['week_precip'][:])
  nc_fid.close()

  # GPCP at resolution of NCEP & ECMWF
  gpcp_file = dir_in+region+'_GPCP_weekly_DJF_NCEP_times_'+str(y)+'.nc'
  nc_fid = Dataset(gpcp_file, 'r')
  if y == years[0]:
    gpcp_ncep_lat = np.array(nc_fid.variables['latitude'][:])
    gpcp_ncep_lon = np.array(nc_fid.variables['longitude'][:])
    all_gpcp4ncep = np.zeros((len(years),nleads,len(starts_ukmo),len(gpcp_ncep_lat),len(gpcp_ncep_lon)))
  all_gpcp4ncep[yc-1,:,:,:,:] = np.array(nc_fid.variables['week_precip'][:])
  nc_fid.close()

'''
2) Mask data for DJF, compute lag/ensemble means, and create a dry mask
'''
# mask data based on forecasts only in December, January and February
date_mask = np.where((week_mid_dates[:,:,:,1] > 2) & (week_mid_dates[:,:,:,1] < 12))
all_ukmo_ens[date_mask[0],date_mask[1],date_mask[2],:,:,:] = np.nan
all_ncep_ens[date_mask[0],date_mask[1],date_mask[2],:,:,:,:] = np.nan
all_ecmf_ens[date_mask[0],date_mask[1],date_mask[2],:,:,:,:] = np.nan
all_gpcp4ukmo[date_mask[0],date_mask[1],date_mask[2],:,:] = np.nan
all_gpcp4ncep[date_mask[0],date_mask[1],date_mask[2],:,:] = np.nan

# Forecast ensemble and lagged means
all_ukmo = np.nanmean(all_ukmo_ens,axis=3)
all_ncep = np.nanmean(all_ncep_ens,axis=(3,4))
all_ecmf = np.nanmean(all_ecmf_ens,axis=(3,4))

# average forecasts over lags only (leaving ensembles)
lag_ncep = np.nanmean(all_ncep_ens,axis=3)
lag_ecmf = np.nanmean(all_ecmf_ens,axis=3)

# Forecast means over years and ensembles (& lags for ncep/ecmf)
week_ukmo_mean = np.nanmean(all_ukmo_ens,axis=(0,3))
week_ncep_mean = np.nanmean(all_ncep_ens,axis=(0,3,4))
week_ecmf_mean = np.nanmean(all_ecmf_ens,axis=(0,3,4))

#  GPCP mean for each week
week_gpcp4ukmo_mean = np.nanmean(all_gpcp4ukmo,axis=0)
week_gpcp4ncep_mean = np.nanmean(all_gpcp4ncep,axis=0)

# Create a dry mask (GPCP seasonal clim < 1 mmd)
# this is used to mask Brier Skill Score later
# as BSS results are dodgy over dry regions
dry_mask_ukmo = np.nanmean(np.copy(week_gpcp4ukmo_mean),axis=(0,1))
dry_mask_ncep = np.nanmean(np.copy(week_gpcp4ncep_mean),axis=(0,1))
dry_mask_ukmo[dry_mask_ukmo < 1] = np.nan
dry_mask_ukmo[dry_mask_ukmo >= 1] = 1
dry_mask_ncep[dry_mask_ncep < 1] = np.nan
dry_mask_ncep[dry_mask_ncep >= 1] = 1


'''
3) Compute the Anomaly correlation coefficient (ACC), bias and RMSE
'''
# set up empty arrays
ukmo_acc_pre = np.zeros((3,nleads,len(ukmo_lat),len(ukmo_lon)))
ncep_acc_pre = np.zeros((3,nleads,len(ncep_lat),len(ncep_lon)))
ecmf_acc_pre = np.zeros((3,nleads,len(ecmf_lat),len(ecmf_lon)))

ukmo_bias = np.zeros((nleads,len(ukmo_lat),len(ukmo_lon)))
ncep_bias = np.zeros((nleads,len(ncep_lat),len(ncep_lon)))
ecmf_bias = np.zeros((nleads,len(ecmf_lat),len(ecmf_lon)))
ukmo_week_bias = np.zeros((nleads,len(starts_ukmo),len(ukmo_lat),len(ukmo_lon)))
ncep_week_bias = np.zeros((nleads,len(starts_ukmo),len(ncep_lat),len(ncep_lon)))
ecmf_week_bias = np.zeros((nleads,len(starts_ukmo),len(ecmf_lat),len(ecmf_lon)))

ukmo_rmse = np.zeros((nleads,len(ukmo_lat),len(ukmo_lon)))
ncep_rmse = np.zeros((nleads,len(ncep_lat),len(ncep_lon)))
ecmf_rmse = np.zeros((nleads,len(ecmf_lat),len(ecmf_lon)))

nsamp = np.zeros((3,nleads)) # count number of samples at each lead time
nsamp_week = np.zeros((3,nleads,len(starts_ukmo)))
nsamp_week_only = np.zeros((nleads))

for y in np.arange(0,len(years)):
  curr_n = np.zeros((nleads)) # count the number of samples for each lead time in each year
  for w in np.arange(0,len(starts_ukmo)): # loop through starts
    for l in np.arange(0,nleads):# loop through lead times
      # Only use forecasts in DJF since data is masked with NaN's for different months:
      if week_mid_dates[y,l,w,1] in [12,1,2]:
        print 'lead time week '+str(l+1)+' mid-week date: '+str(week_mid_dates[y,l,w,0])+'/'+str(week_mid_dates[y,l,w,1])+'/'+str(week_mid_dates[y,l,w,2])
        nsamp_week_only[l] = nsamp_week_only[l] + 1 # count total number of weeks in season at each lead
        curr_n[l] = curr_n[l] + 1

        # Compute current anomalies for GPCP
        curr_gpcp4ukmo_anom = all_gpcp4ukmo[y,l,w,:,:] - week_gpcp4ukmo_mean[l,w,:,:]
        curr_gpcp4ncep_anom = all_gpcp4ncep[y,l,w,:,:] - week_gpcp4ncep_mean[l,w,:,:]

        # overall bias throughout whole season, weekly bias and overall rmse
        for e in np.arange(0,ecmf_nmembers):
          if e < ukmo_nmembers:
            nsamp[0,l] = nsamp[0,l] + 1 # counts the number of samples at each week lead time
            nsamp_week[0,l,w] = nsamp_week[0,l,w] + 1
            ukmo_bias[l,:,:] = ukmo_bias[l,:,:] + (all_ukmo_ens[y,l,w,e,:,:] - all_gpcp4ukmo[y,l,w,:,:])
            ukmo_week_bias[l,w,:,:] = ukmo_week_bias[l,w,:,:] + (all_ukmo_ens[y,l,w,e,:,:] - all_gpcp4ukmo[y,l,w,:,:])
            ukmo_rmse[l,:,:] = ukmo_rmse[l,:,:] + (all_ukmo_ens[y,l,w,e,:,:] - all_gpcp4ukmo[y,l,w,:,:])**2

          if e < ncep_nmembers: # ncep
            for g in np.arange(0,ncep_lags): # loop through NCEP lags (7 days)
              nsamp[1,l] = nsamp[1,l] + 1
              nsamp_week[1,l,w] = nsamp_week[1,l,w] + 1
              ncep_bias[l,:,:] = ncep_bias[l,:,:] + (all_ncep_ens[y,l,w,g,e,:,:] - all_gpcp4ncep[y,l,w,:,:])
              ncep_week_bias[l,w,:,:] = ncep_week_bias[l,w,:,:] + (all_ncep_ens[y,l,w,g,e,:,:] - all_gpcp4ncep[y,l,w,:,:])
              ncep_rmse[l,:,:] = ncep_rmse[l,:,:] + (all_ncep_ens[y,l,w,g,e,:,:] - all_gpcp4ncep[y,l,w,:,:])**2

          # ecmf
          for g in np.arange(0,ecmf_lags): # loop through lags
            nsamp[2,l] = nsamp[2,l] + 1
            nsamp_week[2,l,w] = nsamp_week[2,l,w] + 1
            ecmf_bias[l,:,:] = ecmf_bias[l,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_gpcp4ncep[y,l,w,:,:])
            ecmf_week_bias[l,w,:,:] = ecmf_week_bias[l,w,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_gpcp4ncep[y,l,w,:,:])
            ecmf_rmse[l,:,:] = ecmf_rmse[l,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_gpcp4ncep[y,l,w,:,:])**2

        # for anomaly corelation coefficent just compare lagged, ensemble means
        curr_ukmo_anom = []
        curr_ukmo_anom = all_ukmo[y,l,w,:,:] - week_ukmo_mean[l,w,:,:]
        ukmo_acc_pre[0,l,:,:] = ukmo_acc_pre[0,l,:,:] + (curr_ukmo_anom*curr_gpcp4ukmo_anom)
        ukmo_acc_pre[1,l,:,:] = ukmo_acc_pre[1,l,:,:] + curr_ukmo_anom**2
        ukmo_acc_pre[2,l,:,:] = ukmo_acc_pre[2,l,:,:] + curr_gpcp4ukmo_anom**2

        curr_ncep_anom = []
        curr_ncep_anom = all_ncep[y,l,w,:,:] - week_ncep_mean[l,w,:,:]
        ncep_acc_pre[0,l,:,:] = ncep_acc_pre[0,l,:,:] + (curr_ncep_anom*curr_gpcp4ncep_anom)
        ncep_acc_pre[1,l,:,:] = ncep_acc_pre[1,l,:,:] + curr_ncep_anom**2
        ncep_acc_pre[2,l,:,:] = ncep_acc_pre[2,l,:,:] + curr_gpcp4ncep_anom**2

        curr_ecmf_anom = []
        curr_ecmf_anom = all_ecmf[y,l,w,:,:] - week_ecmf_mean[l,w,:,:]
        ecmf_acc_pre[0,l,:,:] = ecmf_acc_pre[0,l,:,:] + (curr_ecmf_anom*curr_gpcp4ncep_anom)
        ecmf_acc_pre[1,l,:,:] = ecmf_acc_pre[1,l,:,:] + curr_ecmf_anom**2
        ecmf_acc_pre[2,l,:,:] = ecmf_acc_pre[2,l,:,:] + curr_gpcp4ncep_anom**2

      else: # if not DJF
        ukmo_week_bias[l,w,:,:] = np.nan
        ncep_week_bias[l,w,:,:] = np.nan
        ncep_week_bias[l,w,:,:] = np.nan
        nsamp_week[:,l,w] = np.nan
  # print nsamp_week_only
  # print curr_n
  # print calendar.isleap(years[y])

# Compute ACC
ukmo_acc = ukmo_acc_pre[0,:,:,:]/np.sqrt(ukmo_acc_pre[1,:,:,:]*ukmo_acc_pre[2,:,:,:])
ncep_acc = ncep_acc_pre[0,:,:,:]/np.sqrt(ncep_acc_pre[1,:,:,:]*ncep_acc_pre[2,:,:,:])
ecmf_acc = ecmf_acc_pre[0,:,:,:]/np.sqrt(ecmf_acc_pre[1,:,:,:]*ecmf_acc_pre[2,:,:,:])

for l in np.arange(0,nleads):
  ukmo_bias[l,:,:] = ukmo_bias[l,:,:]/nsamp[0,l]
  ncep_bias[l,:,:] = ncep_bias[l,:,:]/nsamp[1,l]
  ecmf_bias[l,:,:] = ecmf_bias[l,:,:]/nsamp[2,l]
  ukmo_rmse[l,:,:] = np.sqrt(ukmo_rmse[l,:,:]/nsamp[0,l])
  ncep_rmse[l,:,:] = np.sqrt(ncep_rmse[l,:,:]/nsamp[1,l])
  ecmf_rmse[l,:,:] = np.sqrt(ecmf_rmse[l,:,:]/nsamp[2,l])
  for w in np.arange(0,len(starts_ukmo)):
    ukmo_week_bias[l,w,:,:] = ukmo_week_bias[l,w,:,:]/nsamp_week[0,l,w]
    ncep_week_bias[l,w,:,:] = ncep_week_bias[l,w,:,:]/nsamp_week[0,l,w]
    ecmf_week_bias[l,w,:,:] = ecmf_week_bias[l,w,:,:]/nsamp_week[0,l,w]

# alternatively compute bias using lagged ensemble means
# (there is a difference to the above method but the difference is v. small)
ukmo_bias_s = np.nanmean(all_ukmo - all_gpcp4ukmo,axis=(0,2))
ncep_bias_s = np.nanmean(all_ncep - all_gpcp4ncep,axis=(0,2))
ecmf_bias_s = np.nanmean(all_ecmf - all_gpcp4ncep,axis=(0,2))

# Average ACC, RMSE and BIAS over a sub-region at each lead time
sublon = [-60,-35]
sublat = [-25,-5]
subregion_ukmo_lat_id = np.where((ukmo_lat >= sublat[0]) & (ukmo_lat <= sublat[1]))[0]
subregion_ukmo_lon_id = np.where((ukmo_lon >= sublon[0]) & (ukmo_lon <= sublon[1]))[0]
subregion_ncep_lat_id = np.where((ncep_lat >= sublat[0]) & (ncep_lat <= sublat[1]))[0]
subregion_ncep_lon_id = np.where((ncep_lon >= sublon[0]) & (ncep_lon <= sublon[1]))[0]

acc_region = np.zeros((3,nleads))
rmse_region = np.zeros((3,nleads))
bias_region = np.zeros((3,nleads))
acc_region[0,:] = np.nanmean(ukmo_acc[:,subregion_ukmo_lat_id,:][:,:,subregion_ukmo_lon_id],axis=(1,2))
acc_region[1,:] = np.nanmean(ncep_acc[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))
acc_region[2,:] = np.nanmean(ecmf_acc[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))

bias_region[0,:] = np.nanmean(ukmo_bias[:,subregion_ukmo_lat_id,:][:,:,subregion_ukmo_lon_id],axis=(1,2))
bias_region[1,:] = np.nanmean(ncep_bias[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))
bias_region[2,:] = np.nanmean(ecmf_bias[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))

rmse_region[0,:] = np.nanmean(ukmo_rmse[:,subregion_ukmo_lat_id,:][:,:,subregion_ukmo_lon_id],axis=(1,2))
rmse_region[1,:] = np.nanmean(ncep_rmse[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))
rmse_region[2,:] = np.nanmean(ecmf_rmse[:,subregion_ncep_lat_id,:][:,:,subregion_ncep_lon_id],axis=(1,2))



'''
4) Compute brier skill score for each tercile
'''
percentiles = [100./3.,200./3.]
# work out the climatological frequency from the observations, for each week
ukmo_bs_clim = np.zeros((3,nleads,len(ukmo_lat),len(ukmo_lon)))
ncep_bs_clim = np.zeros((3,nleads,len(ncep_lat),len(ncep_lon)))

# work out brier score for all forecasts
ukmo_bs = np.zeros((3,nleads,len(ukmo_lat),len(ukmo_lon)))
ukmo_bs_bc = np.zeros((3,nleads,2,len(ukmo_lat),len(ukmo_lon)))
ncep_bs = np.zeros((3,nleads,len(ncep_lat),len(ncep_lon)))
ncep_bs_bc = np.zeros((3,nleads,2,len(ncep_lat),len(ncep_lon)))
ecmf_bs = np.zeros((3,nleads,len(ecmf_lat),len(ecmf_lon)))
ecmf_bs_bc = np.zeros((3,nleads,2,len(ecmf_lat),len(ecmf_lon)))

lag_ncep_bs = np.zeros((3,nleads,len(ncep_lat),len(ncep_lon)))
lag_ncep_bs_bc = np.zeros((3,nleads,2,len(ncep_lat),len(ncep_lon)))
lag_ecmf_bs = np.zeros((3,nleads,len(ecmf_lat),len(ecmf_lon)))
lag_ecmf_bs_bc = np.zeros((3,nleads,2,len(ecmf_lat),len(ecmf_lon)))

nw = 0
for l in np.arange(0,nleads):
  print 'Computing brier scores during lead week '+str(l+1)
  for y in np.arange(0,len(years)):
    for w in np.arange(0,12):
      if l == 0:
        nw = nw + 1
      if week_mid_dates[y,l,w,1] in [12,1,2]:
        if y == 0: # only compute percentiles in first year
          gpcp4ukmo_percentiles=np.nanpercentile(all_gpcp4ukmo[:,l,w,:,:],percentiles,axis=0)
          gpcp4ncep_percentiles=np.nanpercentile(all_gpcp4ncep[:,l,w,:,:],percentiles,axis=0)

        curr_gpcp4ukmo = np.copy(all_gpcp4ukmo[y,l,w,:,:])
        curr_gpcp4ncep = np.copy(all_gpcp4ncep[y,l,w,:,:])

        # loop through each tercile
        for ter in [0,1,2]:
          # empty arrays for counting probabilities
          c_ukmo = np.zeros((len(ukmo_lat),len(ukmo_lon)))
          c_ukmo_bc = np.zeros((2,len(ukmo_lat),len(ukmo_lon)))
          c_ncep = np.zeros((len(ncep_lat),len(ncep_lon)))
          c_ncep_bc = np.zeros((2,len(ncep_lat),len(ncep_lon)))
          c_ecmf = np.zeros((len(ecmf_lat),len(ecmf_lon)))
          c_ecmf_bc = np.zeros((2,len(ecmf_lat),len(ecmf_lon)))
          c_lag_ncep = np.zeros((len(ncep_lat),len(ncep_lon)))
          c_lag_ncep_bc = np.zeros((2,len(ncep_lat),len(ncep_lon)))
          c_lag_ecmf = np.zeros((len(ecmf_lat),len(ecmf_lon)))
          c_lag_ecmf_bc = np.zeros((2,len(ecmf_lat),len(ecmf_lon)))


          if ter == 0: # below 33 percentile
            p_gpcp4ukmo = mask_percentiles(curr_gpcp4ukmo,gpcp4ukmo_percentiles[0,:,:],'below')
            p_gpcp4ncep = mask_percentiles(curr_gpcp4ncep,gpcp4ncep_percentiles[0,:,:],'below')

          elif ter == 1: # between 33 and 66 percentile
            p_gpcp4ukmo = mask_percentiles(curr_gpcp4ukmo,gpcp4ukmo_percentiles,'between')
            p_gpcp4ncep = mask_percentiles(curr_gpcp4ncep,gpcp4ncep_percentiles,'between')

          elif ter == 2: # above 66 percentile
            p_gpcp4ukmo = mask_percentiles(curr_gpcp4ukmo,gpcp4ukmo_percentiles[1,:,:],'above')
            p_gpcp4ncep = mask_percentiles(curr_gpcp4ncep,gpcp4ncep_percentiles[1,:,:],'above')

          ukmo_bs_clim[ter,l,:,:] = ukmo_bs_clim[ter,l,:,:] + (0.33 - p_gpcp4ukmo)**2
          ncep_bs_clim[ter,l,:,:] = ncep_bs_clim[ter,l,:,:] + (0.33 - p_gpcp4ncep)**2
          # compute brier score for each lead time
          # compute the probability that the ensemble was
          # above/below each tercile
          for e in np.arange(0,ecmf_nmembers):

            if e < 7: # 1) do UKMO
              curr_ukmo = np.copy(all_ukmo_ens[y,l,w,e,:,:].squeeze())
              # empty array to store bias corrections
              curr_ukmo_bc = np.zeros((2,len(ukmo_lat),len(ukmo_lon)))

              # DO BIAS CORRECTIONS (UKMO only) for brier score
              # 1. Bias correct for given week in season
              curr_ukmo_bc[0,:,:] = np.copy(all_ukmo_ens[y,l,w,e,:,:].squeeze()) - np.copy(ukmo_week_bias[l,w,:,:].squeeze())
              # 2. Bias correct using overall DJF model bias
              curr_ukmo_bc[1,:,:] = np.copy(all_ukmo_ens[y,l,w,e,:,:].squeeze()) - np.copy(ukmo_bias[l,:,:].squeeze())
              tmp_ukmo = []
              tmp_ukmo_bc = np.zeros(curr_ukmo_bc.shape)
              if ter == 0: # below 33 percentile
                tmp_ukmo = mask_percentiles(curr_ukmo,gpcp4ukmo_percentiles[0,:,:],'below')
                for r in [0,1]: # loop for the 2 different bias corrections
                  tmp_ukmo_bc[r,:,:] = mask_percentiles(curr_ukmo_bc[r,:,:],gpcp4ukmo_percentiles[0,:,:],'below')

              elif ter == 1: # between 33 and 66 percentile
                tmp_ukmo = mask_percentiles(curr_ukmo,gpcp4ukmo_percentiles,'between')
                for r in [0,1]:
                  tmp_ukmo_bc[r,:,:] = mask_percentiles(curr_ukmo_bc[r,:,:],gpcp4ukmo_percentiles,'between')

              elif ter == 2: # above 66 percentile
                tmp_ukmo = mask_percentiles(curr_ukmo,gpcp4ukmo_percentiles[1,:,:],'above')
                for r in [0,1]:
                  tmp_ukmo_bc[r,:,:] = mask_percentiles(curr_ukmo_bc[r,:,:],gpcp4ukmo_percentiles[1,:,:],'above')

              c_ukmo = c_ukmo + tmp_ukmo
              c_ukmo_bc = c_ukmo_bc + tmp_ukmo_bc

            if e < 4: # bias correct NCEP
              # brier skill score for lagged ensembles and ensembles
              for g in np.arange(0,7): # lag loop
                curr_ncep = np.copy(all_ncep_ens[y,l,w,g,e,:,:])
                curr_ncep_bc = np.zeros((2,len(ncep_lat),len(ncep_lon)))
                # just correct bias in lagged ensemble means
                curr_ncep_bc[0,:,:] = np.copy(all_ncep_ens[y,l,w,g,e,:,:].squeeze()) - np.copy(ncep_week_bias[l,w,:,:].squeeze())
                curr_ncep_bc[1,:,:] = np.copy(all_ncep_ens[y,l,w,g,e,:,:].squeeze()) - np.copy(ncep_bias[l,:,:].squeeze())
                tmp_ncep = []
                tmp_ncep_bc = np.zeros(curr_ncep_bc.shape)
                if ter == 0: # below 33 percentile
                  tmp_ncep = mask_percentiles(curr_ncep,gpcp4ncep_percentiles[0,:,:],'below')
                  for r in [0,1]: # loop for the 2 different bias corrections
                    tmp_ncep_bc[r,:,:] = mask_percentiles(curr_ncep_bc[r,:,:],gpcp4ncep_percentiles[0,:,:],'below')

                elif ter == 1: # between 33 and 66 percentile
                  tmp_ncep = mask_percentiles(curr_ncep,gpcp4ncep_percentiles,'between')
                  for r in [0,1]:
                    tmp_ncep_bc[r,:,:] = mask_percentiles(curr_ncep_bc[r,:,:],gpcp4ncep_percentiles,'between')

                elif ter == 2: # above 66 percentile
                  tmp_ncep = mask_percentiles(curr_ncep,gpcp4ncep_percentiles[1,:,:],'above')
                  for r in [0,1]:
                    tmp_ncep_bc[r,:,:] = mask_percentiles(curr_ncep_bc[r,:,:],gpcp4ncep_percentiles[1,:,:],'above')

                c_ncep = c_ncep + tmp_ncep
                c_ncep_bc = c_ncep_bc + tmp_ncep_bc

            # ECMWF
            for g in np.arange(0,3):
              curr_ecmf = np.copy(all_ecmf_ens[y,l,w,g,e,:,:])
              curr_ecmf_bc = np.zeros((2,len(ecmf_lat),len(ecmf_lon)))
              curr_ecmf_bc[0,:,:] = np.copy(all_ecmf_ens[y,l,w,g,e,:,:].squeeze()) - np.copy(ecmf_week_bias[l,w,:,:].squeeze())
              curr_ecmf_bc[1,:,:] = np.copy(all_ecmf_ens[y,l,w,g,e,:,:].squeeze()) - np.copy(ecmf_bias[l,:,:].squeeze())
              tmp_ecmf = []
              tmp_ecmf_bc = np.zeros(curr_ncep_bc.shape)
              if ter == 0: # below 33 percentile
                tmp_ecmf = mask_percentiles(curr_ecmf,gpcp4ncep_percentiles[0,:,:],'below')
                for r in [0,1]: # loop for the 2 different bias corrections
                  tmp_ecmf_bc[r,:,:] = mask_percentiles(curr_ecmf_bc[r,:,:],gpcp4ncep_percentiles[0,:,:],'below')

              elif ter == 1: # between 33 and 66 percentile
                tmp_ecmf = mask_percentiles(curr_ecmf,gpcp4ncep_percentiles,'between')
                for r in [0,1]:
                  tmp_ecmf_bc[r,:,:] = mask_percentiles(curr_ecmf_bc[r,:,:],gpcp4ncep_percentiles,'between')

              elif ter == 2: # above 66 percentile
                tmp_ecmf = mask_percentiles(curr_ecmf,gpcp4ncep_percentiles[1,:,:],'above')
                for r in [0,1]:
                  tmp_ecmf_bc[r,:,:] = mask_percentiles(curr_ecmf_bc[r,:,:],gpcp4ncep_percentiles[1,:,:],'above')

              c_ecmf = c_ecmf + tmp_ecmf
              c_ecmf_bc = c_ecmf_bc + tmp_ecmf_bc

          # Probability and brier score
          p_ukmo = c_ukmo/ukmo_nmembers # probability
          p_ncep = c_ncep/(ncep_nmembers*7) # probability
          p_ecmf = c_ecmf/(ecmf_nmembers*3) # probability
          ukmo_bs[ter,l,:,:] = ukmo_bs[ter,l,:,:] + (p_ukmo - p_gpcp4ukmo)**2 # brier score sum
          ncep_bs[ter,l,:,:] = ncep_bs[ter,l,:,:] + (p_ncep - p_gpcp4ncep)**2 # brier score sum
          ecmf_bs[ter,l,:,:] = ecmf_bs[ter,l,:,:] + (p_ecmf - p_gpcp4ncep)**2 # brier score sum

          # Bias corrected probability and brier score
          p_ukmo_bc = c_ukmo_bc/ukmo_nmembers
          p_ncep_bc = c_ncep_bc/(ncep_nmembers*7)
          p_ecmf_bc = c_ecmf_bc/(ecmf_nmembers*3)
          for r in [0,1]:
            ukmo_bs_bc[ter,l,r,:,:] = ukmo_bs_bc[ter,l,r,:,:] + (p_ukmo_bc[r,:,:] - p_gpcp4ukmo)**2
            ncep_bs_bc[ter,l,r,:,:] = ncep_bs_bc[ter,l,r,:,:] + (p_ncep_bc[r,:,:] - p_gpcp4ncep)**2
            ecmf_bs_bc[ter,l,r,:,:] = ecmf_bs_bc[ter,l,r,:,:] + (p_ecmf_bc[r,:,:] - p_gpcp4ncep)**2

ukmo_bs_clim = ukmo_bs_clim/nw
ncep_bs_clim = ncep_bs_clim/nw
ukmo_bss = (ukmo_bs_clim - (ukmo_bs/nw))/ukmo_bs_clim
ncep_bss = (ncep_bs_clim - (ncep_bs/nw))/ncep_bs_clim
ecmf_bss = (ncep_bs_clim - (ecmf_bs/nw))/ncep_bs_clim

ukmo_bss_bc = np.zeros(ukmo_bs_bc.shape)
ncep_bss_bc = np.zeros(ncep_bs_bc.shape)
ecmf_bss_bc = np.zeros(ecmf_bs_bc.shape)
for r in [0,1]:
  ukmo_bss_bc[:,:,r,:,:] = (ukmo_bs_clim - (ukmo_bs_bc[:,:,r,:,:]/nw))/ukmo_bs_clim
  ncep_bss_bc[:,:,r,:,:] = (ncep_bs_clim - (ncep_bs_bc[:,:,r,:,:]/nw))/ncep_bs_clim
  ecmf_bss_bc[:,:,r,:,:] = (ncep_bs_clim - (ecmf_bs_bc[:,:,r,:,:]/nw))/ncep_bs_clim


'''
5) Plot metrics
- regional averages at each lead time
- spatial maps of mean precipitation, model biases, ACC, BSS
'''
print 'Plotting figs...'

model_ls = ['UKMO','NCEP','ECMWF']
fmt_ls = ['-o','-x','-+']
ms_ls = [8,10,10]
mw_ls =[1,2,2]
col_ls = ['firebrick','dodgerblue','forestgreen']
leads = np.arange(1,nleads+1,1)

fname_plot = dir_out+region+'_S2S_weekly_ACC_DJF_region_average_'+str(years[0])+'_'+str(years[len(years)-1])
size1=[5,3.5]
fig = plt.figure(figsize=(size1[0],size1[1]))
for p in [0,1,2]:
  plt.plot(leads,acc_region[p,:],fmt_ls[p],color=col_ls[p],linewidth=2,ms=ms_ls[p],alpha=0.9,mew=mw_ls[p])
plt.ylim([0,0.62])
plt.xlim([0.8,5.2])
plt.ylabel('ACC')
plt.xlabel('Lead time (weeks)')
plt.legend(model_ls,loc=0)
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
plt.savefig(fname_plot+'.png',bbox_inches='tight')
#plt.show()
plt.close()

fname_plot = dir_out+region+'_S2S_weekly_RMSE_DJF_region_average_'+str(years[0])+'_'+str(years[len(years)-1])
size1=[5,3.5]
fig = plt.figure(figsize=(size1[0],size1[1]))
for p in [0,1,2]:
  plt.plot(leads,rmse_region[p,:],fmt_ls[p],color=col_ls[p],linewidth=2,ms=ms_ls[p],alpha=0.9,mew=mw_ls[p])
plt.ylim([4,8])
plt.xlim([0.8,5.2])
plt.ylabel('RMSE (mm d$^{-1}$)')
plt.xlabel('Lead time (weeks)')
plt.legend(model_ls,loc=0)
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
plt.savefig(fname_plot+'.png',bbox_inches='tight')
#plt.show()
plt.close()

fname_plot = dir_out+region+'_S2S_weekly_Bias_DJF_region_average_'+str(years[0])+'_'+str(years[len(years)-1])
size1=[5,3.5]
fig = plt.figure(figsize=(size1[0],size1[1]))
for p in [0,1,2]:
  plt.plot(leads,bias_region[p,:],fmt_ls[p],color=col_ls[p],linewidth=2,ms=ms_ls[p],alpha=0.9,mew=mw_ls[p])
  plt.plot(leads,np.zeros(len(leads)),'-k')
plt.ylim([-1,2])
plt.xlim([0.8,5.2])
plt.ylabel('Bias (mm d$^{-1}$)')
plt.xlabel('Lead time (weeks)')
plt.legend(model_ls,loc=0,prop={'size': 10})
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
plt.savefig(fname_plot+'.png',bbox_inches='tight')
#plt.show()
plt.close()


# Plot DJF mean precipitation
precip_colors = ["#fdfdfd","#f2f2f2","#bfbfbf","#04e9e7","#019ff4","#0300f4","#02fd02","#01c501","#008e00","#fdf802","#e5bc00","#fd9500","#fd0000", "#d40000","#bc0000","#f800fd","#9854c6"]
precip_colormap = matplotlib.colors.ListedColormap(precip_colors)
cols = precip_colormap
cols = 'PuBu'
cmin = 0
cmax = 16
cspc = 2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'Mean precipitation (mm d$^{-1}$)'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
lw = 1
gl = 20

size1 = [5,3.5]
fname_plot = dir_out+region+'_GPCC_weekly_mean_DJF_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
mymap = Basemap(projection='cyl',resolution='l',\
        llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
        llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
mymap.drawcoastlines(linewidth=lw)
mymap.drawcountries(linewidth=lw)
mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
uncal = mymap.pcolormesh(x,y,np.nanmean(week_gpcp4ukmo_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
plt.title('GPCP')
plt.colorbar(uncal,label=clabel1,extend='max')
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
#plt.show()
plt.close()

nrow = 4
ncol = nleads
size1=[15,12]
fname_plot = dir_out+region+'_S2S_weekly_mean_DJF_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
for n in np.arange(0,nleads):
  if n == 2:
    plt.subplot(nrow,ncol,3)
    mymap = Basemap(projection='cyl',resolution='l',\
            llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
            llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
    uncal = mymap.pcolormesh(x,y,np.nanmean(week_gpcp4ukmo_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title('GPCP')

  plt.subplot(nrow,ncol,n+6)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
  uncal = mymap.pcolormesh(x,y,np.nanmean(week_ukmo_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)

  plt.title('UKMO Week '+str(n+1))

  plt.subplot(nrow,ncol,n+11)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  uncal = mymap.pcolormesh(x,y,np.nanmean(week_ncep_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('NCEP Week '+str(n+1))

  plt.subplot(nrow,ncol,n+16)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  uncal = mymap.pcolormesh(x,y,np.nanmean(week_ncep_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('ECMWF Week '+str(n+1))

plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.6)
fig.subplots_adjust(right=0.90)
cbar_pos = [0.92, 0.10, 0.015, 0.35] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='max')
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
#plt.show()
plt.close()

# plot ACC just at lead times Week 1, Week 3 and Week 5 (save on panels)
cols = 'RdYlBu'
cmin = -1
cmax = 1
cspc = 0.2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'ACC'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)

lw = 1
nrow = 3
ncol = 3
size1=[10,9.5]
gl = 20
fname_plot = dir_out+region+'_S2S_weekly_sub_ACC_DJF_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
nc = 0
for n in [0,2,4]:
  nc = nc + 1
  plt.subplot(nrow,ncol,nc)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if nc in [1,4,7]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
  uncal = mymap.pcolormesh(x,y,ukmo_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)

  plt.title('UKMO Week '+str(n+1))

  nc = nc + 1
  plt.subplot(nrow,ncol,nc)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if nc in [1,4,7]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  uncal = mymap.pcolormesh(x,y,ncep_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('NCEP Week '+str(n+1))

  nc = nc + 1
  plt.subplot(nrow,ncol,nc)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if nc in [1,4,7]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  uncal = mymap.pcolormesh(x,y,ecmf_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('ECMWF Week '+str(n+1))

plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.6)
fig.subplots_adjust(right=0.90)
cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
#plt.savefig(fname_plot+'.png',bbox_inches='tight')
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
#plt.show()
plt.close()


# plot ACC
cols = 'RdYlBu'
cmin = -1
cmax = 1
cspc = 0.2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'ACC'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)

lw = 1
nrow = 3
ncol = nleads
size1=[15,9]
gl = 20
fname_plot = dir_out+region+'_S2S_weekly_ACC_DJF_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
for n in np.arange(0,nleads):
  plt.subplot(nrow,ncol,n+1)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
  uncal = mymap.pcolormesh(x,y,ukmo_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)

  plt.title('UKMO Week '+str(n+1))

  plt.subplot(nrow,ncol,n+6)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  uncal = mymap.pcolormesh(x,y,ncep_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('NCEP Week '+str(n+1))

  plt.subplot(nrow,ncol,n+11)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  uncal = mymap.pcolormesh(x,y,ecmf_acc[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('ECMWF Week '+str(n+1))

plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.6)
fig.subplots_adjust(right=0.90)
cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
#plt.savefig(fname_plot+'.png',bbox_inches='tight')
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
#plt.show()
plt.close()


# plot Bias
cols = 'RdBu'
cmin = -10
cmax = 10
cspc = 2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'Bias (mm d$^{-1}$)'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)

lw = 1
nrow = 3
ncol = nleads
size1=[15,9]
gl = 20
fname_plot = dir_out+region+'_S2S_weekly_Bias_DJF_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
for n in np.arange(0,nleads):
  plt.subplot(nrow,ncol,n+1)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
  tmp_bias = []
  tmp_bias = np.copy(ukmo_bias[n,:,:])
  tmp_bias[tmp_bias == 0] = np.nan
  uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('UKMO Week '+str(n+1))

  plt.subplot(nrow,ncol,n+6)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))

  tmp_bias = []
  tmp_bias = np.copy(ncep_bias[n,:,:])
  tmp_bias[tmp_bias == 0] = np.nan
  uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('NCEP Week '+str(n+1))

  plt.subplot(nrow,ncol,n+11)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  tmp_bias = []
  tmp_bias = np.copy(ecmf_bias[n,:,:])
  tmp_bias[tmp_bias == 0] = np.nan
  uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('ECMWF Week '+str(n+1))

plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.4)
fig.subplots_adjust(right=0.90)
cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
#plt.show()
plt.close()

# plot bias for lead times week 1 and 5 only
nrow = 2
ncol = 3
size1=[10,7]
week_want = [0,4] # index for weeks want
fname_plot = dir_out+region+'_S2S_weekly_1&5_Bias_DJF_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
nc = 0
for n in week_want:
  nc = nc + 1
  plt.subplot(nrow,ncol,nc)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if nc in [1,4]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
  uncal = mymap.pcolormesh(x,y,ukmo_bias[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('UKMO Week '+str(n+1))
  nc = nc + 1

  plt.subplot(nrow,ncol,nc)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)

  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if nc in [1,4]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  #uncal = mymap.contourf(x,y,ukmo_acc[n,:,:],clevs,cmap=cols,extend='both')
  uncal = mymap.pcolormesh(x,y,ncep_bias[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('NCEP Week '+str(n+1))

  nc = nc + 1
  plt.subplot(nrow,ncol,nc)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if nc in [1,4]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  uncal = mymap.pcolormesh(x,y,ecmf_bias[n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('ECMWF Week '+str(n+1))

plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.4)
fig.subplots_adjust(right=0.90)
cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
#plt.show()
plt.close()



# plot RMSE
cols = 'YlOrRd'
cmin = 0
cmax = 16
cspc = 2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'RMSE (mm d$^{-1}$)'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)

lw = 1
nrow = 3
ncol = nleads
size1=[15,9]
gl = 20
fname_plot = dir_out+region+'_S2S_weekly_RMSE_DJF_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
for n in np.arange(0,nleads):
  plt.subplot(nrow,ncol,n+1)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
  tmp_bias = []
  tmp_bias = np.copy(ukmo_rmse[n,:,:])
  tmp_bias[tmp_bias == 0] = np.nan
  uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('UKMO Week '+str(n+1))

  plt.subplot(nrow,ncol,n+6)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  tmp_bias = []
  tmp_bias = np.copy(ncep_rmse[n,:,:])
  tmp_bias[tmp_bias == 0] = np.nan
  uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('NCEP Week '+str(n+1))

  plt.subplot(nrow,ncol,n+11)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
          llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  tmp_bias = []
  tmp_bias = np.copy(ecmf_rmse[n,:,:])
  tmp_bias[tmp_bias == 0] = np.nan
  uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('ECMWF Week '+str(n+1))

plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.4)
fig.subplots_adjust(right=0.90)
cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='max')
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
#plt.show()
plt.close()

# plot BSS
tercile_name = ['Below','Normal','Above']
cols = 'RdYlBu'
clabel1 = 'BSS'
lw = 1
nrow = 3
ncol = 5
size1=[15,9]
gl = 20

cmin = -0.5
cmax = 0.5
cspc = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
norm = BoundaryNorm(boundaries=clevs, ncolors=256)

for tc in [0,1,2]:
  fname_plot = dir_out+region+'_S2S_weekly_BSS_DJF_'+str(years[0])+'_'+str(years[len(years)-1])+'_'+tercile_name[tc]
  fig = plt.figure(figsize=(size1[0],size1[1]))

  for n in np.arange(0,nleads):
    plt.subplot(nrow,ncol,n+1)
    mymap = Basemap(projection='cyl',resolution='l',\
            llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
            llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
    mask_array = np.ma.array(ukmo_bss[tc,n,:,:]*dry_mask_ukmo, mask=np.isnan(dry_mask_ukmo))
    uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title('UKMO Week '+str(n+1))

    plt.subplot(nrow,ncol,n+6)
    mymap = Basemap(projection='cyl',resolution='l',\
            llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
            llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
    mask_array = np.ma.array(ncep_bss[tc,n,:,:]*dry_mask_ncep, mask=np.isnan(dry_mask_ncep))
    uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title('NCEP Week '+str(n+1))

    plt.subplot(nrow,ncol,n+11)
    mymap = Basemap(projection='cyl',resolution='l',\
            llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
            llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
    mask_array = np.ma.array(ecmf_bss[tc,n,:,:]*dry_mask_ncep, mask=np.isnan(dry_mask_ncep))
    uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title('ECMWF Week '+str(n+1))

  plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.4)
  fig.subplots_adjust(right=0.90)
  cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]

  cbar_ax = fig.add_axes(cbar_pos)
  cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
  plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
  #plt.show()
  plt.close()

  fname_plot = dir_out+region+'_S2S_weekly_BSS_BiasCorrected_DJF_'+str(years[0])+'_'+str(years[len(years)-1])+'_'+tercile_name[tc]
  fig = plt.figure(figsize=(size1[0],size1[1]))
  for n in np.arange(0,nleads):
    plt.subplot(nrow,ncol,n+1)
    mymap = Basemap(projection='cyl',resolution='l',\
            llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
            llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
    mask_array = np.ma.array(ukmo_bss_bc[tc,n,bb,:,:]*dry_mask_ukmo, mask=np.isnan(dry_mask_ukmo))
    uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title('UKMO Week '+str(n+1))

    plt.subplot(nrow,ncol,n+6)
    mymap = Basemap(projection='cyl',resolution='l',\
            llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
            llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
    #uncal = mymap.contourf(x,y,ukmo_acc[n,:,:],clevs,cmap=cols,extend='both')
    mask_array = np.ma.array(ncep_bss_bc[tc,n,bb,:,:]*dry_mask_ncep, mask=np.isnan(dry_mask_ncep))
    uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title('NCEP Week '+str(n+1))

    plt.subplot(nrow,ncol,n+11)
    mymap = Basemap(projection='cyl',resolution='l',\
            llcrnrlat=latlim[0],urcrnrlat=latlim[1],\
            llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
    mask_array = np.ma.array(ecmf_bss_bc[tc,n,bb,:,:]*dry_mask_ncep, mask=np.isnan(dry_mask_ncep))
    uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title('ECMWF Week '+str(n+1))

  plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.4)
  fig.subplots_adjust(right=0.90)
  cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
  cbar_ax = fig.add_axes(cbar_pos)
  cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
  plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
  #plt.show()
  plt.close()

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
dir_out	= '/home/users/myoung02/DUBSTEP_paper_results/'
dir_out2 = dir_out+'individual_panels/'
dir_prob	= dir_in+'hindcast_probabilities/'

model_ls= ['BAM','ECMWF','NCEP','UKMO']
# number of ensemble members for each forecast
ukmo_nmembers = 7
ncep_nmembers = 4
ecmf_nmembers = 11
bam_nmembers = 11
# number of lagged ensembles
ncep_lags = 7
ecmf_lags = 3

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

season = season_ls[1]
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
all_ukmo_ens_clim = np.nanmean(all_ukmo_ens,axis=0)
all_ncep_ens_clim = np.nanmean(all_ncep_ens,axis=0)
all_ecmf_ens_clim = np.nanmean(all_ecmf_ens,axis=0)
all_bam_ens_clim = np.nanmean(all_bam_ens,axis=0)

# average forecasts over lags only (leaving ensembles)
lag_ncep = np.nanmean(all_ncep_ens,axis=3)# (11, 5, 20, 4, 40, 47)
lag_ecmf = np.nanmean(all_ecmf_ens,axis=3)# (11, 5, 20, 11, 40, 47)

# Forecast means over years and ensembles (& lags for ncep/ecmf)
week_ukmo_mean = np.nanmean(all_ukmo_ens,axis=(0,3))# (5, 20, 40, 47)
week_ncep_mean = np.nanmean(all_ncep_ens,axis=(0,3,4))
week_ecmf_mean = np.nanmean(all_ecmf_ens,axis=(0,3,4))
week_bam_mean = np.nanmean(all_bam_ens,axis=(0,4))
week_mean = [week_bam_mean,week_ecmf_mean,week_ncep_mean,week_ukmo_mean]

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

for y in np.arange(0,len(years)):
  curr_n = np.zeros((nleads))# count the number of samples for each lead time in each year
  for l in np.arange(0,nleads):
    if season == 'NDJFM':
      for wbam in np.arange(0,5):
        for n in [0,1]:
          if week_mid_dates_bam[y,l,wbam,n,1] in keep_mon:
            nsamp_week_only_bam[l] = nsamp_week_only_bam[l]+1	# count total number of weeks in season at each lead

            curr_chirps_anom = []
            curr_chirps_anom = all_chirps_bam[y,l,wbam,n,:,:] - week_chirps_bam_mean[l,wbam,n,:,:]
            curr_bam_anom = []
            curr_bam_anom = all_bam[y,l,wbam,n,:,:] - week_bam_mean[l,wbam,n,:,:]
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

          else:
            nsamp_week_bam[l,wbam,n] = np.nan
    else:
      acc_pre[0,:,:,:,:] = np.nan
      bias[0,:,:,:] = np.nan
      rmse[0,:,:,:] = np.nan
      rmse_anom[0,:,:,:] = np.nan

    for w in np.arange(0,nweeks):# loop through starts
      # Only use forecasts in DJF since data is masked with NaN's for different months:
      if week_mid_dates[y,l,w,1] in keep_mon:
        #print 'lead time week '+str(l+1)+' mid-week date: '+str(int(week_mid_dates[y,l,w,0])).zfill(2)+'/'+str(int(week_mid_dates[y,l,w,1])).zfill(2)+'/'+str(int(week_mid_dates[y,l,w,2]))
        nsamp_week_only[l] = nsamp_week_only[l] + 1	# count total number of weeks in season at each lead
        curr_n[l] = curr_n[l] + 1
        # Compute current anomalies for chirps
        curr_chirps_anom = all_chirps[y,l,w,:,:] - week_chirps_mean[l,w,:,:]

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

          for g in np.arange(0,ecmf_lags):
            nsamp[1,l] = nsamp[1,l]+1
            nsamp_week[1,l,w] = nsamp_week[1,l,w]+1
            bias[1,l,:,:] = bias[1,l,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_chirps[y,l,w,:,:])
            week_bias[1,l,w,:,:] = week_bias[1,l,w,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_chirps[y,l,w,:,:])
            rmse[1,l,:,:] = rmse[1,l,:,:] + (all_ecmf_ens[y,l,w,g,e,:,:] - all_chirps[y,l,w,:,:])**2
            curr_ecmf_anom = []
            curr_ecmf_anom = all_ecmf_ens[y,l,w,g,e,:,:] - all_ecmf_ens_clim[l,w,g,e,:,:]
            rmse_anom[1,l,:,:] = rmse_anom[1,l,:,:] + (curr_ecmf_anom - curr_chirps_anom)**2

        curr_ecmf_anom = []
        curr_ecmf_anom = all_ecmf[y,l,w,:,:] - week_ecmf_mean[l,w,:,:]
        acc_pre[1,0,l,:,:] = acc_pre[1,0,l,:,:] + (curr_ecmf_anom*curr_chirps_anom)
        acc_pre[1,1,l,:,:] = acc_pre[1,1,l,:,:] + curr_ecmf_anom**2
        acc_pre[1,2,l,:,:] = acc_pre[1,2,l,:,:] + curr_chirps_anom**2
        rmse_mean[1,l,:,:] = rmse_mean[1,l,:,:] + (all_ecmf[y,l,w,:,:]-all_chirps[y,l,w,:,:])**2

        curr_ncep_anom = []
        curr_ncep_anom = all_ncep[y,l,w,:,:] - week_ncep_mean[l,w,:,:]
        acc_pre[2,0,l,:,:] = acc_pre[2,0,l,:,:] + (curr_ncep_anom*curr_chirps_anom)
        acc_pre[2,1,l,:,:] = acc_pre[2,1,l,:,:] + curr_ncep_anom**2
        acc_pre[2,2,l,:,:] = acc_pre[2,2,l,:,:] + curr_chirps_anom**2
        rmse_mean[2,l,:,:] = rmse_mean[2,l,:,:] + (all_ncep[y,l,w,:,:]-all_chirps[y,l,w,:,:])**2

        # for anomaly corelation coefficent just compare lagged, ensemble means
        curr_ukmo_anom = []
        curr_ukmo_anom = all_ukmo[y,l,w,:,:] - week_ukmo_mean[l,w,:,:]
        acc_pre[3,0,l,:,:] = acc_pre[3,0,l,:,:] + (curr_ukmo_anom*curr_chirps_anom)
        acc_pre[3,1,l,:,:] = acc_pre[3,1,l,:,:] + curr_ukmo_anom**2
        acc_pre[3,2,l,:,:] = acc_pre[3,2,l,:,:] + curr_chirps_anom**2
        rmse_mean[3,l,:,:] = rmse_mean[3,l,:,:] + (all_ukmo[y,l,w,:,:]-all_chirps[y,l,w,:,:])**2

      else: # if not DJF
        week_bias[0:4,l,w,:,:] = np.nan
        nsamp_week[:,l,w] = np.nan

# Compute ACC
for m in np.arange(0,4):
  acc[m,:,:,:] = acc_pre[m,0,:,:,:]/np.sqrt(acc_pre[m,1,:,:,:]*acc_pre[m,2,:,:,:])

for l in np.arange(0,nleads):
  bias[0,l,:,:]=bias[0,l,:,:]/nsamp_bam[l]
  rmse[0,l,:,:]=np.sqrt(rmse[0,l,:,:]/nsamp_bam[l])
  rmse_anom[0,l,:,:]=np.sqrt(rmse_anom[0,l,:,:]/nsamp_bam[l])

  rmse_mean[0,l,:,:]=np.sqrt(rmse_mean[0,l,:,:]/nsamp_bam[l])
  for m in [1,2,3]:
    bias[m,l,:,:]=bias[m,l,:,:]/nsamp[m,l]
    rmse[m,l,:,:]=np.sqrt(rmse[m,l,:,:]/nsamp[m,l])
    rmse_anom[m,l,:,:]=np.sqrt(rmse_anom[m,l,:,:]/nsamp[m,l])
    rmse_mean[m,l,:,:]=np.sqrt(rmse_mean[m,l,:,:]/nsamp[m,l])
    for w in np.arange(0,nweeks):
      week_bias[m,l,w,:,:]=week_bias[m,l,w,:,:]/nsamp_week[m,l,w]

# save ACC
fname_acc = dir_in+'ACC_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_acc,acc)

###
# checl effect
ensemble_bin = np.arange(0,7)
acc_ens = np.zeros((nleads,len(ensemble_bin),len(all_lat),len(all_lon)))
acc_pre_ens = np.zeros((3,nleads,len(ensemble_bin),len(all_lat),len(all_lon)))
for e in ensemble_bin:
  print e
  tmp_week_mean = np.nanmean(all_ncep_ens[:,:,:,0:e+1,:,:,:],axis=(0,3,4))
  tmp_all = np.nanmean(all_ncep_ens[:,:,:,0:e+1,:,:,:],axis=(3,4))
  for l in np.arange(0,nleads):
    for y in np.arange(0,len(years)):
      for w in np.arange(0,nweeks):# loop through starts
        # Only use forecasts in DJF since data is masked with NaN's for different months:
        if week_mid_dates[y,l,w,1] in keep_mon:
          # Compute current anomalies for chirps
          curr_chirps_anom = all_chirps[y,l,w,:,:] - week_chirps_mean[l,w,:,:]
          curr_anom  = []
          curr_anom = tmp_all[y,l,w,:,:] - tmp_week_mean[l,w,:,:]
          acc_pre_ens[0,l,e,:,:] = acc_pre_ens[0,l,e,:,:] + (curr_anom*curr_chirps_anom)
          acc_pre_ens[1,l,e,:,:] = acc_pre_ens[1,l,e,:,:] + curr_anom**2
          acc_pre_ens[2,l,e,:,:] = acc_pre_ens[2,l,e,:,:] + curr_chirps_anom**2
acc_ens[:,:,:,:] = acc_pre_ens[0,:,:,:,:]/np.sqrt(acc_pre_ens[1,:,:,:,:]*acc_pre_ens[2,:,:,:,:])


ensemble_bin_ecmf = np.arange(0,ecmf_nmembers*ecmf_lags)
acc_ens_ecmf = np.zeros((nleads,len(ensemble_bin_ecmf),len(all_lat),len(all_lon)))
acc_pre_ens_ecmf = np.zeros((3,nleads,len(ensemble_bin_ecmf),len(all_lat),len(all_lon)))

e = 0
em = 0
for j in np.arange(0,ecmf_lags):
  el =0
  em = em +1
  for i in np.arange(0,ecmf_nmembers):
    el = el+1
    e = e + 1
    print e
    tmp_week_mean = np.nanmean(all_ecmf_ens[:,:,:,0:el,0:em,:,:],axis=(0,3,4)).squeeze()
    tmp_all = np.nanmean(all_ecmf_ens[:,:,:,0:el,0:e+1,:,:],axis=(3,4)).squeeze()
    for l in np.arange(0,nleads):
      for y in np.arange(0,len(years)):
        for w in np.arange(0,nweeks):# loop through starts
          # Only use forecasts in DJF since data is masked with NaN's for different months:
          if week_mid_dates[y,l,w,1] in keep_mon:
            # Compute current anomalies for chirps
            curr_chirps_anom = all_chirps[y,l,w,:,:] - week_chirps_mean[l,w,:,:]
            curr_anom  = []
            curr_anom = tmp_all[y,l,w,:,:] - tmp_week_mean[l,w,:,:]
            acc_pre_ens_ecmf[0,l,e-1,:,:] = acc_pre_ens_ecmf[0,l,e-1,:,:] + (curr_anom*curr_chirps_anom)
            acc_pre_ens_ecmf[1,l,e-1,:,:] = acc_pre_ens_ecmf[1,l,e-1,:,:] + curr_anom**2
            acc_pre_ens_ecmf[2,l,e-1,:,:] = acc_pre_ens_ecmf[2,l,e-1,:,:] + curr_chirps_anom**2

acc_ens_ecmf[:,:,:,:] = acc_pre_ens_ecmf[0,:,:,:,:]/np.sqrt(acc_pre_ens_ecmf[1,:,:,:,:]*acc_pre_ens_ecmf[2,:,:,:,:])

col_ls = ['black','blue','green','orange','red']
#plt.subplot(1,2,1)
fig = plt.figure(figsize=(4,3))
for l in np.arange(0,5):
  #plt.plot(ensemble_bin+1,np.nanmean(acc_ens[l,:,:,:],axis=(1,2)),'--k',linewidth=2,color=col_ls[l],label='Week '+str(l+1))
  plt.plot(ensemble_bin_ecmf+1,np.nanmean(acc_ens_ecmf[l,:,:,:],axis=(1,2)),'-ok',linewidth=2,color=col_ls[l],label='_nolegend_')
#plt.ylim([0,0.5])
#plt.xlim([.7,7.3])
plt.ylabel('ACC')
plt.xlabel('Number of lags')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fname_plot = dir_out+region+'_S2S_weekly_lags_acc_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
plt.savefig(fname_plot+'.pdf',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()

# # alternatively compute bias using lagged ensemble means
# # (there is a difference to the above method but the difference is v. small)
# ukmo_bias_s	= np.nanmean(all_ukmo - all_chirps,axis=(0,2))
# ncep_bias_s	= np.nanmean(all_ncep - all_chirps,axis=(0,2))
# ecmf_bias_s	= np.nanmean(all_ecmf - all_chirps,axis=(0,2))

'''
Compute brier skill score for each tercile
'''

'''
Forecast probabilities
'''
percentiles = [100./3.,200./3.]
p_name = 'tercile' # give categories an appropriate name
p_names = ['dry','normal','wet']
p_cat = 0.333 # probability of falling into tercile

fname_p_chirps = dir_prob+region+'_chirps_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
fname_p_ukmo = dir_prob+region+'_UKMO_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
fname_p_ncep = dir_prob+region+'_NCEP_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
fname_p_ecmf = dir_prob+region+'_ECMWF_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'

fname_p_bam = dir_prob+region+'_BAM_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
fname_p_chirps_bam = dir_prob+region+'_CHIRPS_for_BAM_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'


if os.path.exists(fname_p_ukmo) == True:

  tmp_p_chirps = np.load(fname_p_chirps)
  tmp_p_ukmo = np.load(fname_p_ukmo)
  tmp_p_ncep = np.load(fname_p_ncep)
  tmp_p_ecmf = np.load(fname_p_ecmf)
  if season == 'NDJFM':
    tmp_p_chirps_for_bam = np.load(fname_p_chirps_bam)
    tmp_p_bam = np.load(fname_p_bam)
    p_chirps_for_bam = np.copy(tmp_p_chirps_for_bam)
    p_bam = np.copy(tmp_p_bam)
    n_week_p_bam = np.zeros((nleads))
    nan_week_p_bam = np.zeros((nleads))

  p_chirps = np.copy(tmp_p_chirps)
  p_ukmo = np.copy(tmp_p_ukmo)
  p_ncep = np.copy(tmp_p_ncep)
  p_ecmf = np.copy(tmp_p_ecmf)

  n_week_p = np.zeros((nleads))
  nan_week_p = np.zeros((nleads))

  #  mask the probabilities on where there are actually obs/fcsts!!!!
  for l in np.arange(0,nleads):
    for y in np.arange(0,len(years)):
      for w in np.arange(0,nweeks):# loop through starts
        if week_mid_dates[y,l,w,1] in keep_mon:
          n_week_p[l] = n_week_p[l] +1
          print 'lead time week '+str(l+1)+' mid-week date:'
        else:
          nan_week_p[l] = nan_week_p[l] + 1
          p_chirps[:,l,w,y,:,:] = np.nan
          p_ukmo[:,:,l,w,y,:,:] = np.nan
          p_ncep[:,:,l,w,y,:,:] = np.nan
          p_ecmf[:,:,l,w,y,:,:] = np.nan

      if season == 'NDJFM':
        for wbam in np.arange(0,5):# loop through starts
          for n in [0,1]:
            if week_mid_dates_bam[y,l,wbam,n,1] in keep_mon:
              n_week_p_bam[l] = n_week_p_bam[l] +1
              print 'lead time week '+str(l+1)+' mid-week date:'
            else:
              nan_week_p[l] = nan_week_p[l] + 1
              p_chirps_for_bam[:,l,wbam,n,y,:,:] = np.nan
              p_bam[:,:,l,wbam,n,y,:,:] = np.nan

  # BSS
  # nw = n_week_p[0] # number of data samples
  chirps_bs_clim = np.zeros((len(percentiles)+1,nleads,len(all_lat),len(all_lon)))
  if season == 'NDJFM':
    chirps_bs_clim_bam = np.zeros((len(percentiles)+1,nleads,len(all_lat),len(all_lon)))

  model_bs = np.zeros((4,len(percentiles)+1,2,nleads,len(all_lat),len(all_lon)))
  # brier score
  for l in np.arange(0,nleads):
    for ter in np.arange(0,len(percentiles)+1):
      chirps_bs_clim[ter,l,:,:] = np.nanmean((p_cat - p_chirps[ter,l,:,:,:,:])**2,axis=(0,1))#/n_week_p[l]
      if season == 'NDJFM':
        chirps_bs_clim_bam[ter,l,:,:] = np.nanmean((p_cat - p_chirps_for_bam[ter,l,:,:,:,:])**2,axis=(0,1,2))#/n_week_p_bam[l]

      for r in [0,1]:
        model_bs[3,ter,r,l,:,:] = np.nanmean((p_ukmo[ter,r,l,:,:,:,:] - p_chirps[ter,l,:,:,:,:])**2,axis=(0,1))#/n_week_p[l] # brier score sum
        model_bs[2,ter,r,l,:,:] = np.nanmean((p_ncep[ter,r,l,:,:,:,:] - p_chirps[ter,l,:,:,:,:])**2,axis=(0,1))#/n_week_p[l] # brier score sum
        model_bs[1,ter,r,l,:,:] = np.nanmean((p_ecmf[ter,r,l,:,:,:,:] - p_chirps[ter,l,:,:,:,:])**2,axis=(0,1))#/n_week_p[l] # brier score sum
        if season == 'NDJFM':
          model_bs[0,ter,r,l,:,:] = np.nanmean((p_bam[ter,r,l,:,:,:,:,:] - p_chirps_for_bam[ter,l,:,:,:,:,:])**2,axis=(0,1,2))#/n_week_p_bam[l] # brier score sum
        else:
          model_bs[0,ter,r,l,:,:] = np.nan
  # chirps_bs_clim = chirps_bs_clim/nw	# nw = 220
  model_bss = np.nan*np.zeros((4,len(percentiles)+1,nleads,len(all_lat),len(all_lon)))
  for m in [0,1,2,3]:
    if m == 0:
      if season == 'NDJFM':
        model_bss[m,:,:,:,:] = (chirps_bs_clim_bam - (model_bs[m,:,1,:,:,:])) / chirps_bs_clim_bam
    else:
      model_bss[m,:,:,:,:] = (chirps_bs_clim - (model_bs[m,:,1,:,:,:])) / chirps_bs_clim

# cnt=0
# for l in [0,4]:
#   for m in [0,1,2,3]:
#     cnt = cnt +1
#     plt.subplot(2,4,cnt)
#     plt.pcolor(all_lon,all_lat,model_bss[m,0,l,:,:]*dry_mask,vmin=-0.5,vmax=0.5,cmap='RdBu')
#     plt.colorbar()
# plt.show()
'''
Average ACC, RMSE and BIAS over  sub-regions at each lead time
'''
rg_names = ['NSA','AMZ','NDE','SESA','AND','PAT']
rg_lon_min = [-80,-67,-47,-60,-75,-75]
rg_lon_max = [-50,-47,-34,-48,-67,-60]
rg_lat_min = [0,-15,-15,-35,-40,-50]
rg_lat_max = [12,-5,-5,-22,-15,-40]
acc_region	= np.nan*np.zeros((len(rg_names),4,nleads))
rmse_region	= np.nan*np.zeros((len(rg_names),4,nleads))
rmse_anom_region	= np.nan*np.zeros((len(rg_names),4,nleads))
rmse_mean_region	= np.nan*np.zeros((len(rg_names),4,nleads))
bias_region	= np.nan*np.zeros((len(rg_names),4,nleads))
bss_region = np.nan*np.zeros((len(rg_names),4,3,nleads))
for rg in np.arange(0,len(rg_names)):
  subregion_lon_id = np.where((all_lon >= rg_lon_min[rg]) & (all_lon <= rg_lon_max[rg]))[0]
  subregion_lat_id = np.where((all_lat >= rg_lat_min[rg]) & (all_lat <= rg_lat_max[rg]))[0]
  acc_region[rg,:,:]=np.nanmean(acc[:,:,subregion_lat_id,:][:,:,:,subregion_lon_id],axis=(2,3))
  bias_region[rg,:,:]=np.nanmean(bias[:,:,subregion_lat_id,:][:,:,:,subregion_lon_id],axis=(2,3))
  rmse_region[rg,:,:]=np.nanmean(rmse[:,:,subregion_lat_id,:][:,:,:,subregion_lon_id],axis=(2,3))
  rmse_anom_region[rg,:,:]=np.nanmean(rmse_anom[:,:,subregion_lat_id,:][:,:,:,subregion_lon_id],axis=(2,3))
  rmse_mean_region[rg,:,:]=np.nanmean(rmse_mean[:,:,subregion_lat_id,:][:,:,:,subregion_lon_id],axis=(2,3))

  for m in np.arange(0,4):
    tmp_bss = []
    tmp_bss = model_bss[m,:,:,:,:]*dry_mask
    bss_region[rg,m,:,:] = np.nanmean(tmp_bss[:,:,subregion_lat_id,:][:,:,:,subregion_lon_id],axis=(2,3))



'''
Plot metrics
- regional averages at each lead time
- spatial maps of mean precipitation, model biases, ACC, BSS
'''
print 'Plotting figs...'
label_ls = ['(a)','(b)','(c)','(d)','(e)','(f)']
fmt_ls	= ['-s','-+','-x','-o']
ms_ls	= [6,10,10,8]
mw_ls	= [1,2,2,1]
col_ls	= ['black','forestgreen','dodgerblue','firebrick']
leads	= np.arange(1,nleads+1,1)
nrow = 2
ncol = 3

fname_plot = dir_out+region+'_S2S_weekly_ACC_'+season+'_subregions_average_'+str(years[0])+'_'+str(years[len(years)-1])
size1 = [8,4]
fig = plt.figure(figsize=(size1[0],size1[1]))
for rg in np.arange(0,len(rg_names)):
  plt.subplot(nrow,ncol,rg+1)
  for m in np.arange(0,len(model_ls)):
    # ms - markersize
    # mew - markeredgewidth
    # alpha - 0.0 transparent; 1.0 opaque
    plt.plot(leads,acc_region[rg,m,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.7,mew=mw_ls[m])
  plt.ylim([0,0.8])
  plt.xlim([0.8,5.2])
  plt.xticks(np.arange(1,nleads+1,1))
  plt.text(3.5,0.7,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
  if rg in [0,3]:
    plt.ylabel('ACC')
  if rg > 2:
    plt.xlabel('Lead time (weeks)')
  if rg == len(rg_names)-1:
    plt.legend(model_ls,loc=2,prop={'size':8})
plt.savefig(fname_plot+'.pdf',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()

fname_plot = dir_out+region+'_S2S_weekly_Bias_'+season+'_subregions_average_'+str(years[0])+'_'+str(years[len(years)-1])
size1 = [8,4]
fig = plt.figure(figsize=(size1[0],size1[1]))
for rg in np.arange(0,len(rg_names)):
  plt.subplot(nrow,ncol,rg+1)
  plt.plot(np.arange(0,7),np.zeros((7)),'-k',alpha=0.5,label='_nolegend_')
  for m in np.arange(0,len(model_ls)):
    plt.plot(leads,bias_region[rg,m,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.9,mew=mw_ls[m])
  plt.ylim([-2.5,2.5])
  plt.xlim([0.8,5.2])
  plt.xticks(np.arange(1,nleads+1,1))
  plt.text(3.5,2,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
  if rg in [0,3]:
    plt.ylabel('Bias (mm d$^{-1}$)')
  if rg > 2:
    plt.xlabel('Lead time (weeks)')
  if rg == len(rg_names)-1:
    plt.legend(model_ls,loc=4,prop={'size':8})
plt.savefig(fname_plot+'.pdf',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()

fname_plot = dir_out+region+'_S2S_weekly_RMSE_'+season+'_subregions_average_'+str(years[0])+'_'+str(years[len(years)-1])
size1 = [8,4]
fig = plt.figure(figsize=(size1[0],size1[1]))
for rg in np.arange(0,len(rg_names)):
  plt.subplot(nrow,ncol,rg+1)
  for m in np.arange(0,len(model_ls)):
    plt.plot(leads,rmse_region[rg,m,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.9,mew=mw_ls[m])
  plt.ylim([1.7,8.3])
  plt.xlim([0.8,5.2])
  plt.xticks(np.arange(1,nleads+1,1))
  plt.text(3.5,7.6,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
  if rg in [0,3]:
    plt.ylabel('RMSE (mm d$^{-1}$)')
  if rg > 2:
    plt.xlabel('Lead time (weeks)')
  if rg == len(rg_names)-2:
    plt.legend(model_ls,loc=2,prop={'size':8})
plt.savefig(fname_plot+'.pdf',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()

fname_plot = dir_out+region+'_S2S_weekly_RMSE_anom_'+season+'_subregions_average_'+str(years[0])+'_'+str(years[len(years)-1])
size1 = [8,4]
fig = plt.figure(figsize=(size1[0],size1[1]))
for rg in np.arange(0,len(rg_names)):
  plt.subplot(nrow,ncol,rg+1)
  for m in np.arange(0,len(model_ls)):
    plt.plot(leads,rmse_anom_region[rg,m,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.9,mew=mw_ls[m])
  plt.ylim([1.5,7])
  plt.xlim([0.8,5.2])
  plt.xticks(np.arange(1,nleads+1,1))
  plt.text(3.5,7.6,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
  if rg in [0,3]:
    plt.ylabel('RMSE (mm d$^{-1}$)')
  if rg > 2:
    plt.xlabel('Lead time (weeks)')
  if rg == len(rg_names)-2:
    plt.legend(model_ls,loc=2,prop={'size':8})
plt.savefig(fname_plot+'.pdf',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()

fname_plot = dir_out+region+'_S2S_weekly_RMSE_ensmean_'+season+'_subregions_average_'+str(years[0])+'_'+str(years[len(years)-1])
size1 = [8,4]
fig = plt.figure(figsize=(size1[0],size1[1]))
for rg in np.arange(0,len(rg_names)):
  plt.subplot(nrow,ncol,rg+1)
  for m in np.arange(0,len(model_ls)):
    plt.plot(leads,rmse_mean_region[rg,m,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.9,mew=mw_ls[m])
  if season == 'NDJFM':
    plt.ylim([0.2,2.1])
  elif season == 'MJJAS':
    plt.ylim([0.2,2.3])
  plt.xlim([0.8,5.2])
  plt.xticks(np.arange(1,nleads+1,1))
  plt.text(3.5,7.6,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
  if rg in [0,3]:
    plt.ylabel('RMSE (mm d$^{-1}$)')
  if rg > 2:
    plt.xlabel('Lead time (weeks)')
  if rg == len(rg_names)-2:
    plt.legend(model_ls,loc=2,prop={'size':8})
plt.savefig(fname_plot+'.pdf',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()

nrow = 2
ncol = 3
for ter in np.arange(0,len(p_names)):
  fname_plot = dir_out+region+'_S2S_weekly_BSS_tercile_'+p_names[ter]+'_'+season+'_subregions_average_'+str(years[0])+'_'+str(years[len(years)-1])
  size1 = [8,4]
  fig = plt.figure(figsize=(size1[0],size1[1]))
  for rg in np.arange(0,len(rg_names)):
    plt.subplot(nrow,ncol,rg+1)
    plt.plot(np.arange(0,7),np.zeros((7)),'-k',alpha=0.5,label='_nolegend_')
    for m in np.arange(0,len(model_ls)):
      plt.plot(leads,bss_region[rg,m,ter,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.9,mew=mw_ls[m])
    if (rg == 2) & (season=='MJJAS') & (ter in [0,2]):
      print 'plt'
    elif (rg == 1) & (season=='MJJAS') & (ter in [0]):
      print 'plt'
    else:
      plt.ylim([-0.6,0.3])
    plt.xlim([0.8,5.2])
    plt.xticks(np.arange(1,nleads+1,1))
    plt.text(3.5,0.17,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
    if rg in [0,3]:
      plt.ylabel('BSS')
    if rg > 2:
      plt.xlabel('Lead time (weeks)')
    if rg == len(rg_names)-3:
      plt.legend(model_ls,loc=0,prop={'size':8})
  plt.savefig(fname_plot+'.pdf',bbox_inches='tight',dpi=300)
  #plt.show()
  plt.close()



long_lead = np.arange(0,7)
for rg in np.arange(0,len(rg_names)):
  fname_plot = dir_out2+region+'_S2S_weekly_CC_'+season+'_region_average_'+rg_names[rg]
  fig = plt.figure(figsize=(2.5,2))
  plt.plot(long_lead,0.2*np.ones(len(long_lead)),'--k',alpha=0.5)
  for m in np.arange(0,len(model_ls)):
    # ms - markersize
    # mew - markeredgewidth
    # alpha - 0.0 transparent; 1.0 opaque
    plt.plot(leads,acc_region[rg,m,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.7,mew=mw_ls[m])
  plt.ylim([0,0.8])
  plt.xlim([0.8,5.2])
  plt.xticks(np.arange(1,nleads+1,1))
  #plt.text(3.5,0.7,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
  plt.ylabel('CC')
  plt.xlabel('Lead time (weeks)')
  plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
  #plt.show()
  plt.close()



for rg in np.arange(0,len(rg_names)):
  fname_plot = dir_out2+region+'_S2S_weekly_RMSE_'+season+'_region_average_'+rg_names[rg]
  fig = plt.figure(figsize=(2.5,2))
  for m in np.arange(0,len(model_ls)):
    # ms - markersize
    # mew - markeredgewidth
    # alpha - 0.0 transparent; 1.0 opaque
    plt.plot(leads,rmse_region[rg,m,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.7,mew=mw_ls[m])
  plt.ylim([1.5,7.5])
  plt.xlim([0.8,5.2])
  plt.xticks(np.arange(1,nleads+1,1))
  plt.yticks(np.arange(1,8,1))
  #plt.text(3.5,0.7,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
  plt.ylabel('RMSE (mm d$^{-1}$)')
  plt.xlabel('Lead time (weeks)')
  plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
  #plt.show()
  plt.close()

for rg in np.arange(0,len(rg_names)):
  fname_plot = dir_out2+region+'_S2S_weekly_RMSE_anom_'+season+'_region_average_'+rg_names[rg]
  fig = plt.figure(figsize=(2.5,2))
  for m in np.arange(0,len(model_ls)):
    # ms - markersize
    # mew - markeredgewidth
    # alpha - 0.0 transparent; 1.0 opaque
    plt.plot(leads,rmse_anom_region[rg,m,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.7,mew=mw_ls[m])
  plt.ylim([1.5,7])
  plt.xlim([0.8,5.2])
  plt.xticks(np.arange(1,nleads+1,1))
  plt.yticks(np.arange(1,8,1))
  #plt.text(3.5,0.7,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
  plt.ylabel('RMSE (mm d$^{-1}$)')
  plt.xlabel('Lead time (weeks)')
  plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
  plt.close()

for rg in np.arange(0,len(rg_names)):
  fname_plot = dir_out2+region+'_S2S_weekly_RMSE_ensmean_'+season+'_region_average_'+rg_names[rg]
  fig = plt.figure(figsize=(2.5,2))
  for m in np.arange(0,len(model_ls)):
    # ms - markersize
    # mew - markeredgewidth
    # alpha - 0.0 transparent; 1.0 opaque
    plt.plot(leads,rmse_mean_region[rg,m,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.7,mew=mw_ls[m])
  if season == 'NDJFM':
    plt.ylim([0.2,2.1])
    plt.yticks([0.2,0.5,1.0,1.5,2.0])
  elif season == 'MJJAS':
    plt.ylim([0.2,2.3])
    plt.yticks([0.2,0.5,1.0,1.5,2.0])
  plt.xlim([0.8,5.2])
  plt.xticks(np.arange(1,nleads+1,1))
  #plt.text(3.5,0.7,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
  plt.ylabel('RMSE (mm d$^{-1}$)')
  plt.xlabel('Lead time (weeks)')
  plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
  #plt.show()
  plt.close()

for ter in np.arange(0,len(p_names)):
  for rg in np.arange(0,len(rg_names)):
    fname_plot = dir_out2+region+'_S2S_weekly_BSS_tercile_'+p_names[ter]+'_'+season+'_region_average_'+rg_names[rg]
    fig = plt.figure(figsize=(2.5,2))
    plt.plot(np.arange(0,7),np.zeros((7)),'-k',alpha=0.5,label='_nolegend_')
    for m in np.arange(0,len(model_ls)):
      # ms - markersize
      # mew - markeredgewidth
      # alpha - 0.0 transparent; 1.0 opaque
      plt.plot(leads,bss_region[rg,m,ter,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.7,mew=mw_ls[m])
    if (rg == 2) & (season=='MJJAS') & (ter in [0,2]):
      print 'plt'
    elif (rg == 1) & (season=='MJJAS') & (ter in [0]):
      print 'plt'
    else:
      plt.ylim([-0.6,0.3])

    if ter == 2:
      plt.ylim([-0.3,0.33])
      plt.yticks(np.arange(-0.3,0.3+0.1,0.1))
    # plt.ylim([-0.6,0.3])
    plt.xlim([0.8,5.2])
    plt.xticks(np.arange(1,nleads+1,1))
    #plt.yticks(np.arange(plt.ylim()[0],plt.ylim()[1]+0.2,0.2))
    #plt.text(3.5,0.7,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
    plt.ylabel('BSS')
    plt.xlabel('Lead time (weeks)')
    plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
    #plt.show()
    plt.close()


# Plot seasonal mean precipitation
# web colours
# these are hexadecimal colour codes
# a hex triplet is a six-digit, three-byte hexadecimal number used in many computing applications to represent colours.
# bytes represent the red, green and blue components of the colour.
# one byte represent a number in the range 00 to FF (in hexadecimal notation) or 0 to 255 in decimal notation.
#precip_colors	= ["#fdfdfd","#f2f2f2","#bfbfbf","#04e9e7","#019ff4","#0300f4","#02fd02","#01c501","#008e00","#fdf802","#e5bc00","#fd9500","#fd0000", "#d40000","#bc0000","#f800fd","#9854c6"]
#precip_colormap	= matplotlib.colors.ListedColormap(precip_colors)
#cols		= precip_colormap
# a sequential colormap
# more details can be found https://matplotlib.org/tutorials/colors/colormaps.html
cols = 'PuBu'
cmin = 0
cmax = 16
cspc = 2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'Mean precipitation (mm d$^{-1}$)'
# matplotlib.colors.BoundaryNorm generates a colormap index based on discrete intervals.
# BoundaryNorm maps values to integers.
# BoundaryNorm defines the edges of bins, and data falling within a bin is mapped to the color with the same index.
# if the number of bins doesn't equal ncolors, the color is chosen by linear interpolation of the bin number onto color numbers.
norm = BoundaryNorm(boundaries=clevs,ncolors=256)
lw = 1
gl = 20

'''
size1 = [5,3.5]
fname_plot = dir_out+region+'_GPCC_weekly_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
mymap.drawcoastlines(linewidth=lw)
mymap.drawcountries(linewidth=lw)
# labels list of 4 values that control whether parallels are labelled where they intersect the left, right, top or bottom of the plot.
# +/-, north and south latitudes are labelled with "+" and "-", otherwise they are labelled with "N" and "S".
mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
x,y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
# create a psesudocolor plot with a non-regular rectuanglar grid
# norm scales the data values to the canonocal colormap range [0,1] for mapping to colors
uncal = mymap.pcolormesh(x,y,np.nanmean(week_chirps_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
plt.title('chirps '+season)
plt.colorbar(uncal,label=clabel1,extend='max')
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
plt.show()
plt.close()
'''
cols = 'PuBu'
cmin = 0
cmax = 16
cspc = 2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'Mean precipitation (mm d$^{-1}$)'
# matplotlib.colors.BoundaryNorm generates a colormap index based on discrete intervals.
# BoundaryNorm maps values to integers.
# BoundaryNorm defines the edges of bins, and data falling within a bin is mapped to the color with the same index.
# if the number of bins doesn't equal ncolors, the color is chosen by linear interpolation of the bin number onto color numbers.
norm = BoundaryNorm(boundaries=clevs,ncolors=256)
lw = 1
gl = 20
m_number=[6,11,16,21]
nrow = 5
ncol = nleads
size1 =[10,12]
fname_plot = dir_out+region+'_S2S_weekly_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
for m in [0,1,2,3]:
  for n in np.arange(0,nleads):
    if n == 2:
      plt.subplot(nrow,ncol,3)
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim2[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      x,y = mymap(*np.meshgrid(all_lon,all_lat))
  #   uncal = mymap.pcolormesh(x,y,np.nanmean(week_chirps_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      uncal = mymap.pcolormesh(x,y,np.nanmean(week_chirps_mean[0,:,:,:],axis=0),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title('CHIRPS '+season,fontweight='bold')

    plt.subplot(nrow,ncol,n+m_number[m])
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim2[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    if m == 3:
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x,y = mymap(*np.meshgrid(all_lon,all_lat))
  # uncal = mymap.pcolormesh(x,y,np.nanmean(week_ukmo_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    if m == 3:
      uncal = mymap.pcolormesh(x,y,np.nanmean(week_mean[m][n,:,:,:,:],axis=(0,1))*land_sea,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    else:
      uncal = mymap.pcolormesh(x,y,np.nanmean(week_mean[m][n,:,:,:],axis=0)*land_sea,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    #plt.text(-85,12,model_ls[m]+' Week '+str(n+1),fontsize=14)
    plt.title(model_ls[m]+' Week '+str(n+1))

#plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.5)
fig.subplots_adjust(right=0.90)
cbar_pos = [0.92, 0.05, 0.015, 0.15] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='max')
plt.savefig(fname_plot+'.pdf',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()


from matplotlib.patches import Polygon
def draw_screen_poly(lats,lons,m,lw):
    x, y = m(lons,lats)
    xy = zip(x,y)
    poly = Polygon(xy,linestyle='--',facecolor='none',edgecolor='firebrick',linewidth=lw)
    plt.gca().add_patch(poly)

'''
individual plot of CHIRPS
'''
cols = 'PuBu'
cmap = matplotlib.cm.PuBu
cmin = 0
cmax = 16
cspc = 2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'Mean precipitation (mm d$^{-1}$)'
norm = BoundaryNorm(boundaries=clevs,ncolors=256)
lw = 1
gl = 20
m_number=[6,11,16,21]
nrow = 5
ncol = nleads
size1 =[2,2]
fname_plot = dir_out2+region+'_CHIRPS_clim_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim2[1])
mymap.drawcoastlines(linewidth=lw)
#mymap.drawcountries(linewidth=lw)
mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
x,y = mymap(*np.meshgrid(all_lon,all_lat))
#   uncal = mymap.pcolormesh(x,y,np.nanmean(week_chirps_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
uncal = mymap.pcolormesh(x,y,np.nanmean(week_chirps_mean[0,:,:,:],axis=0),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)

# rg_names = ['NSA','AMZ','NDE','SESA','AND','PAT']
#
# x_coor = [60,41,]
# y_coor = [rg_lat_max[r]+1,-5,-18,-29,-50]
for r in np.arange(0,len(rg_names)):
  lat_box = [rg_lat_min[r],rg_lat_max[r],rg_lat_max[r],rg_lat_min[r]]
  lon_box = [rg_lon_min[r],rg_lon_min[r],rg_lon_max[r],rg_lon_max[r]]
  draw_screen_poly(lat_box,lon_box,mymap,1.8)
  #plt.text(x_coor[r],y_coor[r],rg_names[r],fontsize=10)
#plt.title('CHIRPS '+season,fontweight='bold')
#plt.colorbar(label=clabel1,extend='max')
plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()

fig = plt.figure(figsize=(1,2))
ax = fig.add_axes([0.4, 0.05, 0.07, 0.9])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='max',spacing='uniform',orientation='vertical',label=clabel1)
#cb.set_label(clabel)
plt.savefig(dir_out2+'CHIRPS_clim_colorbar.ps',bbox_inches='tight',dpi=300)
plt.close()

# cols = 'PuBu'
# cmap = matplotlib.cm.PuBu
cols = 'Reds'
cmap= matplotlib.cm.Reds
cmin = 0
cmax = 12
cspc = 2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'S.D. precipitation (mm d$^{-1}$)'
norm = BoundaryNorm(boundaries=clevs,ncolors=256)
lw = 1
gl = 20
m_number=[6,11,16,21]
nrow = 5
ncol = nleads
size1 =[2,2]
fname_plot = dir_out2+region+'_CHIRPS_clim_std_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim2[1])
mymap.drawcoastlines(linewidth=lw)
#mymap.drawcountries(linewidth=lw)
mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
x,y = mymap(*np.meshgrid(all_lon,all_lat))
#   uncal = mymap.pcolormesh(x,y,np.nanmean(week_chirps_mean,axis=(0,1)),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
uncal = mymap.pcolormesh(x,y,np.nanstd(week_chirps_mean[0,:,:,:],axis=0),vmin=cmin,vmax=cmax,cmap=cols,norm=norm)

# rg_names = ['NSA','AMZ','NDE','SESA','AND','PAT']
#
# x_coor = [60,41,]
# y_coor = [rg_lat_max[r]+1,-5,-18,-29,-50]
# for r in np.arange(0,len(rg_names)):
#   lat_box = [rg_lat_min[r],rg_lat_max[r],rg_lat_max[r],rg_lat_min[r]]
#   lon_box = [rg_lon_min[r],rg_lon_min[r],rg_lon_max[r],rg_lon_max[r]]
#   draw_screen_poly(lat_box,lon_box,mymap,1.8)
  #plt.text(x_coor[r],y_coor[r],rg_names[r],fontsize=10)
#plt.title('CHIRPS '+season,fontweight='bold')
#plt.colorbar(label=clabel1,extend='max')
plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()

fig = plt.figure(figsize=(1,2))
ax = fig.add_axes([0.4, 0.05, 0.07, 0.9])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='max',spacing='uniform',orientation='vertical',label=clabel1)
#cb.set_label(clabel)
plt.savefig(dir_out2+'CHIRPS_clim_std_colorbar.ps',bbox_inches='tight',dpi=300)
plt.close()

fig = plt.figure(figsize=(3,1))
ax = fig.add_axes([0.05, 0.5, 0.9, 0.07])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='max',spacing='uniform',orientation='horizontal',label=clabel1)
#cb.set_label(clabel)
plt.savefig(dir_out2+'CHIRPS_clim_std_colorbar_horizontal.ps',bbox_inches='tight',dpi=300)
plt.close()

'''
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
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
plt.show()
plt.close()
'''



# plot ACC
cols = 'RdYlBu'
cmin = -1
cmax = 1
cspc = 0.2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'CC'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
m_number=[1,6,11,16]
lw = 1
nrow = len(model_ls)
ncol = nleads
size1=[8,7.5]
gl = 20
fname_plot = dir_out+region+'_S2S_weekly_CC_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
for m in np.arange(0,len(model_ls)):
  for n in np.arange(0,nleads):
    plt.subplot(nrow,ncol,n+m_number[m])
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    if m == 3:
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(all_lon,all_lat))
    uncal = mymap.pcolormesh(x,y,acc[m,n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title(model_ls[m]+' Week '+str(n+1))

#plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
plt.tight_layout(pad=2,w_pad=0.02,h_pad=0.02)
fig.subplots_adjust(bottom=0.08)
cbar_pos = [0.25, 0.03, 0.5, 0.015] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,orientation='horizontal',extend='both')
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
nrow = len(model_ls)
ncol = nleads
gl = 20
fname_plot = dir_out+region+'_S2S_weekly_Bias_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(8,7.5))
for m in np.arange(0,len(model_ls)):
  for n in np.arange(0,nleads):
    plt.subplot(nrow,ncol,n+m_number[m])#indexes goes from 1 to nrows * ncols
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    if m == 3:
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x,y = mymap(*np.meshgrid(all_lon,all_lat))
    tmp_bias = []
    tmp_bias = np.copy(bias[m,n,:,:])
    tmp_bias[tmp_bias == 0] = np.nan
    uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title(model_ls[m]+' Week '+str(n+1))
#plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.5)
fig.subplots_adjust(bottom=0.08)
cbar_pos = [0.25, 0.03, 0.5, 0.015] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,orientation='horizontal',extend='both')
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
#plt.show()
plt.close()


# plot Bias
cols = 'RdBu'
cmap = matplotlib.cm.RdBu
cmin = -10
cmax = 10
cspc = 2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'Bias (mm d$^{-1}$)'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
lw = 1
nrow = len(model_ls)
ncol = nleads
gl = 20
for m in np.arange(0,len(model_ls)):
  for n in [0,1,2,3,4]:
    fname_plot = dir_out2+region+'_S2S_weekly_Bias_'+season+'_'+model_ls[m]+'_Week'+str(n+1)
    fig = plt.figure(figsize=(2,2))
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x,y = mymap(*np.meshgrid(all_lon,all_lat))
    tmp_bias = []
    tmp_bias = np.copy(bias[m,n,:,:])
    tmp_bias[tmp_bias == 0] = np.nan
    uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    # plt.title(model_ls[m]+' Week '+str(n+1))
    # #plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
    # plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.5)
    # fig.subplots_adjust(bottom=0.08)
    # cbar_pos = [0.25, 0.03, 0.5, 0.015] #[left, bottom, width, height]
    # cbar_ax = fig.add_axes(cbar_pos)
    # cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,orientation='horizontal',extend='both')
    plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
    #plt.show()
    plt.close()

fig = plt.figure(figsize=(1,2))
ax = fig.add_axes([0.4, 0.05, 0.07, 0.9])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='both',spacing='uniform',orientation='vertical',label=clabel1)
#cb.set_label(clabel)
plt.savefig(dir_out2+'Bias_colorbar_vertical.ps',bbox_inches='tight',dpi=300)
plt.close()

fig = plt.figure(figsize=(3,1))
ax = fig.add_axes([0.05, 0.5, 0.9, 0.07])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='both',spacing='uniform',orientation='horizontal',label=clabel1)
#cb.set_label(clabel)
plt.savefig(dir_out2+'Bias_colorbar_horizontal.ps',bbox_inches='tight',dpi=300)
plt.close()


# plot ACC
cols = 'RdYlBu'
cmap= matplotlib.cm.RdYlBu
cmin = -1
cmax = 1
cspc = 0.2
clevs = np.arange(cmin,cmax+cspc,cspc)
# clevs = np.array([-1.0,-0.6,-0.2,0,0.2,0.6,1.0])
clabel1 = 'CC'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
m_number=[1,6,11,16]
lw = 1
nrow = len(model_ls)
ncol = nleads
gl = 20

for m in np.arange(0,len(model_ls)):
  for n in [0,1,2,3,4]:
    fname_plot = dir_out2+region+'_S2S_weekly_CC_'+season+'_'+model_ls[m]+'_Week'+str(n+1)
    fig = plt.figure(figsize=(2,2))
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(all_lon,all_lat))
    uncal = mymap.pcolormesh(x,y,acc[m,n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    # con = mymap.contour(x,y,acc[m,n,:,:],[0.2],colors='grey',linewidth=1,linestyle='dashed')
    plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
    #plt.show()
    plt.close()


fig = plt.figure(figsize=(1,2))
ax = fig.add_axes([0.4, 0.05, 0.07, 0.9])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,spacing='uniform',orientation='vertical',label=clabel1)
cb.set_ticks([-1.0,-0.6,-0.2,0.2,0.6,1.0])
cb.set_ticklabels([-1.0,-0.6,-0.2,0.2,0.6,1.0])
#cb.set_label(clabel)
plt.savefig(dir_out2+'CC_colorbar_vertical.ps',bbox_inches='tight',dpi=300)
plt.close()

fig = plt.figure(figsize=(3,1))
ax = fig.add_axes([0.05, 0.5, 0.9, 0.07])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,spacing='uniform',orientation='horizontal',label=clabel1)
cb.set_ticks([-1.0,-0.6,-0.2,0.2,0.6,1.0])
cb.set_ticklabels([-1.0,-0.6,-0.2,0.2,0.6,1.0])
#cb.set_label(clabel)
plt.savefig(dir_out2+'CC_colorbar_horizontal.ps',bbox_inches='tight',dpi=300)
plt.close()


# plot RMSE
cols = 'Reds'
cmap= matplotlib.cm.Reds
cmin = 0
cmax = 14
cspc = 2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'RMSE (mm d$^{-1}$)'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
m_number=[1,6,11,16]
lw = 1
nrow = len(model_ls)
ncol = nleads
gl = 20

for m in np.arange(0,len(model_ls)):
  for n in [0,1,2,3,4]:
    fname_plot = dir_out2+region+'_S2S_weekly_RMSE_'+season+'_'+model_ls[m]+'_Week'+str(n+1)
    fig = plt.figure(figsize=(2,2))
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(all_lon,all_lat))
    uncal = mymap.pcolormesh(x,y,rmse[m,n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
    #plt.show()
    plt.close()

fig = plt.figure(figsize=(1,2))
ax = fig.add_axes([0.4, 0.05, 0.07, 0.9])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,spacing='uniform',orientation='vertical',label=clabel1,extend='max')
#cb.set_label(clabel)
plt.savefig(dir_out2+'RMSE_colorbar_vertical.ps',bbox_inches='tight',dpi=300)
plt.close()

fig = plt.figure(figsize=(3,1))
ax = fig.add_axes([0.05, 0.5, 0.9, 0.07])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,spacing='uniform',orientation='horizontal',label=clabel1,extend='max')
#cb.set_label(clabel)
plt.savefig(dir_out2+'RMSE_colorbar_horizontal.ps',bbox_inches='tight',dpi=300)
plt.close()


# plot RMSE
cols = 'Reds'
cmap= matplotlib.cm.Reds
cmin = 0
cmax = 12
cspc = 2
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'RMSE (mm d$^{-1}$)'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
m_number=[1,6,11,16]
lw = 1
nrow = len(model_ls)
ncol = nleads
gl = 20

for m in np.arange(0,len(model_ls)):
  for n in [0,1,2,3,4]:
    fname_plot = dir_out2+region+'_S2S_weekly_RMSE_anom_'+season+'_'+model_ls[m]+'_Week'+str(n+1)
    fig = plt.figure(figsize=(2,2))
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(all_lon,all_lat))
    uncal = mymap.pcolormesh(x,y,rmse_anom[m,n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
    #plt.show()
    plt.close()

fig = plt.figure(figsize=(1,2))
ax = fig.add_axes([0.4, 0.05, 0.07, 0.9])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,spacing='uniform',orientation='vertical',label=clabel1,extend='max')
#cb.set_label(clabel)
plt.savefig(dir_out2+'RMSE_anom_colorbar_vertical.ps',bbox_inches='tight',dpi=300)
plt.close()

fig = plt.figure(figsize=(3,1))
ax = fig.add_axes([0.05, 0.5, 0.9, 0.07])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,spacing='uniform',orientation='horizontal',label=clabel1,extend='max')
#cb.set_label(clabel)
plt.savefig(dir_out2+'RMSE_anom_colorbar_horizontal.ps',bbox_inches='tight',dpi=300)
plt.close()

# plot RMSE
cols = 'Reds'
cmap= matplotlib.cm.Reds
cmin = 0
cmax = 5
cspc = 0.5
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'RMSE (mm d$^{-1}$)'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
m_number=[1,6,11,16]
lw = 1
nrow = len(model_ls)
ncol = nleads
gl = 20

for m in np.arange(0,len(model_ls)):
  for n in [0,1,2,3,4]:
    fname_plot = dir_out2+region+'_S2S_weekly_RMSE_ensmean_'+season+'_'+model_ls[m]+'_Week'+str(n+1)
    fig = plt.figure(figsize=(2,2))
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(all_lon,all_lat))
    uncal = mymap.pcolormesh(x,y,rmse_mean[m,n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
    #plt.show()
    plt.close()

fig = plt.figure(figsize=(1,2))
ax = fig.add_axes([0.4, 0.05, 0.07, 0.9])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,spacing='uniform',orientation='vertical',label=clabel1,extend='max')
#cb.set_label(clabel)
plt.savefig(dir_out2+'RMSE_ensmean_colorbar_vertical.ps',bbox_inches='tight',dpi=300)
plt.close()

fig = plt.figure(figsize=(3,1))
ax = fig.add_axes([0.05, 0.5, 0.9, 0.07])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,spacing='uniform',orientation='horizontal',label=clabel1,extend='max')
#cb.set_label(clabel)
plt.savefig(dir_out2+'RMSE_ensmean_colorbar_horizontal.ps',bbox_inches='tight',dpi=300)
plt.close()

# plot BSS
cols = 'RdYlBu'
cmap= matplotlib.cm.RdYlBu
cmin = -0.5
cmax = 0.5
cspc = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'BSS'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
m_number=[1,6,11,16]
lw = 1
nrow = len(model_ls)
ncol = nleads
gl = 20

for ter in np.arange(0,len(p_names)):
  for m in np.arange(0,len(model_ls)):
    for n in [0,1,2,3,4]:
      fname_plot = dir_out2+region+'_S2S_weekly_BSS_tercile_'+p_names[ter]+'_'+season+'_'+model_ls[m]+'_Week'+str(n+1)
      fig = plt.figure(figsize=(2,2))
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(all_lon,all_lat))
      uncal = mymap.pcolormesh(x,y,model_bss[m,ter,n,:,:]*dry_mask,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
      #plt.show()
      plt.close()

fig = plt.figure(figsize=(1,2))
ax = fig.add_axes([0.4, 0.05, 0.07, 0.9])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='both',spacing='uniform',orientation='vertical',label=clabel1)
#cb.set_label(clabel)
plt.savefig(dir_out2+'BSS_colorbar_vertical.ps',bbox_inches='tight',dpi=300)
plt.close()

fig = plt.figure(figsize=(3,1))
ax = fig.add_axes([0.05, 0.5, 0.9, 0.07])
cb = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,extend='both',spacing='uniform',orientation='horizontal',label=clabel1)
#cb.set_label(clabel)
plt.savefig(dir_out2+'BSS_colorbar_horizontal.ps',bbox_inches='tight',dpi=300)
plt.close()


# plot ACC
cols = 'RdYlBu'
cmin = -0.5
cmax = 0.5
cspc = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'BSS'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
m_number=[1,6,11,16]
lw = 1
nrow = len(model_ls)
ncol = nleads
size1=[8,7.5]
gl = 20

for ter in np.arange(0,len(p_names)):
  fname_plot = dir_out+region+'_S2S_weekly_BSS_tercile_'+p_names[ter]+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
  fig = plt.figure(figsize=(size1[0],size1[1]))
  for m in np.arange(0,len(model_ls)):
    for n in np.arange(0,nleads):
      plt.subplot(nrow,ncol,n+m_number[m])
      mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1])
      mymap.drawcoastlines(linewidth=lw)
      #mymap.drawcountries(linewidth=lw)
      mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
      mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
      if n in [0]:
        mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
      if m == 3:
        mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
      x, y = mymap(*np.meshgrid(all_lon,all_lat))
      uncal = mymap.pcolormesh(x,y,model_bss[m,ter,n,:,:]*dry_mask,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
      plt.title(model_ls[m]+' Week '+str(n+1))

  #plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
  plt.tight_layout(pad=2,w_pad=0.02,h_pad=0.02)
  fig.subplots_adjust(bottom=0.08)
  cbar_pos = [0.25, 0.03, 0.5, 0.015] #[left, bottom, width, height]
  cbar_ax = fig.add_axes(cbar_pos)
  cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,orientation='horizontal',extend='both')
  plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
  #plt.show()
  plt.close()


'''
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
'''


'''
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
fname_plot = dir_out+region+'_S2S_weekly_RMSE_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(size1[0],size1[1]))
for n in np.arange(0,nleads):
  plt.subplot(nrow,ncol,n+1)
  mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  #mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
  tmp_bias = []
  tmp_bias = np.copy(ukmo_rmse[n,:,:])
  tmp_bias[tmp_bias == 0] = np.nan
  uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('UKMO Week '+str(n+1))

  plt.subplot(nrow,ncol,n+6)
  mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  #mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x, y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  tmp_bias = []
  tmp_bias = np.copy(ncep_rmse[n,:,:])
  tmp_bias[tmp_bias == 0] = np.nan
  uncal = mymap.pcolormesh(x,y,tmp_bias,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('NCEP Week '+str(n+1))

  plt.subplot(nrow,ncol,n+11)
  mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  #mymap.drawcountries(linewidth=lw)
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

plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
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
size1 = [15,9]
gl = 20

cmin = -0.5
cmax = 0.5
cspc = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
norm = BoundaryNorm(boundaries=clevs,ncolors=256)

for tc in [0,1,2]:
  fname_plot = dir_out+region+'_S2S_weekly_BSS_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'_'+tercile_name[tc]
  fig = plt.figure(figsize=(size1[0],size1[1]))

  for n in np.arange(0,nleads):
    plt.subplot(nrow,ncol,n+1)
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x,y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
    mask_array = np.ma.array(ukmo_bss[tc,n,:,:]*dry_mask, mask=np.isnan(dry_mask))
    uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title('UKMO Week '+str(n+1))

    plt.subplot(nrow,ncol,n+6)
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    #mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
    mask_array = np.ma.array(ncep_bss[tc,n,:,:]*dry_mask, mask=np.isnan(dry_mask))
    uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title('NCEP Week '+str(n+1))

    plt.subplot(nrow,ncol,n+11)
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
    if n in [0]:
      mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
    mask_array = np.ma.array(ecmf_bss[tc,n,:,:]*dry_mask, mask=np.isnan(dry_mask))
    uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    plt.title('ECMWF Week '+str(n+1))

  plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
  plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.4)
  fig.subplots_adjust(right=0.90)
  cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]

  cbar_ax = fig.add_axes(cbar_pos)
  cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
  plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
# plt.show()
  plt.close()

'''

'''
fname_plot = dir_out+region+'_S2S_weekly_BSS_BiasCorrected_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'_'+tercile_name[tc]
fig = plt.figure(figsize=(size1[0],size1[1]))
for n in np.arange(0,nleads):
  plt.subplot(nrow,ncol,n+1)
  mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x,y = mymap(*np.meshgrid(ukmo_lon,ukmo_lat))
  mask_array = np.ma.array(ukmo_bss_bc[tc,n,bb,:,:]*dry_mask, mask=np.isnan(dry_mask))
  uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('UKMO Week '+str(n+1))

  plt.subplot(nrow,ncol,n+6)
  mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  #uncal = mymap.contourf(x,y,ukmo_acc[n,:,:],clevs,cmap=cols,extend='both')
  mask_array = np.ma.array(ncep_bss_bc[tc,n,bb,:,:]*dry_mask, mask=np.isnan(dry_mask))
  uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('NCEP Week '+str(n+1))

  plt.subplot(nrow,ncol,n+11)
  mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim[0],urcrnrlon=lonlim[1])
  mymap.drawcoastlines(linewidth=lw)
  mymap.drawcountries(linewidth=lw)
  mymap.drawparallels(np.arange(-90,90,gl),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,0],labelstyle='+/-')
  if n in [0]:
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
  x,y = mymap(*np.meshgrid(ncep_lon,ncep_lat))
  mask_array = np.ma.array(ecmf_bss_bc[tc,n,bb,:,:]*dry_mask, mask=np.isnan(dry_mask))
  uncal = mymap.pcolormesh(x,y,mask_array,vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title('ECMWF Week '+str(n+1))

plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.4)
fig.subplots_adjust(right=0.90)
cbar_pos = [0.92, 0.22, 0.015, 0.55] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='both')
plt.savefig(fname_plot+'.pdf',bbox_inches='tight')
plt.show()
plt.close()
'''

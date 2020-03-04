'''
Analyse skill of weekly S2S precipitation forecasts
from UKMO, NCEP and ECMWF over a given region

1) Reads in forecasts data saved over SAmerica
using script 'read_hindcasts_for_analysis.py'

2) Computes mean precipitation, bias, RMSE,
anomaly correlation coefficient and brier skill scores

M. Young 03/01/2020
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
from scipy.stats import pearsonr
from scipy.stats import ttest_ind_from_stats
import scipy.stats as stats
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


def t_test_corr(corr,n,n_effective):
  '''
  Function 't_test_corr':
  Compute the two-tailed t-statistic on a correlation
  for a given sample size

  inputs (3):
  - corr = correlation coefficient to test
  - n = original sample size
  - n_effective = adjusted sample size (adjusted for e.g. autocorrelation)

  outputs (4):
  - output1 = 1 if t-test statistically significant at 95% level,  = 0 if not significant, tested using n_effective.
  - output2 = as output1 except tested using n.
  - t_stat = t-statistic computed for the correlation using n_effective
  - t95_val = t-statistic at 95% level for n_effective sample size, using the t-distribution

  M. Young 02/03/2020
  '''
  # compute t-statistic using n samples and correlations
  t_stat = corr*np.sqrt((n_effective-2)/(1-(corr**2)))
  t_stat_original = corr*np.sqrt((n-2)/(1-(corr**2)))
  #Student-t, p<0.05, 2-tail (have to take half p-value for 2-tail,ie.0.025)
  t95_val = stats.t.ppf(1-0.025, df=n_effective-2)
  t95_val_original = stats.t.ppf(1-0.025, df=n-2)
  if t_stat > t95_val:
    output1 = 1
  else:
    output1 = np.nan
  if t_stat_original > t95_val_original:
    output2 = 1
  else:
    output2 = np.nan
  return(output1,output2,t_stat,t95_val)

# dir_in	= '/group_workspaces/jasmin2/ncas_climate/users/myoung02/datasets/S2S_forecasts/weekly_dubstep_style/'
dir_in = '/gws/nopw/j04/ncas_climate_vol1/users/myoung02/datasets/DUBSTEP_data/'
dir_prob	= dir_in+'hindcast_probabilities/'
# dir_out	= '/home/users/myoung02/DUBSTEP_paper_results/'
dir_out	= '/home/users/myoung02/DUBSTEP_paper_results_revisions_March2020/'
dir_out2 = dir_out+'individual_panels/' # extra directory for single panel plots

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
          if week_mid_dates_bam[y,l,wbam,n,1] in keep_mon:
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
                rmse_anom[0,l,:,:] = rmse_anom[0,l,:,:] + (curr_bam_anom - curr_chirps_anom)**2

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
        anom_ecmf[y,l,w,:,:] = np.copy(curr_ecmf_anom)
        acc_pre[1,0,l,:,:] = acc_pre[1,0,l,:,:] + (curr_ecmf_anom*curr_chirps_anom)
        acc_pre[1,1,l,:,:] = acc_pre[1,1,l,:,:] + curr_ecmf_anom**2
        acc_pre[1,2,l,:,:] = acc_pre[1,2,l,:,:] + curr_chirps_anom**2
        rmse_mean[1,l,:,:] = rmse_mean[1,l,:,:] + (all_ecmf[y,l,w,:,:]-all_chirps[y,l,w,:,:])**2

        curr_ncep_anom = []
        curr_ncep_anom = all_ncep[y,l,w,:,:] - week_ncep_mean[l,w,:,:]
        anom_ncep[y,l,w,:,:] = np.copy(curr_ncep_anom)
        acc_pre[2,0,l,:,:] = acc_pre[2,0,l,:,:] + (curr_ncep_anom*curr_chirps_anom)
        acc_pre[2,1,l,:,:] = acc_pre[2,1,l,:,:] + curr_ncep_anom**2
        acc_pre[2,2,l,:,:] = acc_pre[2,2,l,:,:] + curr_chirps_anom**2
        rmse_mean[2,l,:,:] = rmse_mean[2,l,:,:] + (all_ncep[y,l,w,:,:]-all_chirps[y,l,w,:,:])**2

        # for anomaly corelation coefficent just compare lagged, ensemble means
        curr_ukmo_anom = []
        curr_ukmo_anom = all_ukmo[y,l,w,:,:] - week_ukmo_mean[l,w,:,:]
        anom_ukmo[y,l,w,:,:] = np.copy(curr_ukmo_anom)
        acc_pre[3,0,l,:,:] = acc_pre[3,0,l,:,:] + (curr_ukmo_anom*curr_chirps_anom)
        acc_pre[3,1,l,:,:] = acc_pre[3,1,l,:,:] + curr_ukmo_anom**2
        acc_pre[3,2,l,:,:] = acc_pre[3,2,l,:,:] + curr_chirps_anom**2
        rmse_mean[3,l,:,:] = rmse_mean[3,l,:,:] + (all_ukmo[y,l,w,:,:]-all_chirps[y,l,w,:,:])**2

      else: # if not DJF
        week_bias[0:4,l,w,:,:] = np.nan
        nsamp_week[:,l,w] = np.nan

# Compute final ACC
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


all_anom = np.stack((anom_chirps,anom_ecmf,anom_ncep,anom_ukmo))
# col_ls = ['black','dodgerblue','forestgreen','firebrick']
# m_ls = ['CHIRPS','ECMWF', 'NCEP', 'UKMO']
# for a in np.arange(0,4):
#   tmp_anom = all_anom[a,:,:,:,:,:].flatten()
#   nan_id = np.where(np.isnan(tmp_anom)==False)[0]
#   plt.hist(tmp_anom[nan_id],bins=np.arange(-30,32,2),color=col_ls[a],histtype='step',label=m_ls[a])
# plt.legend(loc=0)
# plt.show()


# Compute the autocorrelation
all_cov = np.nan*np.zeros((4,len(years),nleads,len(all_lat),len(all_lon)))
all_cov_mean = np.nan*np.zeros((4,nleads,len(all_lat),len(all_lon)))
all_sd = np.nan*np.zeros((2,4,nleads,len(all_lat),len(all_lon)))
auto_corr = np.nan*np.zeros((4,nleads,len(all_lat),len(all_lon)))
n_eff = np.nan*np.zeros((4,nleads,len(all_lat),len(all_lon)))
t_signif = np.nan*np.zeros((4,4,nleads,len(all_lat),len(all_lon)))
# 1) compute lagged co-variance for each year
# 2) average lagged co-variance across all years
# 3) compute correlation: divide the mean lagged co-variance (from 2)
#    by the products of the standard deviations of the two full (all years)
#    time series (lagged and not lagged)
for p in np.arange(0,3):
  print p
  for l in np.arange(0,nleads):
    for i in np.arange(0,len(all_lat)):
      for j in np.arange(0,len(all_lon)):
        for y in np.arange(0,len(years)):
          nan_id = np.where(np.isnan(all_anom[p+1,y,l,:,i,j])==False)[0]
          x = all_anom[p+1,y,l,nan_id,i,j]
          all_cov[p+1,y,l,i,j] = np.cov(x[:-1],x[1:])[0][1]

          if p == 0:
            tmp_bam = np.copy(anom_bam[y,l,:,:,i,j].flatten())
            nan_id = np.where(np.isnan(tmp_bam)==False)[0]
            yy = tmp_bam[nan_id]
            all_cov[p,y,l,i,j] = np.cov(yy[:-1],yy[1:])[0][1]
            # tmp_chp = np.copy(anom_chirps_bam[y,l,:,:,i,j].flatten())
            # nan_id = np.where(np.isnan(tmp_chp)==False)[0]
            # yy = tmp_chp[nan_id]
            # all_cov[5,y,l,i,j] = np.cov(yy[:-1],yy[1:])[0][1]

        all_sd[0,p+1,l,i,j] = np.nanstd(all_anom[p+1,:,l,:-1,i,j])
        all_sd[1,p+1,l,i,j] = np.nanstd(all_anom[p+1,:,l,1:,i,j])
        all_cov_mean[p+1,l,i,j] = np.mean(all_cov[p+1,:,l,i,j])
        auto_corr[p+1,l,i,j] = all_cov_mean[p+1,l,i,j]/(all_sd[0,p+1,l,i,j]*all_sd[1,p+1,l,i,j])
        n_eff[p+1,l,i,j] = nsamp_week_only[l]*(1-auto_corr[p+1,l,i,j])/(1+auto_corr[p+1,l,i,j])
        t_signif[:,p+1,l,i,j] = t_test_corr(acc[p+1,l,i,j],nsamp_week_only[l],n_eff[p+1,l,i,j])

        if p == 0:
          tmp_bam = []
          tmp_bam = np.copy(anom_bam[:,l,:,:,i,j].flatten())
          all_sd[0,p,l,i,j] = np.nanstd(tmp_bam[:-1])
          all_sd[1,p,l,i,j] = np.nanstd(tmp_bam[1:])
          all_cov_mean[p,l,i,j] = np.mean(all_cov[p,:,l,i,j])
          auto_corr[p,l,i,j] = all_cov_mean[p,l,i,j]/(all_sd[0,p,l,i,j]*all_sd[1,p,l,i,j])
          n_eff[p,l,i,j] = nsamp_week_only_bam[l]*(1-auto_corr[p,l,i,j])/(1+auto_corr[p,l,i,j])
          t_signif[:,0,l,i,j] = t_test_corr(acc[0,l,i,j],nsamp_week_only_bam[l],n_eff[p,l,i,j])

          # tmp_chp = []
          # tmp_chp = np.copy(anom_chirps_bam[:,l,:,:,i,j].flatten())
          # all_sd[0,5,l,i,j] = np.nanstd(tmp_chp[:-1])
          # all_sd[1,5,l,i,j] = np.nanstd(tmp_chp[1:])
          # all_cov_mean[5,l,i,j] = np.mean(all_cov[5,:,l,i,j])
          # auto_corr[5,l,i,j] = all_cov_mean[5,l,i,j]/(all_sd[0,5,l,i,j]*all_sd[1,5,l,i,j])
          # n_eff[5,l,i,j] = nsamp_week_only_bam[l]*(1-auto_corr[5,l,i,j])/(1+auto_corr[5,l,i,j])

'''
###
# check dependance of skill on ensemble size
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
'''


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

'''
Save evaluation statistics
'''
print 'saving stats...'
# save mean
fname_mean_chp = dir_in+'CHIRPS_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_mean_chp,week_chirps_mean)

# save bias
fname_bias = dir_in+'Bias_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_bias,bias)
# save ACC
fname_acc = dir_in+'ACC_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_acc,acc)
# save significance
fname_t = dir_in+'ttest_significance_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_t,t_signif)
# save RMSE
fname_rmse = dir_in+'RMSE_anom_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_rmse,rmse_anom)
# save RMSE
fname_rmse1 = dir_in+'RMSE_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_rmse1,rmse)
# save RMSE
fname_rmse2 = dir_in+'RMSE_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_rmse2,rmse_mean)
# save BSS
fname_bss = dir_in+'BSS_tercile_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_bss,model_bss)

week_ukmo_mean = np.nanmean(all_ukmo_ens,axis=(0,3))# (5, 20, 40, 47)
week_ncep_mean = np.nanmean(all_ncep_ens,axis=(0,3,4))
week_ecmf_mean = np.nanmean(all_ecmf_ens,axis=(0,3,4))
week_bam_mean = np.nanmean(all_bam_ens,axis=(0,4))
# save mean
fname_ukmo = dir_in+'UKMO_week_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_ukmo,week_ukmo_mean)
fname_ncep = dir_in+'NCEP_week_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_ncep,week_ncep_mean)
fname_emcf = dir_in+'ECMWF_week_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_ecmf,week_ecmf_mean)
fname_bam = dir_in+'BAM_week_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
np.save(fname_bam,week_bam_mean)

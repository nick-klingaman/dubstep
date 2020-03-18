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

def read_hindcasts(season):
  '''
  read_hindcasts
  M. Young, April 2019
  '''
  dir_in = '/gws/nopw/j04/ncas_climate_vol1/users/myoung02/datasets/DUBSTEP_data/'
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
  if season == 'NDJFM':
    # for DJF, prepend five weeks starting 1025
    starts_ukmo = ['0901','0909','0917','0925',
                   '1001','1009','1017','1025',
                   '1101','1109','1117','1125',
                   '1201','1209','1217','1225',
                   '0101','0109','0117','0125',
                   '0201','0209','0217','0225',
                   '0301','0309','0317','0325']
    starts_mon = ['11','12','01','02','03']
    starts_day = ['01','09','17','25']
    starts_mon_bam = ['11','12','01','02','03']
    keep_mon = [11,12,1,2,3]
    end_year_starts_mon = ['11','12']

  if season == 'MJJAS':
    starts_ukmo = ['0201','0209','0217','0225',
                   '0301','0309','0317','0325'
                   '0401','0409','0417','0425',
                   '0501','0509','0517','0525',
                   '0601','0609','0617','0625',
                   '0701','0709','0717','0725',
                   '0801','0809','0817','0825',
                   '0901','0909','0917','0925']
    starts_mon = ['02','03','04','05','06','07','08','09']
    starts_day = ['01','09','17','25']
    keep_mon = [5,6,7,8,9]
    starts_mon_bam = ['0']

  # create empty arrays for the data
  week_mid_dates = np.nan*np.zeros((len(years),nleads,len(starts_mon)*len(starts_day),3))
  week_mid_dates_bam = np.nan*np.zeros((len(years),nleads,len(starts_mon_bam),2,3))

  all_chirps = np.nan*np.zeros((len(years),nleads,len(starts_mon)*len(starts_day),len(all_lat),len(all_lon)))
  all_ukmo_ens = np.nan*np.zeros((len(years),nleads,len(starts_mon)*len(starts_day),ukmo_nmembers,len(all_lat),len(all_lon)))
  all_ncep_ens = np.nan*np.zeros((len(years),nleads,len(starts_mon)*len(starts_day),ncep_lags,ncep_nmembers,len(all_lat),len(all_lon)))
  all_ecmf_ens = np.nan*np.zeros((len(years),nleads,len(starts_mon)*len(starts_day),ecmf_lags,ecmf_nmembers,len(all_lat),len(all_lon)))
  all_bam_ens = np.nan*np.zeros((len(years),nleads,len(starts_mon_bam),2,bam_nmembers,len(all_lat),len(all_lon)))
  all_chirps_bam = np.nan*np.zeros((len(years),nleads,len(starts_mon_bam),2,len(all_lat),len(all_lon)))

  #Read data for one season and for all years
  yc = 0
  for y in years:
    print y
    yc = yc + 1
    # create a mask to mask out all dates not within DJF
    ukmo_file	= []
    ncep_file	= []
    ecmf_file	= []
    chirps_file	= []
    bam_file = []
    chirps_file_bam = []

    i = 0
    for mo in np.arange(0,len(starts_mon)):
      if season == 'NDJFM' and starts_mon[mo] in ['11','12']:
        y_ec	= y-1
      else:
        y_ec	= np.copy(y)

      ukmo_file.append(dir_in+region+'_UKMO_ensemble_weekly_hindcasts_lead1-5_'+str(y_ec)+starts_mon[mo]+'.nc')
      print ukmo_file[mo]+' '+str(os.path.exists(ukmo_file[mo]))
      ncep_file.append(dir_in+region+'_NCEP_ensemble_weekly_hindcasts_lead1-5_'+str(y_ec)+starts_mon[mo]+'.nc')
      print ncep_file[mo]+''+str(os.path.exists(ncep_file[mo]))
      ecmf_file.append(dir_in+region+'_ECMWF_ensemble_weekly_hindcasts_lead1-5_'+str(y_ec)+starts_mon[mo]+'.nc')
      print ecmf_file[mo]+' '+str(os.path.exists(ecmf_file[mo]))
      chirps_file.append(dir_in+region+'_CHIRPS_weekly_'+str(y_ec)+starts_mon[mo]+'.nc')
      print chirps_file[mo]+' '+str(os.path.exists(chirps_file[mo]))

      for dd in starts_day:
        # Get start date as date object
        curr_ec_datetime = datetime.strptime(str(y_ec)+starts_mon[mo]+dd,'%Y%m%d')
        print curr_ec_datetime
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
        # curr_mid = str(int(week_mid_dates[yc-1,l,i,0])).zfill(2)+'/'+str(int(week_mid_dates[yc-1,l,i,1])).zfill(2)+'/'+str(int(week_mid_dates[yc-1,l,i,2]))
        # print 'lead time week '+str(l+1)+' '+curr_mid
        i = i + 1

    # chirps
    nc_fid = MFDataset(chirps_file)
    all_chirps[yc-1,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_precip'][:]),0,1)
    nc_fid.close()

    # UKMO Data, units is kg/m2/day, which needs converted to mm/d. However, nothing need to be done.
    # UKMO has 7 ensemble members
    nc_fid = MFDataset(ukmo_file)
    all_ukmo_ens[yc-1,:,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_precip'][:]),0,1)
    nc_fid.close()
    # NCEP units don't need converted. NCEP has 7 lags and 4 ensemble members
    nc_fid = MFDataset(ncep_file)
    all_ncep_ens[yc-1,:,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_precip'][:]),0,1)
    nc_fid.close()
    # ECMWF units don't need converted. There are 3 lags and 11 ensemble members
    nc_fid = MFDataset(ecmf_file)
    all_ecmf_ens[yc-1,:,:,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_precip'][:]),0,1)
    nc_fid.close()


    # Open BAM data separately to UKMO/NCEP/ECMWF because the start dates are different
    if season == 'NDJFM':
      for mo in np.arange(0,len(starts_mon_bam)):
        if season == 'NDJFM' and starts_mon_bam[mo] in ['11','12']:
          y_ec	= y-1
        else:
          y_ec	= np.copy(y)

        curr_bam_file = dir_in+region+'_BAM_ensemble_weekly_hindcasts_lead1-5_'+str(y_ec)+starts_mon_bam[mo]+'.nc'
        nc_fid = Dataset(curr_bam_file, 'r')
        tmp_bam_dates = np.array(nc_fid.variables['week'][:])
        t_units = str(nc_fid.variables['week'].units)
        all_bam_ens[yc-1,:,mo,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_precip'][:]),0,1)
        nc_fid.close()
        curr_bam_dates = nc4.num2date(tmp_bam_dates,t_units)

        curr_chirps_file = dir_in+region+'_CHIRPS_weekly_for_BAM_'+str(y_ec)+starts_mon_bam[mo]+'.nc'
        nc_fid = Dataset(curr_chirps_file, 'r')
        all_chirps_bam[yc-1,:,mo,:,:,:] = np.swapaxes(np.array(nc_fid.variables['week_precip'][:]),0,1)
        nc_fid.close()

        for dd in [0,1]: # only two bam starts per month
          # Get start date as date object
          for l in np.arange(0,nleads):
            # get dates of the 4th day of each week (this keeps samples consistent at each lead time)
            # Note:
            #  1) if you use the date of the first day in the week this excludes forecasts
            #     with start-dates e.g. 30th November, causing different samples at each lead time.
            #  2) if you use the date of the last day in the week this excludes forecasts where
            #     the end of the week stretches over into the next month.
            week_mid_dates_bam[yc-1,l,mo,dd,0] = (curr_bam_dates[dd,l]+timedelta(days=3)).day
            week_mid_dates_bam[yc-1,l,mo,dd,1] = (curr_bam_dates[dd,l]+timedelta(days=3)).month
            week_mid_dates_bam[yc-1,l,mo,dd,2] = (curr_bam_dates[dd,l]+timedelta(days=3)).year
          # curr_mid = str(int(week_mid_dates[yc-1,l,i,0])).zfill(2)+'/'+str(int(week_mid_dates[yc-1,l,i,1])).zfill(2)+'/'+str(int(week_mid_dates[yc-1,l,i,2]))
          # print 'lead time week '+str(l+1)+' '+curr_mid


  '''
  Mask data for a season, compute lag/ensemble means, and create a dry mask
  '''
  # set negatives to 0
  # for DJF, any other month is masked
  # the size is (440, 3), becasuse date_mask is tuple, therefore, it needs 3 dimensions to identify one position
  # the total is 1,100, 8 out of 20 months have been masked out.
  if season == 'NDJFM':
    date_mask = np.where((week_mid_dates[:,:,:,1] > 5) & (week_mid_dates[:,:,:,1] < 11))
    date_mask_bam = np.where((week_mid_dates_bam[:,:,:,:,1] > 5) & (week_mid_dates_bam[:,:,:,:,1] < 11))
    all_bam_ens[date_mask_bam[0],date_mask_bam[1],date_mask_bam[2],date_mask_bam[3],:,:,:] = np.nan
    all_chirps_bam[date_mask_bam[0],date_mask_bam[1],date_mask_bam[2],date_mask_bam[3],:,:] = np.nan
  if season == 'MJJAS':
    date_mask = np.where((week_mid_dates[:,:,:,1] > 9) & (week_mid_dates[:,:,:,1] < 5))
    date_mask_bam = []

  all_ukmo_ens[date_mask[0],date_mask[1],date_mask[2],:,:,:] = np.nan
  all_ncep_ens[date_mask[0],date_mask[1],date_mask[2],:,:,:,:] = np.nan
  all_ecmf_ens[date_mask[0],date_mask[1],date_mask[2],:,:,:,:] = np.nan
  all_chirps[date_mask[0],date_mask[1],date_mask[2],:,:] = np.nan

  # test the mask
  for y in np.arange(0,len(years)):
    for l in np.arange(0,5):
      for i in np.arange(0,20):
        curr_mid = str(int(week_mid_dates[y,l,i,0])).zfill(2)+'/'+str(int(week_mid_dates[y,l,i,1])).zfill(2)+'/'+str(int(week_mid_dates[y,l,i,2]))
        print str(years[y])+' lead week '+str(l)+' '+starts_ukmo[i]+' '+curr_mid+' '+str(all_chirps[y,l,i,20,23])
  if season == 'NDJFM':
    for y in np.arange(0,len(years)):
      for l in np.arange(0,5):
        for i in np.arange(0,5):
          for s in [0,1]:
            curr_mid = str(int(week_mid_dates_bam[y,l,i,s,0])).zfill(2)+'/'+str(int(week_mid_dates_bam[y,l,i,s,1])).zfill(2)+'/'+str(int(week_mid_dates_bam[y,l,i,s,2]))
            print str(years[y])+' lead week '+str(l)+' '+' '+curr_mid+' '+str(all_chirps_bam[y,l,i,s,20,23])
  return (all_chirps,all_ukmo_ens,all_ncep_ens,all_ecmf_ens,all_chirps_bam,all_bam_ens,week_mid_dates,week_mid_dates_bam)

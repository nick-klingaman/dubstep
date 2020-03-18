'''
Save S2S precipitation forecasts over a given region
from UKMO (initialised weekly) and NCEP
(initialised daily), w.r.t weekly forecasts
of UKMO (intialised once per week)

Saves lagged ensembles of NCEP forecasts
i.e. saves the 7 NCEP forecasts intialised
on and before each UKMO initialisation date

UKMO data are 1 degree resolution
NCEP data are 1.5 degree resolution

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
from mpl_toolkits import basemap

execfile('date_str.py')
execfile('grab_data.py')

# input directories
dir_ukmo = '/group_workspaces/jasmin2/klingaman/datasets/S2S/UKMO/ENSO/'
dir_ncep = '/group_workspaces/jasmin2/klingaman/datasets/S2S/NCEP/ENSO/'
dir_ecmf = '/group_workspaces/jasmin2/klingaman/datasets/S2S/ECMWF/ENSO/'
dir_gpcp = '/group_workspaces/jasmin2/klingaman/datasets/GPCP/one_degree/native_resolution/'

# output directories to save forecast data
dir_out = '/group_workspaces/jasmin2/ncas_climate/users/myoung02/datasets/DUBSTEP/'

nleads = 5 # number of lead times (weeks)
nmembers_ukmo = 7 # number of UKMO ensemble members
nmembers_ncep = 4 # number of NCEP ensemble members
nlags_ncep = 7 # number of lags for NCEP - for lagged ensemble means
# Years over which to
# Years available, UKMO 1996-2017, NCEP 1999-2010, ECMWF 1998-2016
years = np.arange(2000,2010+1,1) #years to read and save

# Define region for analysis over Brazil
region = 'Brazil'
latlim = [-40,20]
lonlim = [-90,-20]
ukmo_lon = np.arange(lonlim[0],lonlim[1]+1,1)
ukmo_lat = np.arange(latlim[0],latlim[1]+1,1)
lons_ukmo_rg,lats_ukmo_rg = np.meshgrid(ukmo_lon,ukmo_lat)

# list start dates for ukmo (just for DJF)
starts_ukmo = ['1025','1101','1109','1117','1125','1201',
               '1209','1217','1225','0101','0109','0117',
               '0125','0201','0209','0217','0225']
# Use version 2017 for 25 March to 25 December and version 18 for 1Jan-17Mar
version_list1 = ['1025','1101','1109','1117','1125','1201','1209','1217','1225']
version_list2 = ['0101','0109','0117','0125','0201','0209','0217']
# start dates for ecmwf
starts_ecmf =np.array(['1016','1019','1023','1026','1030','1102','1106','1109','1113',
              '1116','1120','1123','1127','1130','1204','1207','1211',
              '1214','1218','1221','1225','1228','0101','0104','0108',
              '0111','0115','0118','0122','0125','0129','0205','0208',
              '0212','0215','0219','0222','0226','0301'])


yc = 0
for y in years:
  yc = yc + 1
  print y
  # loop through UKMO start dates
  for s in np.arange(0,len(starts_ukmo)):
    print str(y)+'-'+starts_ukmo[s]
    curr_start = starts_ukmo[s] # get current start date
    curr_d = int(curr_start[2:4]) # day from start date
    curr_m = int(curr_start[0:2]) # month from start date

    # Version of UKMO to use
    if curr_start in version_list1:
      v_year = '2017'
    elif curr_start in version_list2:
      v_year = '2018'

    # year to use (Y-1 for July-December)
    if curr_m in [7,8,9,10,11,12]:
      y_use = y-1
    else:
      y_use = np.copy(y)

    # create date time array for current date
    curr_datetime =datetime.strptime(str(y_use)+curr_start,'%Y%m%d')
    curr_dir_ukmo = dir_ukmo+'hindcasts_'+curr_start+'_'+v_year+'/'
    ukmo_fstamp = 'ukmo_enfh_atmos_day_sfc_GloSea5-GC2_'+str(y_use)+curr_start+'00*pr*.nc'
    # Loop through the UKMO ensemble members for year and lead time
    f_ls = glob.glob(curr_dir_ukmo+ukmo_fstamp)  #-2015080500_pr_r201.nc
    for f in np.arange(0,len(f_ls)):
      (ukmo_rfe,tmp_ukmo_lon,tmp_ukmo_lat,ukmo_dates) = grab_ukmo_region(f_ls[f],lonlim,latlim)
      # create an empty array to store all the data
      # if (yt == years[0]) & (ltime == 0) & (week == 0) & (f == 0):
      #   all_ukmo = np.zeros((5,len(years),12,len(ukmo_lat),len(ukmo_lon)))
      if (s == 0) & (f == 0):
        all_ukmo = np.zeros((nleads,len(starts_ukmo),nmembers_ukmo,len(ukmo_lat),len(ukmo_lon)))

      # loop through the lead times
      for lead in np.arange(0,5):
        if lead == 0: # week 1
          week_id = np.arange(1,7+1)
        elif lead == 1: # week 2
          week_id = np.arange(8,14+1)
        elif lead == 2: # week 3
          week_id = np.arange(15,21+1)
        elif lead == 3: # week 4
          week_id = np.arange(22,28+1)
        elif lead == 4: # week 5
          week_id = np.arange(29,35+1)

        curr_dates = ukmo_dates[week_id-1,:]
        sub_ukmo_rfe = ukmo_rfe[week_id-1,:,:]

        # if grid res for UKMO < 1 degree, regrid
        if abs(tmp_ukmo_lon[0]-tmp_ukmo_lon[1]) < 1:
          # regrid UKMO data
          for gg in np.arange(0,7):
            tmp_ukmo_rfe = np.zeros((nmembers_ukmo,len(ukmo_lat),len(ukmo_lon)))
            tmp_ukmo_rfe[gg,:,:] = basemap.interp(sub_ukmo_rfe[gg,:,:],tmp_ukmo_lon,tmp_ukmo_lat,lons_ukmo_rg,lats_ukmo_rg,order=1)
          sub_ukmo_rfe = []
          sub_ukmo_rfe = np.copy(tmp_ukmo_rfe)
        # print curr_dates
        # compute average of ukmo data for the current week
        all_ukmo[lead,s,f,:,:] = np.mean(sub_ukmo_rfe,axis=0)

        if f == 0: # just do for first UKMO file
          if lead == 0: # only open the ncep,ecmwf and gpcp if lead == 0
            '''
            Open NCEP & compute -6d lagged ensemble mean
            '''
            print 'opening NCEP'
            # get lagged dates from day 0 to 6 days before
            for l in np.arange(0,nlags_ncep):
              print 'reading lead '+str(l+1)+' of '+str(nlags_ncep)
              # subtract 'l' days from current date to get start
              # date to open lagged NCEP initialised forecasts
              lag_date = curr_datetime - timedelta(days=l)
              str_lag_start = lag_date.strftime('%m%d')
              str_lag_year = lag_date.strftime('%Y')

              curr_dir_ncep = dir_ncep+'hindcasts_'+str_lag_start+'_6hr/'
              ncep_fstamp = curr_dir_ncep+'ncep_precip_ens_'+str_lag_year+str_lag_start+'.nc'
              date = lag_date.strftime('%Y-%m-%d')
              (ncep_rfe,ncep_lon,ncep_lat,ncep_dates) = grab_ncep_region(date,ncep_fstamp,lonlim,latlim)

              # create an empty array to store all the data
              if (s == 0) & (lead == 0) & (l == 0):
                all_ncep = np.zeros((nleads,len(starts_ukmo),nlags_ncep,nmembers_ncep,len(ncep_lat),len(ncep_lon)))
                lons_rg, lats_rg = np.meshgrid(ncep_lon,ncep_lat)

              if l == 0:
                tmp_lens = np.zeros((nlags_ncep,nmembers_ncep,len(ncep_dates),len(ncep_lat),len(ncep_lon)))
                tmp_dates = np.zeros((nlags_ncep,len(ncep_dates),4))
              # store lags and lagged dates
              tmp_lens[l,:,:,:,:] = np.copy(ncep_rfe)
              tmp_dates[l,:,:] = np.copy(ncep_dates)

          '''
          Find lagged ensembles for the current week
          '''
          # find lagged ensemble for current lead time date NCEP
          for l in np.arange(0,nlags_ncep):
            ncep_y1_id = np.where(tmp_dates[l,:,0] == curr_dates[0,0])[0]
            ncep_y2_id = np.where(tmp_dates[l,:,0] == curr_dates[-1,0])[0]
            ncep_m1_id = np.where(tmp_dates[l,:,1] == curr_dates[0,1])[0]
            ncep_m2_id = np.where(tmp_dates[l,:,1] == curr_dates[-1,1])[0]
            ncep_d1_id = np.where(tmp_dates[l,:,2] == curr_dates[0,2])[0]
            ncep_d2_id = np.where(tmp_dates[l,:,2] == curr_dates[-1,2])[0]
             # always start accumlation from hour 6 (assuming this is the accumulation between 0z-6z)
            ncep_h1_id = np.where(tmp_dates[l,:,3] == 6)[0]
            ncep_ym1_id = np.intersect1d(ncep_y1_id,ncep_m1_id)
            ncep_ym2_id = np.intersect1d(ncep_y2_id,ncep_m2_id)
            ncep_ymd1_id = np.intersect1d(ncep_ym1_id,ncep_d1_id)
            ncep_ymdh1_id = np.intersect1d(ncep_ymd1_id,ncep_h1_id)
            ncep_ymd2_id = np.intersect1d(ncep_ym2_id,ncep_d2_id)
            # have to add an extra index onto the end of time id to
            # rainfall accumulated upto zero Z
            ncep_time_id = np.arange(ncep_ymdh1_id[0],ncep_ymd2_id[len(ncep_ymd2_id)-1]+2)

            # store lagged ensemble
            all_ncep[lead,s,l,:,:,:] = np.sum(tmp_lens[l,:,ncep_time_id,:,:],axis=0)/7

          '''
          Open GPCP
          '''
          for d in np.arange(0,len(week_id)):
            day_str = day_string(int(curr_dates[d,2]))
            mon_str = mon_string(int(curr_dates[d,1]))
            yr_str = str(int(curr_dates[d,0]))

            # open GPCP data for each day in the week
            curr_gpcp = dir_gpcp+yr_str+'/gpcp_v01r03_daily_d'+yr_str+mon_str+day_str+'_c20170530.nc'
            (gpcp_p,gpcp_lon,gpcp_lat) = grab_gpcp_region(curr_gpcp,lonlim,latlim)
            gpcp_p[gpcp_p <0] = 0

            # regrid to ncep grid
            gpcp_rg = []
            gpcp_rg = basemap.interp(gpcp_p,gpcp_lon,gpcp_lat,lons_rg,lats_rg,order=1)

            if (s == 0) & (d == 0):
              gpcp4ukmo = np.zeros((nleads,len(starts_ukmo),len(gpcp_lat),len(gpcp_lon)))
              gpcp4ncep = np.zeros((nleads,len(starts_ukmo),len(ncep_lat),len(ncep_lon)))
            gpcp4ukmo[lead,s,:,:] = gpcp4ukmo[lead,s,:,:] + gpcp_p
            gpcp4ncep[lead,s,:,:] = gpcp4ncep[lead,s,:,:] + gpcp_rg

  # At end of each year
  # save all data as netcdf
  nlats = len(ukmo_lat)
  nlons = len(ukmo_lon)
  fname_ukmo = dir_out+region+'_UKMO_ensemble_weekly_hindcasts_lead1-5_DJF_'+str(y)+'.nc'
  dataset = Dataset(fname_ukmo,'w',format='NETCDF3_CLASSIC')
  lat = dataset.createDimension('lat', nlats) # create lat (dims depend on region)
  lon = dataset.createDimension('lon', nlons) # create lon
  lead_nc = dataset.createDimension('lead',nleads) # create time
  week_nc = dataset.createDimension('week',len(starts_ukmo)) # create time
  ens_nc = dataset.createDimension('ensemble',nmembers_ukmo)

  # create variables
  latitudes = dataset.createVariable('latitude','f',('lat',))
  longitudes = dataset.createVariable('longitude','f',('lon',))
  weeks = dataset.createVariable('week',np.float64,('week',))
  leads = dataset.createVariable('lead',np.float64,('lead',))
  ensembles = dataset.createVariable('ensemble',np.float64,('ensemble',))
  p_mean = dataset.createVariable('week_precip','d',('lead','week','ensemble','lat','lon'))

  # Global Attributes (will need modified accordingly)
  dataset.description = 'UKMO GloSea5-GC2 ensembles weekly mean precipitation'
  dataset.history = 'Created ' + tt.ctime(tt.time())
  dataset.source = 'Subset by M. Young'
  # Variable Attributes
  latitudes.units = 'degrees_north'
  longitudes.units = 'degrees_east'
  f_week = starts_ukmo[0]
  weeks.units = 'weeks since'+str(y-1)+'-'+f_week[0:2]+'-'+f_week[2:4]
  weeks.calendar = 'gregorian'
  p_mean.units = 'kg m-2 s-1'
  p_mean.long_name = 'Weekly mean precipitation'

  # Fill variables with data
  latitudes[:] = ukmo_lat
  longitudes[:] = ukmo_lon
  p_mean[:] = all_ukmo
  weeks[:] = np.arange(1,len(starts_ukmo)+1)
  leads[:] = np.arange(1,nleads+1)
  ensembles[:] = np.arange(1,nmembers_ukmo+1)
  dataset.close()
  print 'saved UKMO file for '+str(y)

  # save NCEP as netcdf
  nlats = len(ncep_lat)
  nlons = len(ncep_lon)
  fname_ncep = dir_out+region+'_NCEP_ensemble_weekly_hindcasts_lead1-5_DJF_'+str(y)+'.nc'
  dataset = Dataset(fname_ncep,'w',format='NETCDF3_CLASSIC')
  lat = dataset.createDimension('lat', nlats) # create lat (dims depend on region)
  lon = dataset.createDimension('lon', nlons) # create lon
  lead_nc = dataset.createDimension('lead',nleads) # create time
  week_nc = dataset.createDimension('week',len(starts_ukmo)) # create time
  lag_nc = dataset.createDimension('lag',nlags_ncep)
  ens_nc = dataset.createDimension('ensemble',nmembers_ncep)

  # create variables
  latitudes = dataset.createVariable('latitude','f',('lat',))
  longitudes = dataset.createVariable('longitude','f',('lon',))
  leads = dataset.createVariable('lead',np.float64,('lead',))
  weeks = dataset.createVariable('week',np.float64,('week',))
  lags = dataset.createVariable('lag',np.float64,('lag',))
  ensembles = dataset.createVariable('ensemble',np.float64,('ensemble',))
  p_mean = dataset.createVariable('week_precip','d',('lead','week','lag','ensemble','lat','lon'))

  # Global Attributes (will need modified accordingly)
  dataset.description = 'NCEP lagged ensembles weekly mean precipitation'
  dataset.history = 'Created ' + tt.ctime(tt.time())
  dataset.source = 'Subset by M. Young'
  # Variable Attributes
  latitudes.units = 'degrees_north'
  longitudes.units = 'degrees_east'
  f_week = starts_ukmo[0]
  weeks.units = 'weeks since'+str(y-1)+'-'+f_week[0:2]+'-'+f_week[2:4]
  weeks.calendar = 'gregorian'
  p_mean.units = 'kg m-2'
  p_mean.long_name = 'Weekly mean precipitation'

  # Fill variables with data
  latitudes[:] = ncep_lat
  longitudes[:] = ncep_lon
  p_mean[:] = all_ncep
  leads[:] = np.arange(1,nleads+1)
  weeks[:] = np.arange(1,len(starts_ukmo)+1)
  lags[:] = np.arange(0,nlags_ncep)
  ensembles[:] = np.arange(1,nmembers_ncep+1)
  dataset.close()
  print 'saved ncep file for '+str(y)

  # save gpcp data as netcdf
  nlats = len(gpcp_lat)
  nlons = len(gpcp_lon)
  fname_gpcp = dir_out+region+'_GPCP_weekly_DJF_UKMO_times_'+str(y)+'.nc'
  dataset = Dataset(fname_gpcp,'w',format='NETCDF3_CLASSIC')
  lat = dataset.createDimension('lat', nlats) # create lat (dims depend on region)
  lon = dataset.createDimension('lon', nlons) # create lon
  lead_nc = dataset.createDimension('lead',nleads) # create time
  week_nc = dataset.createDimension('week',len(starts_ukmo)) # create time
  # create variables
  latitudes = dataset.createVariable('latitude','f',('lat',))
  longitudes = dataset.createVariable('longitude','f',('lon',))
  weeks = dataset.createVariable('week',np.float64,('week',))
  leads = dataset.createVariable('lead',np.float64,('lead',))
  p_mean = dataset.createVariable('week_precip','d',('lead','week','lat','lon'))

  # Global Attributes (will need modified accordingly)
  dataset.description = 'GPCP weekly mean precipitation at 1 degree res, based on UKMO intialisation dates'
  dataset.history = 'Created ' + tt.ctime(tt.time())
  dataset.source = 'Subset by M. Young'
  # Variable Attributes
  latitudes.units = 'degrees_north'
  longitudes.units = 'degrees_east'
  f_week = starts_ukmo[0]
  weeks.units = 'weeks since'+str(y-1)+'-'+f_week[0:2]+'-'+f_week[2:4]
  weeks.calendar = 'gregorian'
  p_mean.units = 'mm d-1'
  p_mean.long_name = 'Weekly mean precipitation'

  # Fill variables with data
  latitudes[:] = gpcp_lat
  longitudes[:] = gpcp_lon
  p_mean[:] = gpcp4ukmo/7 # convert weekly mean of mm d rather than weekly total
  weeks[:] = np.arange(1,len(starts_ukmo)+1)
  leads[:] = np.arange(1,nleads+1)
  dataset.close()
  print 'saved gpcp file for '+str(y)

  # save all data as netcdf
  nlats = len(ncep_lat)
  nlons = len(ncep_lon)
  fname_gpcp = dir_out+region+'_GPCP_weekly_DJF_NCEP_times_'+str(y)+'.nc'
  dataset = Dataset(fname_gpcp,'w',format='NETCDF3_CLASSIC')
  lat = dataset.createDimension('lat', nlats) # create lat (dims depend on region)
  lon = dataset.createDimension('lon', nlons) # create lon
  lead_nc = dataset.createDimension('lead',nleads) # create time
  week_nc = dataset.createDimension('week',len(starts_ukmo)) # create time
  # create variables
  latitudes = dataset.createVariable('latitude','f',('lat',))
  longitudes = dataset.createVariable('longitude','f',('lon',))
  weeks = dataset.createVariable('week',np.float64,('week',))
  leads = dataset.createVariable('lead',np.float64,('lead',))
  p_mean = dataset.createVariable('week_precip','d',('lead','week','lat','lon'))

  # Global Attributes (will need modified accordingly)
  dataset.description = 'GPCP weekly mean precipitation at NCEP forecast res, based on UKMO intialisation dates'
  dataset.history = 'Created ' + tt.ctime(tt.time())
  dataset.source = 'Subset by M. Young'
  # Variable Attributes
  latitudes.units = 'degrees_north'
  longitudes.units = 'degrees_east'
  f_week = starts_ukmo[0]
  weeks.units = 'weeks since'+str(y-1)+'-'+f_week[0:2]+'-'+f_week[2:4]
  weeks.calendar = 'gregorian'
  p_mean.units = 'mm d-1'
  p_mean.long_name = 'Weekly mean precipitation'

  # Fill variables with data
  latitudes[:] = ncep_lat
  longitudes[:] = ncep_lon
  p_mean[:] = gpcp4ncep/7 # convert weekly mean of mm d rather than weekly total
  weeks[:] = np.arange(1,len(starts_ukmo)+1)
  leads[:] = np.arange(1,nleads+1)
  dataset.close()
  print 'saved gpcp ncep file for '+str(y)

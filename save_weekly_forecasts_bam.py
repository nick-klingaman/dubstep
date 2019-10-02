#!/usr/bin/python2.7
#BSUB -o %J.o
#BSUB -e %J.e
#BSUB -q short-serial
#BSUB -W 24:00
##BSUB -R "rusage[mem=16000]"
##BSUB -M 16000

'''
Save S2S precipitation hindcasts over a given region
from the Brazil Atmospheric Model (BAM)
BAM hindcasts are initialised twice a month.

M. Young, April 2019
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
import os.path
import pdb
from datetime import date

# execfile() is a python built-in function
# it is different from the import statement in that it does not use the module administation - it reads the file unconditionally and does not create a new module.
execfile('date_str.py')
execfile('grab_data.py')

# input directories
# dir_ukmo = '/gws/nopw/j04/klingaman/datasets/S2S/UKMO/ENSO/'
# dir_ncep = '/gws/nopw/j04/klingaman/datasets/S2S/NCEP/ENSO/'
# dir_ecmf = '/gws/nopw/j04/klingaman/datasets/S2S/ECMWF/ENSO/'
# dir_gpcp = '/gws/nopw/j04/klingaman/datasets/GPCP/one_degree/native_resolution/'
dir_bam = '/gws/nopw/j04/ncas_climate_vol1/users/myoung02/datasets/bam_prec/'
dir_chps = '/gws/nopw/j04/ncas_climate_vol2/users/myoung02/datasets/CHIRPS/global_daily/SAmerica_daily_1.5d/'

# output directories to save forecast data
# dir_out = '../output/'
dir_out = '/gws/nopw/j04/ncas_climate_vol1/users/myoung02/datasets/DUBSTEP_data/'

nleads = 5 # number of lead times (weeks)
# nmembers_ukmo = 7 # number of UKMO ensemble members
# nmembers_ncep = 4 # number of NCEP ensemble members
nmembers_bam = 11 # number of BAM ensemble members
# nlags_ncep = 7 # number of lags for NCEP - for lagged ensemble means
# Years over which to
# Years available, UKMO 1996-2017, NCEP 1999-2010, ECMWF 1998-2016
years = np.arange(2010,2011+1,1) #years to read and save

# Define region for analysis over Brazil
region = 'SAmerica' #'Brazil'
# latlim = [-50,20] # [-40,20]
latlim = [-49.5,19.5] # [-40,20]
lonlim = [-90,-20]
#ukmo_lon = np.arange(lonlim[0],lonlim[1]+1,1)
#ukmo_lat = np.arange(latlim[0],latlim[1]+1,1)
#lons_ukmo_rg,lats_ukmo_rg = np.meshgrid(ukmo_lon,ukmo_lat)

# Note, following variables are (string) list
# list start dates for ukmo (just for DJF)
# change months for another seasons
month = ['01','02','03','04','05','06','07','08','09','10','11','12']

yc = 0
for y in years:
  yc = yc + 1

  # Get the BAM data for the two start dates in each month
  for m in np.arange(0,len(month)):
    # list bam files for the month
    bam_fstamp = dir_bam+'bam_precip_ens_'+str(y)+month[m]
    f_ls = glob.glob(bam_fstamp+'*.nc')
    if len(f_ls) > 0:
      weeks_d = np.zeros((2,5)) # for storing the dates of the weeks
      start_date = date(y,int(month[m]),1)
      time_units = 'days since '+str(y)+'-'+month[m]+'-01 00:00:00'
      # time_units = 'days since '+str(y)+'-'+month[m]+'-'+curr_datetime.strftime('%d')+' 00:00:00'
      for s in [0,1]:
        # create date time array for current date
        # datetime.strptime() creates a datetime object from a string representing a date and time and a corresponding format string.
        # curr_datetime = datetime.strptime(str(y)+curr_start,'%Y%m%d')
        curr_datetime = datetime.strptime(f_ls[s][80:88],'%Y%m%d')

        print '------------'
        print 'opening BAM'
        print '------------'

        print f_ls[s]+"	"+str(os.path.exists(f_ls[s]))
        (bam_rfe,bam_lon,bam_lat,bam_dates) = grab_bam_region(f_ls[s],lonlim,latlim)
  #     print ukmo_rfe.shape			# (7,239,40,47), units: kg/m2 per 6-hour?
  #     print ukmo_lon
  #     print ukmo_lat
  #     print ukmo_dates.shape			# (239,3)
        # create an empty array to store all the data
        if (s == 0):
          all_bam = np.zeros((2,nleads,nmembers_bam,len(bam_lat),len(bam_lon)))# (5, 4, 7, 40, 47)

        # loop through the lead times, 6-hourly data
        for lead in np.arange(0,5):
          if lead == 0:
            week_id = np.arange(1,28)
          elif lead == 1:
            week_id = np.arange(29,56)
          elif lead == 2:
            week_id = np.arange(57,84)
          elif lead == 3:
            week_id = np.arange(85,112)
          elif lead == 4:
            week_id = np.arange(113,140)

          curr_date = date(int(bam_dates[week_id[0]-1,0]),int(bam_dates[week_id[0]-1,1]),int(bam_dates[week_id[0]-1,2]))
          weeks_d[s,lead] = (curr_date - start_date).days

          curr_dates = bam_dates[week_id-1,:]
          sub_bam_rfe = bam_rfe[:,week_id-1,:,:]

          # compute average of ukmo data for the current week, kg/m2 per day?
          # multiply by 4 to convert mean in mm/6hr to mm/d
          all_bam[s,lead,:,:,:] = np.mean(sub_bam_rfe,axis=1)*4.

          print '------------'
          print 'opening CHIRPS'
          print '------------'
          for d in np.arange(0,len(week_id),4):
            day_str = str(int(curr_dates[d,2])).zfill(2)
            mon_str = str(int(curr_dates[d,1])).zfill(2)
            yr_str = str(int(curr_dates[d,0]))
            '''
            # open GPCP data for each day in the week
            curr_gpcp = dir_gpcp+yr_str+'/gpcp_v01r03_daily_d'+yr_str+mon_str+day_str+'_c20170530.nc'
            print curr_gpcp+'	'+str(os.path.exists(curr_gpcp))+'	'+str(lead)+'	'+str(s)+'	'+str(d)
            (gpcp_p,gpcp_lon,gpcp_lat) = grab_gpcp_region(curr_gpcp,lonlim,latlim)
            gpcp_p[gpcp_p <0] = 0

            # regrid to ncep grid
            gpcp_rg = []	# it is a list, but what does this line do?
            # interpolate data on a rectilinear grid (with x = gpcp_lon, y = gpcp_lat) to a grid with x =lons_rg, y = lats_rg.
            # gpcp_p is a rank-2 array with 1st dimension corresponding to y, 2nd dimension x.
            # gpcp_lon, gpcp_lat are rank-1 arrays containing x and y of gpcp_p grid in INCREASING order.
            # lons_rg, lats_rg are rand-2 arrays containing x and y of desired output grid.
            gpcp_rg = basemap.interp(gpcp_p,gpcp_lon,gpcp_lat,lons_rg,lats_rg[::-1,:],order=1) # (40, 47)
            '''
            # Open CHIRPS data for each day in the week
            curr_chirps = dir_chps+'CHIRPSv2_SAmerica_1.5d_'+yr_str+mon_str+day_str+'.nc'
            print curr_chirps+'	'+str(os.path.exists(curr_chirps))+'	'+str(lead)+'	'+str(s)+'	'+str(d)

            # (gpcp_p,gpcp_lon,gpcp_lat) = grab_gpcp_region(curr_gpcp,lonlim,latlim)
            # gpcp_p[gpcp_p <0] = 0
            chirps_p,chirps_lon,chirps_lat = grab_chirps(curr_chirps)
            if (lead == 0) & (s == 0) & (d == 0):
              # gpcp4ncep = np.zeros((nleads,len(starts_ukmo),len(ncep_lat),len(ncep_lon))) # (5, 4, 40, 47)
              chirps_weekly = np.zeros((2,nleads,len(bam_lat),len(bam_lon)))

            # gpcp4ncep[lead,s,:,:] = gpcp4ncep[lead,s,:,:] + gpcp_rg[::-1,:] # units: mm/week, turn latitude to [20,-40]
            chirps_weekly[s,lead,:,:] = chirps_weekly[s,lead,:,:] + chirps_p[0,:,:] # units: mm/week, turn latitude to [20,-40]
      # save as netcdf
      nlats	= len(bam_lat)
      nlons	= len(bam_lon)
      fname_bam = dir_out+region+'_BAM_ensemble_weekly_hindcasts_lead1-5_'+str(y)+month[m]+'.nc'
      dataset	= Dataset(fname_bam,'w',format='NETCDF3_CLASSIC')
      lat	= dataset.createDimension('lat',nlats)
      lon	= dataset.createDimension('lon',nlons)
      lead_nc	= dataset.createDimension('lead',nleads)
      week_nc	= dataset.createDimension('week',2)
      ens_nc	= dataset.createDimension('ensemble',nmembers_bam)

      # create variables
      latitudes	= dataset.createVariable('latitude','f',('lat',))
      longitudes= dataset.createVariable('longitude','f',('lon',))
      leads	= dataset.createVariable('lead',np.float64,('lead',))
      weeks	= dataset.createVariable('week',np.float64,('week','lead'))
      ensembles	= dataset.createVariable('ensemble',np.float64,('ensemble',))
      p_mean	= dataset.createVariable('week_precip','d',('week','lead','ensemble','lat','lon'))

      # Global Attributes (will need modified accordingly)
      dataset.description	= 'BAM weekly mean precipitation'
      dataset.history	= 'Created ' + tt.ctime(tt.time())
      dataset.source	= 'Subset by M. Young'
      # Variable Attributes
      latitudes.units	= 'degrees_north'
      longitudes.units = 'degrees_east'
      weeks.units = time_units
      weeks.calendar = 'none'
      p_mean.units = 'kg m-2 / day'
      p_mean.long_name = 'Weekly mean precipitation'

      # Fill variables with data
      latitudes[:]	= bam_lat
      longitudes[:]	= bam_lon
      p_mean[:] = all_bam
      leads[:] = np.arange(1,nleads+1)
      weeks[:] = weeks_d
      ensembles[:] = np.arange(1,nmembers_bam+1)
      dataset.close()
      print ''
      print 'saved ncep file for '+str(y)+'-'+month[m]


      # save GPCP  as netcdf
      nlats	= len(bam_lat)
      nlons	= len(bam_lon)
      fname_gpcp	= dir_out+region+'_CHIRPS_weekly_for_BAM_'+str(y)+month[m]+'.nc'
      dataset	= Dataset(fname_gpcp,'w',format='NETCDF3_CLASSIC')
      lat = dataset.createDimension('lat',nlats)
      lon = dataset.createDimension('lon',nlons)
      lead_nc	= dataset.createDimension('lead',nleads)
      week_nc	= dataset.createDimension('week',None)

      # create variables
      latitudes	= dataset.createVariable('latitude','f',('lat',))
      longitudes	= dataset.createVariable('longitude','f',('lon',))
      leads	= dataset.createVariable('lead',np.float64,('lead',))
      weeks	= dataset.createVariable('week',np.float64,('week','lead'))
      p_mean	= dataset.createVariable('week_precip','d',('week','lead','lat','lon'))

      # Global Attributes (will need modified accordingly)
      dataset.description	= 'CHIRPS weekly mean precipitation at forecast 1.5 degrees, based on BAM intialisation dates'
      dataset.history	= 'Created ' + tt.ctime(tt.time())
      dataset.source	= 'Subset by M. Young'
      # Variable Attributes
      latitudes.units	= 'degrees_north'
      longitudes.units	= 'degrees_east'
      weeks.units		= time_units
      weeks.calendar	= 'none'
      p_mean.units	= 'mm d-1'
      p_mean.long_name	= 'Weekly mean precipitation'

      # Fill variables with data
      latitudes[:] = chirps_lat
      longitudes[:]	= chirps_lon
      p_mean[:]		= chirps_weekly/7.			# convert weekly mean of mm d rather than weekly total
      leads[:] = np.arange(1,nleads+1)
      weeks[:] = weeks_d
      dataset.close()
      print ''
      print 'saved chirps file for '+str(y)+'-'+month[m]
    else:
      print 'no files for '+str(y)+'-'+month[m]

# # save gpcp data as netcdf
# nlats = len(gpcp_lat)
# nlons = len(gpcp_lon)
# fname_gpcp = dir_out+region+'_GPCP_weekly_DJF_UKMO_times_'+str(y)+'.nc'
# dataset = Dataset(fname_gpcp,'w',format='NETCDF3_CLASSIC')
# lat = dataset.createDimension('lat', nlats) # create lat (dims depend on region)
# lon = dataset.createDimension('lon', nlons) # create lon
# lead_nc = dataset.createDimension('lead',nleads) # create time
# week_nc = dataset.createDimension('week',len(starts_ukmo)) # create time
# # create variables
# latitudes = dataset.createVariable('latitude','f',('lat',))
# longitudes = dataset.createVariable('longitude','f',('lon',))
# weeks = dataset.createVariable('week',np.float64,('week',))
# leads = dataset.createVariable('lead',np.float64,('lead',))
# p_mean = dataset.createVariable('week_precip','d',('lead','week','lat','lon'))

# # Global Attributes (will need modified accordingly)
# dataset.description = 'GPCP weekly mean precipitation at 1 degree res, based on UKMO intialisation dates'
# dataset.history = 'Created ' + tt.ctime(tt.time())
# dataset.source = 'Subset by M. Young'
# # Variable Attributes
# latitudes.units = 'degrees_north'
# longitudes.units = 'degrees_east'
# f_week = starts_ukmo[0]
# weeks.units = 'weeks since'+str(y-1)+'-'+f_week[0:2]+'-'+f_week[2:4]
# weeks.calendar = 'gregorian'
# p_mean.units = 'mm d-1'
# p_mean.long_name = 'Weekly mean precipitation'

# # Fill variables with data
# latitudes[:] = gpcp_lat
# longitudes[:] = gpcp_lon
# p_mean[:] = gpcp4ukmo/7 # convert weekly mean of mm d rather than weekly total
# weeks[:] = np.arange(1,len(starts_ukmo)+1)
# leads[:] = np.arange(1,nleads+1)
# dataset.close()
# print 'saved gpcp file for '+str(y)

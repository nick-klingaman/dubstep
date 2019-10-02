'''
Script to subset CHIRPS into daily files for S2S analysis over Brazil
M. Young, April 2019
'''
from netCDF4 import Dataset
import netCDF4 as nc4
import datetime as dt
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import basemap
import os
import time as tt
from datetime import date
execfile('date_str.py')

# dir_chp = '/gws/nopw/j04/ncas_climate_vol2/users/myoung02/datasets/CHIRPS/global_monthly/'
dir_chp = '/gws/nopw/j04/ncas_climate_vol2/users/myoung02/datasets/CHIRPS/global_daily/'
dir_out = '/gws/nopw/j04/ncas_climate_vol2/users/myoung02/datasets/CHIRPS/global_daily/SAmerica_daily_1.5d/'


# Define region for analysis over Brazil
region = 'SAmerica' # Brazil
latlim = [-49.5,19.5] # [-40,20]
lonlim = [-90,-20]
ukmo_lon = np.arange(lonlim[0],lonlim[1],1.5)
ukmo_lat = np.arange(latlim[0],latlim[1]+1.5,1.5)
lon_rg,lat_rg = np.meshgrid(ukmo_lon,ukmo_lat)

years = np.arange(2012,2018+1,1)#[2011] #np.arange(1999,2011+1,1)
for y in np.arange(0,len(years)):
  print years[y]
  # open chirps data and regrid
  chp_file = dir_chp+'chirps-v2.0.'+str(years[y])+'.days_p05_'+region+'.nc'
  nc_fid = Dataset(chp_file,'r')
  chp_lat = np.array(nc_fid.variables['latitude'][:])
  chp_lon = np.array(nc_fid.variables['longitude'][:])
  t_time = nc_fid.variables['time'][:]
  t_units = str(nc_fid.variables['time'].units)
  all_dates = nc4.num2date(t_time,t_units)
  chp_yrs = []
  chp_mons = []
  chp_days = []
  for i in range(0,len(all_dates)):
    curr_day = all_dates[i]
    start_date = date(curr_day.year,curr_day.month,curr_day.day)
    time_units = 'days since '+str(curr_day.year)+'-'+mon_string(curr_day.month)+'-'+day_string(curr_day.day)+' 00:00:00'
    date_stamp = str(curr_day.year)+mon_string(curr_day.month)+day_string(curr_day.day)
    print date_stamp
    dayn = (start_date - start_date).days
    chp_rfe = []
    chp_rfe = np.array(nc_fid.variables['precip'][i,:,:],dtype='f').squeeze()
    chp_rfe[chp_rfe < 0] = np.nan
    chp_rg = []
    chp_rg = basemap.interp(chp_rfe,chp_lon,chp_lat,lon_rg,lat_rg,order=1)

    # save re-grids as netcdfs
    nc_outfile = dir_out+'CHIRPSv2_'+region+'_1.5d_'+date_stamp+'.nc'
    dataset = Dataset(nc_outfile,'w',format='NETCDF4')
    time = dataset.createDimension('time', 1) # create time
    lat = dataset.createDimension('lat', len(ukmo_lat)) # create lat (dims depend on region)
    lon = dataset.createDimension('lon', len(ukmo_lon)) # create lon
    # create variables
    rfe_out = dataset.createVariable('rfe','d',('time','lat','lon'))
    latitudes = dataset.createVariable('latitude','f',('lat',))
    longitudes = dataset.createVariable('longitude','f',('lon',))
    times = dataset.createVariable('time', np.float64, ('time',))
    # Global Attributes (will need modified accordingly)
    dataset.description = 'CHIRPS v2.0 daily rainfall, regridded to 1.5 degree.'
    dataset.history = 'Created ' + tt.ctime(tt.time())
    dataset.source = 'Subset by M. Young'
    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    rfe_out.units = 'mm'
    times.units = time_units
    times.calendar = 'standard'
    # Fill variables with data
    latitudes[:] = ukmo_lat
    longitudes[:] = ukmo_lon
    rfe_out[:] = chp_rg
    times[:] = 0
    dataset.close()
  nc_fid.close()

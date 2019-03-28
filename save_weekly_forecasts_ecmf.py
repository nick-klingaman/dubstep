#!/usr/bin/python2.7
#BSUB -o %J.o
#BSUB -e %J.e
#BSUB -q short-serial
#BSUB -W 24:00
##BSUB -R "rusage[mem=16000]"
##BSUB -M 16000

'''
Save bi-weekly initialised S2S precipitation forecasts
from ECMWF over a region w.r.t weekly forecasts
of UKMO (intialised once per week)

Saves lagged ensembles of ECMWF forecasts
- saves the 3 ECMWF forecasts intialised
before each UKMO initialisation date

ECMWF data are 1.5 degree resolution

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
import os.path
import pdb
execfile('date_str.py')
execfile('grab_data.py')

# input directories
dir_ukmo = '/gws/nopw/j04/klingaman/datasets/S2S/UKMO/ENSO/'
dir_ncep = '/gws/nopw/j04/klingaman/datasets/S2S/NCEP/ENSO/'
dir_ecmf = '/gws/nopw/j04/klingaman/datasets/S2S/ECMWF/ENSO/'
dir_gpcp = '/gws/nopw/j04/klingaman/datasets/GPCP/one_degree/native_resolution/'

# output directories to save forecast data
dir_out = '../output/'

nleads = 5 # number of lead times (weeks)
nmembers_ecmf = 11 # number of ECMWF ensemble members
nlags_ecmf = 3 # number of ECMWF forecasts lags to save
# Years available, UKMO 1996-2017, NCEP 1999-2010, ECMWF 1998-2016
years = np.arange(1999,2010+1,1) # years to read and save

# Define region for analysis over Brazil
region = 'Brazil'
latlim = [-40,20]
lonlim = [-90,-20]

# shorter list of start dates for ukmo (just for DJF)
month = ['01','02','03','04','05','06','07','08','09','10','11','12'] 
starts_ukmo = ['01','09','17','25']
# Use version 2017 for 25 March to 25 December and version 18 for 1Jan-17Mar
version_list1 = ['0801','0809','0817','0825',
                 '0901','0909','0917','0925',
                 '1001','1009','1017','1025',
                 '1101','1109','1117','1125',
                 '1201','1209','1217','1225']
version_list2 = ['0101','0109','0117','0125',
                 '0201','0209','0217','0225',
                 '0301','0309','0317','0325',
                 '0401','0409','0417','0425',
                 '0501','0509','0517','0525',
                 '0601','0609','0617','0625',
                 '0701','0709','0717','0725']

# start dates for ecmwf
starts_ecmf =np.array(['0101','0104','0108','0111','0115','0118','0122','0125','0129',
                       '0205','0208','0212','0215','0219','0222','0226',
                       '0301','0305','0308','0312','0315','0319','0322','0326','0329',
                       '0402','0405','0409','0412','0416','0419','0423','0426','0430',
                       '0501','0504','0508','0511','0515','0518','0522','0525','0529',
                       '0601','0605','0608','0612','0615','0619','0622','0626','0629',
                       '0703','0706','0710','0713','0717','0720','0724','0727','0731',
                       '0803','0807','0810','0814','0817','0821','0824','0828','0831',
                       '0904','0907','0911','0914','0918','0921','0925','0928',
                       '1002','1005','1009','1012','1016','1019','1023','1026','1030',
                       '1102','1106','1109','1113','1116','1120','1123','1127','1130',
                       '1204','1207','1211','1214','1218','1221','1225','1228'])

yc = 0
for y in years:
  yc = yc + 1

  # loop through UKMO start dates
  for m in np.arange(0,len(month)):
    for s in np.arange(0,len(starts_ukmo)):

      print '======================='
      print str(y)+'-'+month[m]+'-'+starts_ukmo[s]
      curr_start = month[m]+starts_ukmo[s]
      curr_m = int(month[m])
      curr_d = int(starts_ukmo[s])

      # Version of UKMO to use
      if curr_start in version_list1:
        v_year = '2017'
      elif curr_start in version_list2:
        v_year = '2018'

      # create date time array for current date
      curr_datetime = datetime.strptime(str(y)+curr_start,'%Y%m%d')

      print '========='
      print 'open ukmo'
      print '========='
      curr_dir_ukmo = dir_ukmo+'hindcasts_'+curr_start+'_'+v_year+'_6hr/'
      ukmo_fstamp = 'ukmo_precip_ens_'+str(y)+str(curr_start)+'.nc'
      f_ls = glob.glob(curr_dir_ukmo+ukmo_fstamp)
      print f_ls[0]+'	'+str(os.path.exists(f_ls[0]))

      nc_fid = Dataset(f_ls[0], 'r')
      t_time = np.array(nc_fid.variables['t'][:])
      t_units = str(nc_fid.variables['t'].units)
      uk_dates = nc4.num2date(t_time,t_units)
      nc_fid.close()

      # split dates into year, month and day columns
      ukmo_dates = np.zeros((len(uk_dates),3))
      for t in np.arange(0,len(uk_dates)):
        curr_time = uk_dates[t]
        ukmo_dates[t,0]	= curr_time.year
        ukmo_dates[t,1]	= curr_time.month
        ukmo_dates[t,2]	= curr_time.day

      # loop through the lead times
      for lead in np.arange(0,nleads):
        if lead == 0:
          week_id = np.arange(  1, 28)
        elif lead == 1:
          week_id = np.arange( 29, 56)
        elif lead == 2:
          week_id = np.arange( 57, 84)
        elif lead == 3:
          week_id = np.arange( 85,112)
        elif lead == 4:
          week_id = np.arange(113,140)

        curr_dates = ukmo_dates[week_id-1,:]

        if lead == 0: # only open ecmwf and gpcp if lead == 0
          print '============='
          print 'opening ecmwf'
          print '============='

          #Find the nearest start date to UKMO
          ec_tdiff = np.zeros((len(starts_ecmf)))
          for e in np.arange(0,len(starts_ecmf)):
            c_st = starts_ecmf[e]
            c_mo = int(c_st[0:2])
            y_ec = np.copy(y)
            curr_ec_datetime = datetime.strptime(str(y_ec)+c_st,'%Y%m%d')
            # find time difference between UKMO and ECMWF startdates
            ec_tdiff[e]	= (curr_ec_datetime - curr_datetime).days	# +: ECMWF ahead, -: ecmwf behind UKMO
          # find ecmf start dates <= 0
          id_neg = np.where(ec_tdiff <= 0)[0]
          # get the last 3 initialisation dates (can also try last 2)
          if len(id_neg) == 1:
            ec_date_id = [102,103,0]
          elif len(id_neg) == 2:
            ec_date_id = [103,0,1]
          else:
            ec_date_id = np.arange(id_neg[-3],id_neg[-1]+1)
          ec_lag_dates = starts_ecmf[ec_date_id]

          # get last 3 initialisations for ECMWF
          for l in np.arange(0,nlags_ecmf):
            print 'reading lag '+str(l+1)+' of '+str(nlags_ecmf)
            # subtract 'l' days from current date to get start
            # date to open lagged NCEP initialised forecasts
            ec_lag_date = ec_lag_dates[l]
            ec_lag_mo = int(ec_lag_date[0:2])
            ec_lag_year = str(y)

            # gets correct ecmf model version
            if ec_lag_mo in np.arange(5,12+1):
              v_ec_year = '2017'
            elif ec_lag_mo in np.arange(1,4+1):
              v_ec_year = '2018'

            if len(id_neg) == 1 and (l == 0 or l == 1):
              ec_lag_year = str(int(ec_lag_year) - 1)
            elif len(id_neg) == 2 and l == 0:
              ec_lag_year = str(int(ec_lag_year) - 1)
            ec_f_date = datetime.strptime(ec_lag_year+ec_lag_date,'%Y%m%d')
            date = ec_f_date.strftime('%Y-%m-%d')

            curr_dir_ecmf = dir_ecmf+'hindcasts_'+ec_lag_date+'_'+v_ec_year+'/'
            ecmf_fstamp = curr_dir_ecmf+'ecmwf_precip_ens_'+ec_lag_year+ec_lag_date+'.nc'
            print ecmf_fstamp+'	'+str(os.path.exists(ecmf_fstamp))
            (ecmf_rfe,ecmf_lon,ecmf_lat,ecmf_dates) = grab_ecmf_region(date,ecmf_fstamp,lonlim,latlim)
	    #print ecmf_rfe.shape			# (11, 188, 40, 47)
            # create an empty array to store all the data
            if (s == 0) & (l == 0):
              all_ecmf = np.zeros((nleads,len(starts_ukmo),nlags_ecmf,nmembers_ecmf,len(ecmf_lat),len(ecmf_lon))) # (5, 4, 3, 11, 40, 47)
              lons_rg_ec, lats_rg_ec = np.meshgrid(ecmf_lon,ecmf_lat)

            if l == 0:
              tmp_ec_lens = np.zeros((nlags_ecmf,nmembers_ecmf,len(ecmf_dates),len(ecmf_lat),len(ecmf_lon))) # (3, 11, 188, 40, 47)
              tmp_ec_dates = np.zeros((nlags_ecmf,len(ecmf_dates),4))
            # store lags and lagged dates
            tmp_ec_lens[l,:,:,:,:] = np.copy(ecmf_rfe)
            tmp_ec_dates[l,:,:] = np.copy(ecmf_dates)

        '''
        Find lagged ensembles for the current week
        '''
        # find lagged ensemble for current lead time date ECMWF
        for l in np.arange(0,nlags_ecmf):
          ecmf_y1_id = np.where(tmp_ec_dates[l,:,0] == curr_dates[ 0,0])[0]
          ecmf_y2_id = np.where(tmp_ec_dates[l,:,0] == curr_dates[-1,0])[0]
          ecmf_m1_id = np.where(tmp_ec_dates[l,:,1] == curr_dates[ 0,1])[0]
          ecmf_m2_id = np.where(tmp_ec_dates[l,:,1] == curr_dates[-1,1])[0]
          ecmf_d1_id = np.where(tmp_ec_dates[l,:,2] == curr_dates[ 0,2])[0]
          ecmf_d2_id = np.where(tmp_ec_dates[l,:,2] == curr_dates[-1,2])[0]
          # always start accumlation from hour 6 (assuming this is the accumulation between 0z-6z)
          ecmf_h1_id = np.where(tmp_ec_dates[l,:,3] == 6)[0]
          ecmf_ym1_id = np.intersect1d(ecmf_y1_id  ,ecmf_m1_id)
          ecmf_ym2_id = np.intersect1d(ecmf_y2_id  ,ecmf_m2_id)
          ecmf_ymd1_id = np.intersect1d(ecmf_ym1_id ,ecmf_d1_id)
          ecmf_ymdh1_id = np.intersect1d(ecmf_ymd1_id,ecmf_h1_id)
          ecmf_ymd2_id = np.intersect1d(ecmf_ym2_id ,ecmf_d2_id)
          # have to add an extra index onto the end of time id to
          # rainfall accumulated upto zero Z
          ecmf_time_id = np.arange(ecmf_ymdh1_id[0],ecmf_ymd2_id[len(ecmf_ymd2_id)-1]+2)		# (28)
	  #print ecmf_time_id
          # store lagged ensemble
          all_ecmf[lead,s,l,:,:,:] = np.sum(tmp_ec_lens[l,:,ecmf_time_id,:,:],axis=0)/7					# kg/m2/day

    # save ECMWF as netcdf
    nlats	= len(ecmf_lat)
    nlons	= len(ecmf_lon)
    fname_ncep	= dir_out+region+'_ECMWF_ensemble_weekly_hindcasts_lead1-5_'+str(y)+month[m]+'.nc'
    dataset	= Dataset(fname_ncep,'w',format='NETCDF3_CLASSIC')
    lat		= dataset.createDimension('lat'     , nlats           )
    lon		= dataset.createDimension('lon'     , nlons           )
    lead_nc	= dataset.createDimension('lead'    , nleads          )
    week_nc	= dataset.createDimension('week'    , None            )
    lag_nc	= dataset.createDimension('lag'     , nlags_ecmf      )
    ens_nc	= dataset.createDimension('ensemble', nmembers_ecmf   )

    # create variables
    latitudes	= dataset.createVariable('latitude' ,'f'       ,('lat',)     )
    longitudes	= dataset.createVariable('longitude','f'       ,('lon',)     )
    leads	= dataset.createVariable('lead'     ,np.float64,('lead',)    )
    weeks	= dataset.createVariable('week'     ,np.float64,('week',)    )
    lags	= dataset.createVariable('lag'      ,np.float64,('lag',)     )
    ensembles	= dataset.createVariable('ensemble' ,np.float64,('ensemble',))
    p_mean	= dataset.createVariable('week_precip','d',('week','lead','lag','ensemble','lat','lon'))

    # Global Attributes (will need modified accordingly)
    dataset.description	= 'ECMWF lagged ensembles weekly mean precipitation'
    dataset.history	= 'Created ' + tt.ctime(tt.time())
    dataset.source	= 'Subset by L. Guo'
    # Variable Attributes
    latitudes.units	= 'degrees_north'
    longitudes.units	= 'degrees_east'
    f_week		= starts_ukmo[0]
    weeks.units		= 'weeks since '+str(y)+'-'+month[m]+'-01'
    weeks.calendar	= 'none'
    p_mean.units	= 'kg m-2 / day'
    p_mean.long_name	= 'Weekly mean precipitation'

    # Fill variables with data
    latitudes[:]	= ecmf_lat
    longitudes[:]	= ecmf_lon
    p_mean[:]		= np.swapaxes(all_ecmf,0,1)
    leads[:]		= np.arange(1,nleads+1)
    weeks[:]		= np.arange(1,len(starts_ukmo)+1)
    lags[:]		= np.arange(0,nlags_ecmf)
    ensembles[:]	= np.arange(1,nmembers_ecmf+1)
    dataset.close()
    print 'saved ecmf file for '+str(y)+'-'+month[m]+'-'+starts_ukmo[s]

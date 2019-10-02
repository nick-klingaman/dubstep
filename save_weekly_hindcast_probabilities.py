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
dir_out	= dir_in+'hindcast_probabilities/'


# number of lagged ensembles
ncep_lags = 7
ecmf_lags = 3

# number of ensemble members for each forecast
ukmo_nmembers = 7
ncep_nmembers = 4
ecmf_nmembers = 11
bam_nmembers = 11

ncep_tmembers = ncep_nmembers*ncep_lags
ecmf_tmembers= ecmf_nmembers*ecmf_lags

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

percentiles = [100./3.,200./3.]
p_name = 'tercile' # give categories an appropriate name

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

'''
Compute forecast probabilities
'''
fname_p_chirps = dir_out+region+'_chirps_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
fname_p_ukmo = dir_out+region+'_UKMO_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
fname_p_ncep = dir_out+region+'_NCEP_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
fname_p_ecmf = dir_out+region+'_ECMWF_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'

fname_p_bam = dir_out+region+'_BAM_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'
fname_p_chirps_bam = dir_out+region+'_CHIRPS_for_BAM_p_'+p_name+'_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])+'.npy'


if os.path.exists(fname_p_ukmo) == True:

  p_chirps = np.zeros((len(percentiles)+1,nleads,nweeks,len(years),len(all_lat),len(all_lon)))
  p_ukmo = np.zeros((len(percentiles)+1,2,nleads,nweeks,len(years),len(all_lat),len(all_lon)))
  p_ncep = np.zeros((len(percentiles)+1,2,nleads,nweeks,len(years),len(all_lat),len(all_lon)))
  p_ecmf = np.zeros((len(percentiles)+1,2,nleads,nweeks,len(years),len(all_lat),len(all_lon)))

  if season == 'NDJFM':
    p_chirps_for_bam = np.zeros((len(percentiles)+1,nleads,5,2,len(years),len(all_lat),len(all_lon)))
    p_bam = np.zeros((len(percentiles)+1,2,nleads,5,2,len(years),len(all_lat),len(all_lon)))

  for l in np.arange(nleads):
    for y in np.arange(0,len(years)):
      if season == 'NDJFM':
        for wbam in np.arange(0,5):
          for n in [0,1]:
            if week_mid_dates_bam[y,l,wbam,n,1] in keep_mon:
              tmp_chirps = []
              tmp_chirps = np.delete(all_chirps_bam[:,l,wbam,n,:,:],y,axis=0)
              chirps_percentiles_bam = np.nanpercentile(tmp_chirps,percentiles,axis=0)
              curr_chirps_bam = np.copy(all_chirps_bam[y,l,wbam,n,:,:])
              nnan1 = np.sum(np.isnan(curr_chirps_bam))

              if (nnan1 < (len(all_lon)*len(all_lat))):
                # 1. Compute observed probability (will be 1 or 0)
                for nn in np.arange(0,len(percentiles)+1):
                  if nn == 0:
                    p_chirps_for_bam[nn,l,wbam,n,y,:,:] = mask_percentiles(curr_chirps_bam,chirps_percentiles_bam[nn,:,:],'below')
                  elif nn == len(percentiles):
                    p_chirps_for_bam[nn,l,wbam,n,y,:,:] = mask_percentiles(curr_chirps_bam,chirps_percentiles_bam[nn-1,:,:],'above')
                  else:
                    p_chirps_for_bam[nn,l,wbam,n,y,:,:] = mask_percentiles(curr_chirps_bam,chirps_percentiles_bam[nn-1:nn+1,:,:],'between')


              # compute the probability that the ensemble was
              # above/below each tercile
              for e in np.arange(0,bam_nmembers):
                #get percentiles for UKMO
                tmp_bam = []
                tmp_bam = np.delete(all_bam_ens[:,l,wbam,n,e,:,:],y,axis=0)
                bam_percentiles = np.nanpercentile(tmp_bam,percentiles,axis=0)# all_chirps in (11, 5, 20, 40, 47); chirps_percentiles in (2, 40, 47)
                curr_bam = np.copy(all_bam_ens[y,l,wbam,n,e,:,:].squeeze())
                nnan2 = np.sum(np.isnan(curr_bam))
                if (nnan2 < (len(all_lon)*len(all_lat))):
                  for nn in np.arange(0,len(percentiles)+1):
                    if nn == 0:
                      p_bam[nn,0,l,wbam,n,y,:,:] = p_bam[nn,0,l,wbam,n,y,:,:] + mask_percentiles(curr_bam,chirps_percentiles_bam[nn,:,:],'below')
                      p_bam[nn,1,l,wbam,n,y,:,:] = p_bam[nn,1,l,wbam,n,y,:,:] + mask_percentiles(curr_bam,bam_percentiles[nn,:,:],'below')
                    elif nn == len(percentiles):
                      p_bam[nn,0,l,wbam,n,y,:,:] = p_bam[nn,0,l,wbam,n,y,:,:] + mask_percentiles(curr_bam,chirps_percentiles_bam[nn-1,:,:],'above')
                      p_bam[nn,1,l,wbam,n,y,:,:] = p_bam[nn,1,l,wbam,n,y,:,:] + mask_percentiles(curr_bam,bam_percentiles[nn-1,:,:],'above')
                    else:
                      p_bam[nn,0,l,wbam,n,y,:,:] = p_bam[nn,0,l,wbam,n,y,:,:] + mask_percentiles(curr_bam,chirps_percentiles_bam[nn-1:nn+1,:,:],'between')
                      p_bam[nn,1,l,wbam,n,y,:,:] = p_bam[nn,1,l,wbam,n,y,:,:] + mask_percentiles(curr_bam,bam_percentiles[nn-1:nn+1,:,:],'between')


      nw = 0
      for w in np.arange(0,nweeks):
        print 'week '+str(w+1)+' of '+str(nweeks)
        if l == 0:
          nw = nw + 1
        if week_mid_dates[y,l,w,1] in keep_mon:
          # compute the qth percentile of the data along the specified axis, while ignoring nan values.
          # here, it means that find 33% and 66% of chirps weekly mean precipitation value at each grid cell.

          # for BSS analysis, we should be removing the current year being assessed
          # from the calculation of percentiles
          tmp_chirps = []
          tmp_chirps = np.delete(all_chirps[:,l,w,:,:],y,axis=0)
          chirps_percentiles = np.nanpercentile(tmp_chirps,percentiles,axis=0)
          curr_chirps = np.copy(all_chirps[y,l,w,:,:])

          # 1. Compute observed probability (will be 1 or 0)
          for nn in np.arange(0,len(percentiles)+1):
            if nn == 0:
              p_chirps[nn,l,w,y,:,:] = mask_percentiles(curr_chirps,chirps_percentiles[nn,:,:],'below')
            elif nn == len(percentiles):
              p_chirps[nn,l,w,y,:,:] = mask_percentiles(curr_chirps,chirps_percentiles[nn-1,:,:],'above')
            else:
              p_chirps[nn,l,w,y,:,:] = mask_percentiles(curr_chirps,chirps_percentiles[nn-1:nn+1,:,:],'between')

          # compute the probability that the ensemble was
          # above/below each tercile
          for e in np.arange(0,ecmf_nmembers):
            if e < 7: # 1) do UKMO
              #get percentiles for UKMO
              tmp_ukmo = []
              tmp_ukmo = np.delete(all_ukmo_ens[:,l,w,e,:,:],y,axis=0)
              ukmo_percentiles = np.nanpercentile(tmp_ukmo,percentiles,axis=0)# all_chirps in (11, 5, 20, 40, 47); chirps_percentiles in (2, 40, 47)
              curr_ukmo = np.copy(all_ukmo_ens[y,l,w,e,:,:].squeeze())

              for nn in np.arange(0,len(percentiles)+1):
                if nn == 0:
                  p_ukmo[nn,0,l,w,y,:,:] = p_ukmo[nn,0,l,w,y,:,:] + mask_percentiles(curr_ukmo,chirps_percentiles[nn,:,:],'below')
                  p_ukmo[nn,1,l,w,y,:,:] = p_ukmo[nn,1,l,w,y,:,:] + mask_percentiles(curr_ukmo,ukmo_percentiles[nn,:,:],'below')
                elif nn == len(percentiles):
                  p_ukmo[nn,0,l,w,y,:,:] = p_ukmo[nn,0,l,w,y,:,:] + mask_percentiles(curr_ukmo,chirps_percentiles[nn-1,:,:],'above')
                  p_ukmo[nn,1,l,w,y,:,:] = p_ukmo[nn,1,l,w,y,:,:] + mask_percentiles(curr_ukmo,ukmo_percentiles[nn-1,:,:],'above')
                else:
                  p_ukmo[nn,0,l,w,y,:,:] = p_ukmo[nn,0,l,w,y,:,:] + mask_percentiles(curr_ukmo,chirps_percentiles[nn-1:nn+1,:,:],'between')
                  p_ukmo[nn,1,l,w,y,:,:] = p_ukmo[nn,1,l,w,y,:,:] + mask_percentiles(curr_ukmo,ukmo_percentiles[nn-1:nn+1,:,:],'between')

            if e < 4:
              for nl in np.arange(0,ncep_lags):
                # loop through ncep lags
                #get percentiles for NCEP
                tmp_ncep = []
                tmp_ncep = np.delete(all_ncep_ens[:,l,w,nl,e,:,:],y,axis=0)
                ncep_percentiles = np.nanpercentile(tmp_ncep,percentiles,axis=0)
                curr_ncep = np.copy(all_ncep_ens[y,l,w,nl,e,:,:])
                for nn in np.arange(0,len(percentiles)+1):
                  if nn == 0:
                    p_ncep[nn,0,l,w,y,:,:] = p_ncep[nn,0,l,w,y,:,:] + mask_percentiles(curr_ncep,chirps_percentiles[nn,:,:],'below')
                    p_ncep[nn,1,l,w,y,:,:] = p_ncep[nn,1,l,w,y,:,:] + mask_percentiles(curr_ncep,ncep_percentiles[nn,:,:],'below')
                  elif nn == len(percentiles):
                    p_ncep[nn,0,l,w,y,:,:] = p_ncep[nn,0,l,w,y,:,:] + mask_percentiles(curr_ncep,chirps_percentiles[nn-1,:,:],'above')
                    p_ncep[nn,1,l,w,y,:,:] = p_ncep[nn,1,l,w,y,:,:] + mask_percentiles(curr_ncep,ncep_percentiles[nn-1,:,:],'above')
                  else:
                    p_ncep[nn,0,l,w,y,:,:] = p_ncep[nn,0,l,w,y,:,:] + mask_percentiles(curr_ncep,chirps_percentiles[nn-1:nn+1,:,:],'between')
                    p_ncep[nn,1,l,w,y,:,:] = p_ncep[nn,1,l,w,y,:,:] + mask_percentiles(curr_ncep,ncep_percentiles[nn-1:nn+1,:,:],'between')

            for el in np.arange(0,ecmf_lags):
              #get percentiles for ECMWF
              tmp_ecmf = []
              tmp_ecmf = np.delete(all_ecmf_ens[:,l,w,el,e,:,:],y,axis=0)
              ecmf_percentiles = np.nanpercentile(tmp_ecmf,percentiles,axis=0)
              curr_ecmf = np.copy(all_ecmf_ens[y,l,w,el,e,:,:])
              for nn in np.arange(0,len(percentiles)+1):
                if nn == 0:
                  p_ecmf[nn,0,l,w,y,:,:] = p_ecmf[nn,0,l,w,y,:,:] + mask_percentiles(curr_ecmf,chirps_percentiles[nn,:,:],'below')
                  p_ecmf[nn,1,l,w,y,:,:] = p_ecmf[nn,1,l,w,y,:,:] + mask_percentiles(curr_ecmf,ecmf_percentiles[nn,:,:],'below')
                elif nn == len(percentiles):
                  p_ecmf[nn,0,l,w,y,:,:] = p_ecmf[nn,0,l,w,y,:,:] + mask_percentiles(curr_ecmf,chirps_percentiles[nn-1,:,:],'above')
                  p_ecmf[nn,1,l,w,y,:,:] = p_ecmf[nn,1,l,w,y,:,:] + mask_percentiles(curr_ecmf,ecmf_percentiles[nn-1,:,:],'above')
                else:
                  p_ecmf[nn,0,l,w,y,:,:] = p_ecmf[nn,0,l,w,y,:,:] + mask_percentiles(curr_ecmf,chirps_percentiles[nn-1:nn+1,:,:],'between')
                  p_ecmf[nn,1,l,w,y,:,:] = p_ecmf[nn,1,l,w,y,:,:] + mask_percentiles(curr_ecmf,ecmf_percentiles[nn-1:nn+1,:,:],'between')


  # At end compute actual Probability
  pf_ukmo = p_ukmo/float(ukmo_nmembers) # probability
  pf_ncep = p_ncep/float(ncep_tmembers) # probability, normalise by N ensemble members
  pf_ecmf = p_ecmf/float(ecmf_tmembers) # probability

  if season == 'NDJFM':
    pf_bam = p_bam/float(bam_nmembers) # bam probability

  # save probabilities
  np.save(fname_p_chirps,p_chirps)
  np.save(fname_p_ukmo,pf_ukmo)
  np.save(fname_p_ncep,pf_ncep)
  np.save(fname_p_ecmf,pf_ecmf)

  if season == 'NDJFM':
    np.save(fname_p_chirps_bam,p_chirps_for_bam)
    np.save(fname_p_bam,pf_bam)

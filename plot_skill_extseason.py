'''
Makes figures illustrating skill of weekly
S2S precipitation forecasts from BAM, UKMO,
NCEP and ECMWF over a given region
evaluated against CHIRPS precipitation

The script
1) Reads in model skill metrics computed from 'analyse_skill_extseason.py'
2) Plots figures
  - regional average skill vs lead-time
  - spatial maps of skill scores

Skill scores plotted are:
-Bias
-Mean precip in obs and models
-Anomaly correlation coefficient
-Separate RMSE's of ensemble members and ensemble mean
-Brier Skill Score for upper tercile precipitation events

M. Young 03/03/2020
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


def matt_t_test_corr(corr,n,n_effective):
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
  # p = 1 - stats.t.cdf(t_stat,df=n_effective-2)
  # p_original = 1 - stats.t.cdf(t_stat_original,df=n-2)

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


'''
Forecast probabilities
'''
percentiles = [100./3.,200./3.]
p_name = 'tercile' # give categories an appropriate name
p_names = ['dry','normal','wet']
p_cat = 0.333 # probability of falling into tercile

'''
Load evaluation statistics
'''
# save mean
fname_mean_chp = dir_in+'CHIRPS_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
week_chirps_mean = np.load(fname_mean_chp+'.npy')
# save bias
fname_bias = dir_in+'Bias_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
bias = np.load(fname_bias+'.npy')
# save ACC
fname_acc = dir_in+'ACC_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
acc = np.load(fname_acc+'.npy')
# save significance
fname_t = dir_in+'ttest_significance_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
t_signif = np.load(fname_t+'.npy')
# save RMSE
fname_rmse = dir_in+'RMSE_anom_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
rmse_anom = np.load(fname_rmse+'.npy')
# save RMSE
fname_rmse1 = dir_in+'RMSE_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
rmse = np.load(fname_rmse1+'.npy')
# save RMSE
fname_rmse2 = dir_in+'RMSE_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
rmse_mean = np.load(fname_rmse2+'.npy')
# save BSS
fname_bss = dir_in+'BSS_tercile_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
model_bss = np.load(fname_bss+'.npy')


fname_ukmo = dir_in+'UKMO_week_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
week_ukmo_mean = np.load(fname_ukmo+'.npy')
fname_ncep = dir_in+'NCEP_week_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
week_ncep_mean = np.load(fname_ncep+'.npy')
fname_emcf = dir_in+'ECMWF_week_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
week_ecmf_mean=np.load(fname_ecmf+'.npy')
fname_bam = dir_in+'BAM_week_mean_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
week_bam_mean = np.load(fname_bam+'.npy')
week_mean = [week_bam_mean,week_ecmf_mean,week_ncep_mean,week_ukmo_mean]


land_sea = np.nanmean(week_chirps_mean,axis=(0,1))
land_sea[land_sea>=0] = 1
land_sea[land_sea<0] = np.nan
# Create a dry mask (chirps seasonal clim < 1 mmd)
# this is used to mask Brier Skill Score later
# as BSS results are dodgy over dry regions
dry_mask = np.nanmean(np.copy(week_chirps_mean),axis=(0,1))
dry_mask[dry_mask <  1] = np.nan
dry_mask[dry_mask >= 1] = 1


'''
Average ACC, RMSE and BIAS over  sub-regions at each lead time
'''
rg_names = ['NSA','AMZ','NDE','SESA','AND','PAT']
rg_lon_min = [-80,-67,-47,-60,-75,-75]
rg_lon_max = [-50,-47,-34,-48,-67,-60]
rg_lat_min = [0,-15,-15,-35,-40,-50]
rg_lat_max = [12,-5,-5,-22,-15,-40]
acc_region	= np.nan*np.zeros((len(rg_names),4,nleads))
sig_region	= np.nan*np.zeros((len(rg_names),4,4,nleads))
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
  sig_region[rg,:,:,:] = np.nanmean(t_signif[:,:,:,subregion_lat_id,:][:,:,:,:,subregion_lon_id],axis=(3,4))
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


fmt_ls  = ['-s','-P','-X','-o']
ms_ls   = [7,9,9,8]
mw_ls   = [1,1,1,1]
fname_plot = dir_out+region+'_S2S_weekly_ACC_sig_'+season+'_subregions_average_'+str(years[0])+'_'+str(years[len(years)-1])
size1 = [8,4]
fig = plt.figure(figsize=(size1[0],size1[1]))
for rg in np.arange(0,len(rg_names)):
  plt.subplot(nrow,ncol,rg+1)
  for m in np.arange(0,len(model_ls)):
    tmp_sig = []
    tmp_sig = np.copy(sig_region[rg,:,m,:])
    sig_msk = np.nan*np.zeros(len(tmp_sig[0,:]))
    msk = tmp_sig[2,:] < tmp_sig[3,:]
    sig_msk[msk==True] = 1

    plt.plot(leads,acc_region[rg,m,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.7,mew=mw_ls[m])
    plt.plot(leads,acc_region[rg,m,:]*sig_msk,fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m], markeredgecolor=col_ls[m], markerfacecolor='w',label='_nolegend_')
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


fmt_ls	= ['-s','-+','-x','-o']
ms_ls	= [6,10,10,8]
mw_ls	= [1,2,2,1]
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


fmt_ls	= ['-s','-+','-x','-o']
ms_ls	= [6,10,10,8]
mw_ls	= [1,2,2,1]
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

fmt_ls  = ['-s','-P','-X','-o']
ms_ls   = [7,9,9,8]
mw_ls   = [1,1,1,1]
long_lead = np.arange(0,7)
for rg in np.arange(0,len(rg_names)):
  fname_plot = dir_out2+region+'_S2S_weekly_CC_sig_'+season+'_region_average_'+rg_names[rg]
  fig = plt.figure(figsize=(2.5,2))
  # plt.plot(long_lead,0.2*np.ones(len(long_lead)),'--k',alpha=0.5)
  for m in np.arange(0,len(model_ls)):
    # ms - markersize
    # mew - markeredgewidth
    # alpha - 0.0 transparent; 1.0 opaque
    tmp_sig = []
    tmp_sig = np.copy(sig_region[rg,:,m,:])
    sig_msk = np.nan*np.zeros(len(tmp_sig[0,:]))
    msk = tmp_sig[2,:] < tmp_sig[3,:]
    sig_msk[msk==True] = 1
    plt.plot(leads,acc_region[rg,m,:],fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m],alpha=0.7,mew=mw_ls[m])
    plt.plot(leads,acc_region[rg,m,:]*sig_msk,fmt_ls[m],color=col_ls[m],linewidth=2,ms=ms_ls[m], markeredgecolor=col_ls[m], markerfacecolor='w',label='_nolegend_')

  plt.ylim([0,0.8])
  plt.xlim([0.8,5.2])
  plt.xticks(np.arange(1,nleads+1,1))
  #plt.text(3.5,0.7,label_ls[rg]+' '+rg_names[rg],fontsize=12)#fontweight='bold'
  plt.ylabel('CC')
  plt.xlabel('Lead time (weeks)')
  plt.savefig(fname_plot+'.ps',bbox_inches='tight',dpi=300)
  #plt.show()
  plt.close()



fmt_ls	= ['-s','-+','-x','-o']
ms_ls	= [6,10,10,8]
mw_ls	= [1,2,2,1]
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
  plt.ylim([1.5,8])
  plt.xlim([0.8,5.2])
  plt.xticks(np.arange(1,nleads+1,1))
  plt.yticks(np.arange(1,9,1))
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

'''
# plot Autocorrelation
cols = 'RdYlBu'
cmin = -0.8
cmax = 0.8
cspc = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'Autocorrelation'
norm = BoundaryNorm(boundaries=clevs, ncolors=256)
m_number=[1,6,11,16]
lw = 1
nrow = len(model_ls)
ncol = nleads
size1=[8,7.5]
gl = 20
auto_ls = [4,1,2,3]#[3,0,1,2]

fname_plot = dir_out+region+'_S2S_weekly_auto_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
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
    uncal = mymap.pcolormesh(x,y,auto_corr[auto_ls[m],n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
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

chp_ls = ['CHIRPS','CHIRPS(BAM)']
fname_plot = dir_out+region+'_S2S_weekly_auto_chirps_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
fig = plt.figure(figsize=(8,6))
for n in np.arange(0,nleads):
  plt.subplot(2,4,n+1)
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
  uncal = mymap.pcolormesh(x,y,auto_corr[0,n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title(chp_ls[0]+' Week '+str(n+1))

  plt.subplot(2,4,n+4)
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
  uncal = mymap.pcolormesh(x,y,auto_corr[5,n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
  plt.title(chp_ls[1]+' Week '+str(n+1))

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
fname_plot = dir_out+region+'_S2S_weekly_CC_sig_auto_'+season+'_'+str(years[0])+'_'+str(years[len(years)-1])
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
    # for i, (x, y) in enumerate([(x,y) for x in range(LON) for y in range(LAT)]):
    for i in np.arange(0,len(all_lat)):
      for j in np.arange(0,len(all_lon)):
        if t_signif[0,m,n,i,j] == 1:
          plt.annotate('x',(all_lon[j],all_lat[i]), xytext=(all_lon[j],all_lat[i]), fontsize=2, color='k')
    plt.title(model_ls[m]+' Week '+str(n+1))

#plt.suptitle(season,fontweight='bold',fontsize=14,horizontalalignment='left')
plt.tight_layout(pad=2,w_pad=0.02,h_pad=0.02)
fig.subplots_adjust(bottom=0.08)
cbar_pos = [0.25, 0.03, 0.5, 0.015] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,orientation='horizontal',extend='both')
plt.savefig(fname_plot+'.png',bbox_inches='tight')
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
    fname_plot = dir_out2+region+'_S2S_weekly_CC_sig_auto_'+season+'_'+model_ls[m]+'_Week'+str(n+1)
    fig = plt.figure(figsize=(2,2))
    mymap = Basemap(projection='cyl',resolution='l',llcrnrlat=latlim[0],urcrnrlat=latlim[1],llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1])
    mymap.drawcoastlines(linewidth=lw)
    #mymap.drawcountries(linewidth=lw)
    mymap.drawparallels(np.arange(-90,90,gl),labels=[1,0,0,0],labelstyle='+/-')
    mymap.drawmeridians(np.arange(0,360,gl),labels=[0,0,0,1],labelstyle='+/-')
    x, y = mymap(*np.meshgrid(all_lon,all_lat))
    uncal = mymap.pcolormesh(x,y,acc[m,n,:,:],vmin=cmin,vmax=cmax,cmap=cols,norm=norm)
    for i in np.arange(0,len(all_lat)):
      for j in np.arange(0,len(all_lon)):
        if t_signif[0,m,n,i,j] == 1:
          plt.annotate('x',(all_lon[j],all_lat[i]), xytext=(all_lon[j],all_lat[i]), fontsize=3, color='k')
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

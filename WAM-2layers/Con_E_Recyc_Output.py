# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:24:45 2016

@author: Ent00002
"""

# delayed runs, comment the 2 lines below when unused
# import time
# time.sleep(100)

#%% Import libraries

import numpy as np
import scipy.io as sio
import calendar
import datetime
from getconstants import getconstants
from timeit import default_timer as timer
import os

#%%BEGIN OF INPUT1 (FILL THIS IN)
years		= np.arange(1979,2016)	# Note that this can be forward as tracking was already finished, whereas the masterscript must be in backward order
yearpart	= np.arange(0,366)	# for a full (leap)year fill in np.arange(0,366)
daily		= 0			# 1 for writing out daily data, 0 for only monthly data
timetracking	= 1			# 0 for not tracking time and 1 for tracking time

# Manage the extent of your dataset (FILL THIS IN)
# Define the latitude and longitude cell numbers to consider and corresponding lakes that should be considered part of the land
latnrs = np.arange(7,114)
lonnrs = np.arange(0,240)

# the lake numbers below belong to the ERA-Interim data on 1.5 degree starting at Northern latitude 79.5 and longitude -180
lake_mask_1 = np.array([9,9,9,12,12,21,21,22,22,23,24,25,23,23,25,25,53,54,61,23,24,23,24,25,27,22,23,24,25,26,27,28,22,25,26,27,28,23,23,12,18])
lake_mask_2 = np.array([120+19,120+40,120+41,120+43,120+44,120+61,120+62,120+62,120+63,120+62,120+62,120+62,120+65,120+66,120+65,120+66,142-120,142-120,143-120,152-120,152-120,153-120,153-120,153-120,153-120,154-120,154-120,154-120,154-120,154-120,154-120,154-120,155-120,155-120,155-120,155-120,155-120,159-120,160-120,144-120,120+55])
lake_mask = np.transpose(np.vstack((lake_mask_1,lake_mask_2))) #recreate the arrays of the matlab model

# obtain the constants
invariant_data = '/home/users/lguo/CSSP/WAM2layersPython-master/download_scripts/invariants.nc'	#invariants
latitude,longitude,lsm,g,density_water,timestep,A_gridcell,L_N_gridcell,L_S_gridcell,L_EW_gridcell,gridcell = getconstants(latnrs,lonnrs,lake_mask,invariant_data)

interdata_folder = r'/home/users/lguo/CSSP/WAM2layersPython-master/interdata'
output_folder = r'/home/users/lguo/CSSP/WAM2layersPython-master/output'
sub_interdata_folder = os.path.join(interdata_folder, 'cn1_backward')

#END OF INPUT

#%% Datapaths (FILL THIS IN)


def data_path(y,a,years,timetracking):
    load_Sa_track = os.path.join(sub_interdata_folder, str(y) + '-' + str(a) + 'Sa_track.mat')
    load_Sa_time = os.path.join(sub_interdata_folder, str(y) + '-' + str(a) + 'Sa_time.mat')
    load_fluxes_and_storages = os.path.join(interdata_folder, str(y) + '-' + str(a) + 'fluxes_storages.mat')

    save_path = os.path.join(output_folder, 'E_track_cn1_full' + str(years[0]) + '-' + str(years[-1]) + '-timetracking' + str(timetracking) + '.mat')
    save_path_daily = os.path.join(output_folder, 'E_track_cn1_daily_full' + str(y) + '-timetracking' + str(timetracking) + '.mat')

    return load_Sa_track,load_Sa_time,load_fluxes_and_storages,save_path,save_path_daily


#%% Runtime & Results
start1 = timer()
startyear = years[0]

E_per_year_per_month		= np.zeros((len(years),12,len(latitude),len(longitude)))
E_track_per_year_per_month	= np.zeros((len(years),12,len(latitude),len(longitude)))
P_per_year_per_month		= np.zeros((len(years),12,len(latitude),len(longitude)))
Sa_track_down_per_year_per_month= np.zeros((len(years),12,len(latitude),len(longitude)))
Sa_track_top_per_year_per_month = np.zeros((len(years),12,len(latitude),len(longitude)))
W_down_per_year_per_month	= np.zeros((len(years),12,len(latitude),len(longitude)))
W_top_per_year_per_month	= np.zeros((len(years),12,len(latitude),len(longitude)))
north_loss_per_year_per_month	= np.zeros((len(years),12,1,len(longitude)))
south_loss_per_year_per_month	= np.zeros((len(years),12,1,len(longitude)))
down_to_top_per_year_per_month	= np.zeros((len(years),12,len(latitude),len(longitude)))
top_to_down_per_year_per_month	= np.zeros((len(years),12,len(latitude),len(longitude)))
water_lost_per_year_per_month	= np.zeros((len(years),12,len(latitude),len(longitude)))

if timetracking == 1:
    Sa_time_down_per_year_per_month	= np.zeros((len(years),12,len(latitude),len(longitude)))
    Sa_time_top_per_year_per_month	= np.zeros((len(years),12,len(latitude),len(longitude)))
    E_time_per_year_per_month		= np.zeros((len(years),12,len(latitude),len(longitude)))
    Sa_dist_down_per_year_per_month	= np.zeros((len(years),12,len(latitude),len(longitude)))
    Sa_dist_top_per_year_per_month	= np.zeros((len(years),12,len(latitude),len(longitude)))
    E_dist_per_year_per_month		= np.zeros((len(years),12,len(latitude),len(longitude)))

for i in range(len(years)):
    y = years[i]
    ly = int(calendar.isleap(y))
    final_time = 364+ly
    
    E_per_day			= np.zeros((365+ly,len(latitude),len(longitude)))
    E_track_per_day		= np.zeros((365+ly,len(latitude),len(longitude)))
    P_per_day			= np.zeros((365+ly,len(latitude),len(longitude)))
    Sa_track_down_per_day	= np.zeros((365+ly,len(latitude),len(longitude)))
    Sa_track_top_per_day	= np.zeros((365+ly,len(latitude),len(longitude)))
    W_down_per_day		= np.zeros((365+ly,len(latitude),len(longitude)))
    W_top_per_day		= np.zeros((365+ly,len(latitude),len(longitude)))
    north_loss_per_day		= np.zeros((365+ly,1,len(longitude)))
    south_loss_per_day		= np.zeros((365+ly,1,len(longitude)))
    down_to_top_per_day		= np.zeros((365+ly,len(latitude),len(longitude)))
    top_to_down_per_day		= np.zeros((365+ly,len(latitude),len(longitude)))
    water_lost_per_day		= np.zeros((365+ly,len(latitude),len(longitude)))
    if timetracking == 1:
        Sa_time_down_per_day	= np.zeros((365+ly,len(latitude),len(longitude)))
        Sa_time_top_per_day	= np.zeros((365+ly,len(latitude),len(longitude)))
        E_time_per_day		= np.zeros((365+ly,len(latitude),len(longitude)))
        Sa_dist_down_per_day	= np.zeros((365+ly,len(latitude),len(longitude)))
        Sa_dist_top_per_day	= np.zeros((365+ly,len(latitude),len(longitude)))
        E_dist_per_day		= np.zeros((365+ly,len(latitude),len(longitude)))
    
    for j in range(len(yearpart)):
        start = timer()
        a = yearpart[j]
        datapath = data_path(y,a,years,timetracking)
        if a > final_time: # a = 365 (366th index) and not a leapyear\
            pass
        else:
            # load tracked data [97,107,240]
            loading_ST = sio.loadmat(datapath[0],verify_compressed_data_integrity=False)
            Sa_track_top	= loading_ST['Sa_track_top']
            Sa_track_down	= loading_ST['Sa_track_down']
            north_loss		= loading_ST['north_loss']
            south_loss		= loading_ST['south_loss']
            down_to_top		= loading_ST['down_to_top']
            top_to_down		= loading_ST['top_to_down']
            water_lost		= loading_ST['water_lost']
            Sa_track		= Sa_track_top + Sa_track_down
            if timetracking == 1:
                loading_STT = sio.loadmat(datapath[1],verify_compressed_data_integrity=False)
                Sa_time_top	= loading_STT['Sa_time_top']
                Sa_time_down	= loading_STT['Sa_time_down']
                Sa_dist_top	= loading_STT['Sa_dist_top']
                Sa_dist_down	= loading_STT['Sa_dist_down']
            
            # load the total moisture data [96/97,107,240]
            loading_FS = sio.loadmat(datapath[2],verify_compressed_data_integrity=False)
            Fa_E_top	= loading_FS['Fa_E_top']
            Fa_N_top	= loading_FS['Fa_N_top']
            Fa_E_down	= loading_FS['Fa_E_down']
            Fa_N_down	= loading_FS['Fa_N_down']
            Fa_Vert	= loading_FS['Fa_Vert']
            E		= loading_FS['E']
            P		= loading_FS['P']
            W_top	= loading_FS['W_top']
            W_down	= loading_FS['W_down']
            W = W_top + W_down
            
            # compute tracked evaporation [96,107,240]
	    # LG: E_track is calculated via the bottom layer Sa.
            E_track = E[:,:,:] * (Sa_track_down[1:,:,:] / W_down[1:,:,:])
            
            # save per day
            E_per_day[a,:,:]			= np.sum(E, axis =0)
            E_track_per_day[a,:,:]		= np.sum(E_track, axis =0)
            P_per_day[a,:,:]			= np.sum(P, axis =0)
            Sa_track_down_per_day[a,:,:]	= np.mean(Sa_track_down[1:,:,:], axis =0)
            Sa_track_top_per_day[a,:,:]		= np.mean(Sa_track_top[1:,:,:], axis =0)
            W_down_per_day[a,:,:]		= np.mean(W_down[1:,:,:], axis =0)
            W_top_per_day[a,:,:]		= np.mean(W_top[1:,:,:], axis =0)
            
            north_loss_per_day[a,:,:]		= np.sum(north_loss, axis =0)
            south_loss_per_day[a,:,:]		= np.sum(south_loss, axis =0)
            down_to_top_per_day[a,:,:]		= np.sum(down_to_top, axis =0)
            top_to_down_per_day[a,:,:]		= np.sum(top_to_down, axis =0)
            water_lost_per_day[a,:,:]		= np.sum(water_lost, axis =0)
            
            if timetracking == 1:
                # compute tracked evaporation time [96,107,240]
		# LG: E_time is mean of two concecutive timesteps
                E_time = 0.5 * ( Sa_time_down[:-1,:,:] + Sa_time_down[1:,:,:] )					# seconds
                # save per day
                Sa_time_down_per_day[a,:,:]	= np.mean(Sa_time_down[:-1,:,:], axis=0)			# seconds
                Sa_time_top_per_day[a,:,:]	= np.mean(Sa_time_top[:-1,:,:], axis=0)				# seconds
                E_time_per_day[a,:,:]		= np.sum((E_time * E_track), axis = 0) / E_track_per_day[a,:,:]	# seconds
                # remove nans                
                where_are_NaNs = np.isnan(E_time_per_day)
                E_time_per_day[where_are_NaNs] = 0

                # compute tracked evaporation distance [96,107,240]
                E_dist = 0.5 * ( Sa_dist_down[:-1,:,:] + Sa_dist_down[1:,:,:] )					# m
                # save per day
                Sa_dist_down_per_day[a,:,:]	= np.mean(Sa_dist_down[:-1,:,:], axis=0)			# m
                Sa_dist_top_per_day[a,:,:]	= np.mean(Sa_dist_top[:-1,:,:], axis=0)				# m
                E_dist_per_day[a,:,:]		= np.sum((E_dist * E_track), axis = 0) / E_track_per_day[a,:,:]	# m
                # remove nans                
                where_are_NaNs = np.isnan(E_dist_per_day)
                E_dist_per_day[where_are_NaNs] = 0
        
        end = timer()
        print 'Runtime output for day ' + str(a+1) + ' in year ' + str(y) + ' is',(end - start),' seconds.'
    
    if daily == 1:
        if timetracking == 0: # create dummy values
            Sa_time_down_per_day	= 0
            Sa_time_top_per_day		= 0
            E_time_per_day		= 0
    
        sio.savemat(datapath[4],
                    {'E_per_day':E_per_day,'E_track_per_day':E_track_per_day,'P_per_day':P_per_day,
                     'Sa_track_down_per_day':Sa_track_down_per_day,'Sa_track_top_per_day':Sa_track_top_per_day, 
                     'Sa_time_down_per_day':Sa_time_down_per_day,'Sa_time_top_per_day':Sa_time_top_per_day, 
                     'Sa_dist_down_per_day':Sa_dist_down_per_day,'Sa_dist_top_per_day':Sa_dist_top_per_day, 
                     'W_down_per_day':W_down_per_day,'W_top_per_day':W_top_per_day,
                     'E_time_per_day':E_time_per_day,'E_dist_per_day':E_dist_per_day},do_compression=True)

    # values per month        
    for m in range(12):
        first_day = int(datetime.date(y,m+1,1).strftime("%j"))
        last_day = int(datetime.date(y,m+1,calendar.monthrange(y,m+1)[1]).strftime("%j"))
        days = np.arange(first_day,last_day+1)-1						# -1 because Python is zero-based
        
        E_per_year_per_month[y-startyear,m,:,:]			= (np.squeeze(np.sum(E_per_day[days,:,:], axis = 0)))
        E_track_per_year_per_month[y-startyear,m,:,:]		= (np.squeeze(np.sum(E_track_per_day[days,:,:], axis = 0)))
        P_per_year_per_month[y-startyear,m,:,:]			= (np.squeeze(np.sum(P_per_day[days,:,:], axis = 0)))
        Sa_track_down_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.mean(Sa_track_down_per_day[days,:,:], axis = 0)))
        Sa_track_top_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.mean(Sa_track_top_per_day[days,:,:], axis = 0)))
        W_down_per_year_per_month[y-startyear,m,:,:]		= (np.squeeze(np.mean(W_down_per_day[days,:,:], axis = 0)))
        W_top_per_year_per_month[y-startyear,m,:,:]		= (np.squeeze(np.mean(W_top_per_day[days,:,:], axis = 0)))
        north_loss_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.sum(north_loss_per_day[days,:,:], axis = 0)))
        south_loss_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.sum(south_loss_per_day[days,:,:], axis = 0)))
        down_to_top_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.sum(down_to_top_per_day[days,:,:], axis = 0)))
        top_to_down_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.sum(top_to_down_per_day[days,:,:], axis = 0)))
        water_lost_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.sum(water_lost_per_day[days,:,:], axis = 0)))
        
        if timetracking == 1:
            # tracked time
            Sa_time_down_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.mean(Sa_time_down_per_day[days,:,:], axis = 0)))
            Sa_time_top_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.mean(Sa_time_top_per_day[days,:,:], axis =0)))
            E_time_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.sum( E_time_per_day[days,:,:] 
                * E_track_per_day[days,:,:], axis=0)) / np.squeeze(E_track_per_year_per_month[y-startyear,m,:,:]))
                
            # remove nans                
            where_are_NaNs = np.isnan(E_time_per_year_per_month)
            E_time_per_year_per_month[where_are_NaNs] = 0

            # tracked distance
            Sa_dist_down_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.mean(Sa_dist_down_per_day[days,:,:], axis = 0)))
            Sa_dist_top_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.mean(Sa_dist_top_per_day[days,:,:], axis =0)))
            E_dist_per_year_per_month[y-startyear,m,:,:]	= (np.squeeze(np.sum( E_dist_per_day[days,:,:] 
                * E_track_per_day[days,:,:], axis=0)) / np.squeeze(E_track_per_year_per_month[y-startyear,m,:,:]))
                
            # remove nans                
            where_are_NaNs = np.isnan(E_dist_per_year_per_month)
            E_dist_per_year_per_month[where_are_NaNs] = 0
                
        elif timetracking == 0:
            Sa_time_down_per_year_per_month	= 0
            Sa_time_top_per_year_per_month	= 0
            E_time_per_year_per_month		= 0
            Sa_dist_down_per_year_per_month	= 0
            Sa_dist_top_per_year_per_month	= 0
            E_dist_per_year_per_month		= 0

# save monthly data
sio.savemat(datapath[3],
           {'E_per_year_per_month':E_per_year_per_month,'E_track_per_year_per_month':E_track_per_year_per_month,'P_per_year_per_month':P_per_year_per_month,
            'Sa_track_down_per_year_per_month':Sa_track_down_per_year_per_month,'Sa_track_top_per_year_per_month':Sa_track_top_per_year_per_month, 
            'Sa_time_down_per_year_per_month':Sa_time_down_per_year_per_month,'Sa_time_top_per_year_per_month':Sa_time_top_per_year_per_month, 
            'E_time_per_year_per_month':E_time_per_year_per_month, 'Sa_dist_down_per_year_per_month':Sa_dist_down_per_year_per_month,'Sa_dist_top_per_year_per_month':Sa_dist_top_per_year_per_month, 'E_dist_per_year_per_month':E_dist_per_year_per_month, 'W_down_per_year_per_month':W_down_per_year_per_month,'W_top_per_year_per_month':W_top_per_year_per_month,
            'north_loss_per_year_per_month':north_loss_per_year_per_month,'south_loss_per_year_per_month':south_loss_per_year_per_month,
            'down_to_top_per_year_per_month':down_to_top_per_year_per_month,'top_to_down_per_year_per_month':top_to_down_per_year_per_month,
            'water_lost_per_year_per_month':water_lost_per_year_per_month})

end1 = timer()
print 'The total runtime of Con_E_Recyc_Output is',(end1-start1),' seconds.'

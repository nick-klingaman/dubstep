"""

Analysing Scales of Precipitation diagnostics
Space-time Coherence package

These functions allow the user to compute and plot diagnostics of the
spatial and temporal coherence of precipitation in observations and/or
models.  These functions should work for any resolution or timescale.

These functions assume that the input precipitation dataset is in netCDF.

Written by Nicholas Klingaman
nicholas.klingaman@ncas.ac.uk

(C) The author 2017

"""

import numpy as np
import iris
import cf,cfplot as cfp
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

def default_options():
    
    """
    Dummy function to return initial values of key parameters that must be defined.
    """
    
    region_size=0
    input_box_size=0
    lag_length=0
    autocorr_length=0
    time_type='None'
    grid_type='None'
    time_desc='None'
    grid_desc='None'
    model_desc='None'
    return(region_size,input_box_size,lag_length,autocorr_length,time_type,grid_type,time_desc,grid_desc,model_desc)

def parameters():
    
    """
    Parameters that control the size of arrays used in other functions. 
    """

    max_box_distance=100
    max_timesteps=100
    max_boxes=100
    return(max_box_distance,max_boxes,max_timesteps)

def read_precip(model_dict):

    """
    Use iris to read precipitation data into a Cube.
    The area to read is controlled by the "region" element of the dataset dictionary (model_dict)
    "region" should be defined as min_lat, max_lat, min_lon, max_lon.

    Arguments:
    * model_dict: The dictionary containing details of the dataset to read.
      (See the create_dict function).

    Returns:
    * precip: An iris cube of precipitation data.
    """

    constraint = iris.Constraint(model_dict['constraint'],
                                     latitude=lambda cell: model_dict['region'][0] <= cell <= model_dict['region'][1],
                                     longitude=lambda cell: model_dict['region'][2] <= cell <= model_dict['region'][3])
    precip = iris.load_cube(model_dict['infile'],constraint)*model_dict['scale_factor']
    return(precip)

def compute_histogram(precip,bins):

    """
    Computes 1D and 2D histograms of precipitation from an input iris cube of precipitation.  
    The 2D histogram is the histogram of precipitation on consecutive timesteps. 
    
    Arguments:
    * precip : 
       An iris cube of precipitation
    * bins:
       A numpy array of the edges of the bins for which to compute the histogram,
       in the same units as "precip" above.

    Returns:
    * oned_hist:
       The one-dimensional histogram, as a numpy array, normalised so that the sum of the values is one.
    * twod_hist:
       The two-dimensional histogram, as a numpy array, normalised so that the sum of the values is one.
    """

    oned_hist, bin_edges = np.histogram(precip.data,bins)
    nbins=len(oned_hist)
    twod_hist=np.zeros((nbins,nbins))
    print len(oned_hist)
    print '---> Computing 2D histogram'
    for t, t_slice in enumerate(precip.slices(['time'])):
        next_slice = t_slice.copy()
        next_slice.data = np.roll(t_slice.data,1,0)
        twod_hist_temp,xedges,yedges = np.histogram2d(t_slice.data,next_slice.data,bins)
        twod_hist=twod_hist+twod_hist_temp
    twod_hist=twod_hist/np.sum(twod_hist)        
    oned_hist=oned_hist/np.float(np.sum(oned_hist))
    return(oned_hist,twod_hist)

def plot_histogram(oned_hist,twod_hist,model_dict,bins,title=True,colorbar=True):

    """
    Creates a PostScript plot of the 1D and 2D histograms calculated in compute_histogram.

    Arguments:
    * oned_hist: 
       The 1D histogram computed in compute_histogram.
    * twod_hist:
       The 2D histogram computed in compute_histogram.
    * model_dict:
       The dictioary containing details of this model.
    * bins:
       The edges of the histogram bins, as a numpy array.
    
    Optional arguments:
    * title:
       A logical to control whether the title is printed at the top of the plot (useful for creating multi-panel plots).
    * colorbar:
       A logical to control whether the colorbar is printed at the bottom of the plot (useful for creating multi-panel plots).
    """

    print '---> Plotting 2D histogram'
    nbins = np.size(oned_hist)
    hist_con_levs=[1e-5,2e-5,4e-5,7e-5,1e-4,2e-4,4e-4,7e-4,1e-3,2e-3,4e-3,7e-3,1e-2,2e-2,4e-2,7e-2,1e-1]
    fig = plt.figure(figsize=[9,8])
    ax = fig.add_subplot(111)
    cmap = plt.cm.get_cmap("viridis_r")
    norm = BoundaryNorm(hist_con_levs,ncolors=cmap.N,clip=True)
    contour = ax.pcolormesh(np.arange(nbins+1),np.arange(nbins+1),twod_hist,cmap=cmap,norm=norm)
    if colorbar == True:
        cbar = fig.colorbar(contour,orientation='horizontal',ticks=hist_con_levs)
        cbar.ax.set_xlabel('Probability',fontsize=18)
        cbar.ax.set_xticklabels(['1e-5','2e-5','3e-5','4e-5','7e-5','1e-4','2e-4','4e-4','7e-4','1e-3','2e-3','4e-3','7e-3','1e-2','7e-2','1e-1'])
    ax.set_xlabel('Precipitation at time t (mm day$^{-1}$)',fontsize=18)
    ax.set_ylabel('Precipitation at time t+1 (mm day$^{-1}$)',fontsize=18)
    ticklabels=['< '+str(bins[1])]
    for bin in xrange(1,nbins):
        ticklabels.append(str(bins[bin]))
    ticklabels.append(' > '+str(bins[nbins-1]))
    ax.set_xticks(np.arange(nbins+1))
    ax.set_xticklabels(ticklabels)
    ax.set_yticks(np.arange(nbins+1))
    ax.set_yticklabels(ticklabels)
    if title == True:
        title_string = '2D histogram for '+model_dict['legend_name']
        if 'time_desc' in model_dict:
            title_string = title_string+' '+model_dict['time_desc']
        if 'grid_desc' in model_dict:
            title_string = title_string+' '+model_dict['grid_desc']
        title_string = title_string + ' data'
        ax.set_title(title_string)
    ax.axis([0,nbins,0,nbins])
    ax.set_xlim(xmin=0,xmax=nbins)
    ax2 = ax.twinx()
    ax2.plot(np.arange(nbins)+0.5,oned_hist,'k--',marker='o',markersize=8)
    ax2.set_yscale('log',nonposy='clip')
    ax2.set_ylim(ymin=0.0009,ymax=1.0)
    ax2.set_ylabel('Probability of precipitation in bin',fontsize=18)
    ax2.set_yticks([1e-3,1.4e-3,2e-3,3e-3,4.5e-3,7e-3,1e-2,1.4e-2,2e-2,3e-2,4.5e-2,7e-2,1e-1,1.4e-1,2e-1,3e-1,4.5e-1,7e-1,1])
    ax2.set_yticklabels(['1.0e-3','1.4e-3','2.0e-3','3.0e-3','4.5e-3','7.0e-3','1.0e-2','1.4e-2','2.0e-2','3.0e-2','4.5e-2',
                         '7.0e-2','1.0e-1','1.4e-1','2.0e-1','3.0e-1','4.5e-1','7.0e-1','1.0e0'])
    ax2.set_xlim(xmin=0,xmax=nbins)
    plot_name='asop_coherence.'+model_dict['name']
    if 'grid_type' in model_dict:
        plot_name=plot_name+'_'+model_dict['grid_type']
    if 'time_type' in model_dict:
        plot_name=plot_name+'_'+model_dict['time_type']
    plot_name=plot_name+'_precip_twodpdf.ps'
    plt.savefig(plot_name,bbox_inches='tight')

def compute_equalgrid_corr(precip,model_dict):
    
    """
    Compute correlations in space and time, using the native spatial and temporal 
    resolutions of the input data.

    Method:
    The input spatial domain is broken down into non-overlapping regions of N x N gridpoints, where
    N is controlled by the "region_size" input parameter.

    Correlations in both space and time are computed with respect to lag-0 at the central 
    point in each region.  

    Composite (average) correlations are computed over all regions in the domain.

    To express the average correlations as a function of distance from the central point, bins of
    distance from the central point are created, with width delta_x (the x spacing of the input data, 
    taken from the dataset dictionary), starting from 0.5*delta_x.  Correlations are averaged within each bin,
    at each lag.
     
    Arguments:
     * precip:
        An iris cube of precipitation data
     * model_dict:
        The dictionary containing information about this dataset.

    Returns:
     * corr_map (lag_length,region_size,region_size):
        A 'map' of the composite correlations over all regions in the domain, at each lag.
     * lag_vs_distance (lag_length,region_size):
        The composite correlation over all regions in the domain, expressed a function of time (lag)
        and distance from the central point, in bins of delta_x (x spacing of the input data) starting
        from 0.5*delta_x.
     * autocorr (lag_length):
        The composite auto-correlation of precipitation at the central point, averaged over all regions
        in the domain.
     * npts (region_size):
        The number of gridpoints in each distance bin of the lag_vs_distance array.  Used by the 
        corresponding plotting function to determine whether there are any points in this bin.
    """
    
    lag_length = model_dict['lag_length']
    region_size = model_dict['region_size']
    print '---> Computing correlations for '+str(region_size)+'x'+str(region_size)+' sub-regions'
    nlon=len(precip.coord('longitude').points)
    nlat=len(precip.coord('latitude').points)
    print '----> Info: Size of domain in native gridpoints: '+str(nlon)+' longitude x '+str(nlat)+' latitude.'
    nregions=0
    npts=np.zeros(region_size,dtype=np.int32)
    distance_x=np.zeros(region_size)
    distance_y=np.zeros(region_size)
    corr_map=np.zeros((lag_length,region_size,region_size))
    lag_vs_distance=np.zeros((lag_length,region_size))
    autocorr=np.zeros((lag_length))
    for region_xstart in xrange(0,nlon-region_size+1,region_size):
        region_xcentre = region_xstart+region_size//2
        for region_ystart in xrange(0,nlat-region_size+1,region_size):
            region_ycentre = region_ystart+region_size//2
            central_precip = precip[:,region_ycentre,region_xcentre]
            for region_x in xrange(region_size):
                distance_x[region_x]=np.abs(region_xstart+region_x-region_xcentre)
            for region_y in xrange(region_size):
                distance_y[region_y]=np.abs(region_ystart+region_y-region_ycentre)*model_dict['dy']/float(model_dict['dx'])
            for region_x in xrange(region_size):
                for region_y in xrange(region_size):
                    distance=np.int(round(np.sqrt(distance_x[region_x]*distance_x[region_x]+distance_y[region_y]*distance_y[region_y])))-1
                    remote_precip=precip[:,region_y+region_ystart,region_x+region_xstart]
                    corr = np.corrcoef([central_precip.data,remote_precip.data])[1,0]
                    if corr == np.nan:
                        print 'NaN detected',region_xcentre,region_ycentre,region_x,region_y,central_precip.data,remote_precip.data
                    corr_map[0,region_y,region_x]=corr+corr_map[0,region_y,region_x]
                    if (region_x + region_xstart == region_xcentre) and (region_y + region_ystart == region_ycentre):
                        autocorr[0]=autocorr[0]+1
                        for lag in xrange(1,lag_length):
                            autocorr[lag]=np.corrcoef(central_precip.data,np.roll(central_precip.data,lag,0))[1,0]+autocorr[lag]
                    else:
                        lag_vs_distance[0,distance]=lag_vs_distance[0,distance]+corr
                        for lag in xrange(1,lag_length):
                            corr = np.corrcoef([central_precip.data,np.roll(remote_precip.data,lag,0)])[1,0]
                            if corr == np.nan:
                                print 'NaN detected',region_xcentre,region_ycentre,region_x,region_y,central_precip.data,remote_precip.data
                            corr_map[lag,region_y,region_x] = corr+corr_map[lag,region_y,region_x]
                            lag_vs_distance[lag,distance] = corr+lag_vs_distance[lag,distance]
                        npts[distance]=npts[distance]+1                                                
            nregions = nregions+1
    corr_map = corr_map/nregions
    print '----> Info: There are '+str(nregions)+' '+str(region_size)+'x'+str(region_size)+' sub-regions in your input data.'
    for lag in xrange(lag_length):
        for dist in xrange(region_size):
            if npts[dist] > 0:
                lag_vs_distance[lag,dist]=lag_vs_distance[lag,dist]/npts[dist]
            # If there are no gridpoints in range, set correlation to a missing value
            if npts[dist] == 0:
                lag_vs_distance[lag,dist]=-999
        autocorr[lag]=autocorr[lag]/(nregions)
    return(corr_map,lag_vs_distance,autocorr,npts)

def plot_equalgrid_corr(corr_map,lag_vs_distance,autocorr,npts,model_dict,title=True,colorbar=True):

    """
    Plots correlations as functions of space and time, which were first computed
    using compute_equalgrid_corr.
    
    Two types of plots are created:

        1. For each lag (from 0 to lag_length defined in model_dict), a 2D map
        of the composite correlations against the central point (at lag 0) for all
        points in the region (of length region_size).
        
        2. A single lag vs. distance plot showing the composite correlations against
        the central point (at lag 0), averaged over all points in each region in each
        distance bin (in steps of dx starting at 0.5dx), as well as the auto-correlation
        at the central point.
        
    See Fig. 2 in Klingaman et al. (2017, GMD, doi:10.5194/gmd-10-57-2017) for examples
    of these diagrams.

    Arguments:
      * corr_map (lag_length,region_size,region_size): 
         Composite maps of correlations at each lag, returned from compute_equalgrid_corr
      * lag_vs_distance (lag_length,region_size):
         Composite correlations over all regions in the domain, expressed as a function of
         time (lag) and distance from the central point, returned from compute_equalgrid_corr
      * autocorr (lag_length):
         The composite auto-correlation of precipitation at the central point, averaged over
         all regions in the domain.
      * npts (region_size):
         The number of gridpoints in each distance bin of the lag_vs_distance array.  Used to
         to determine whether there are any points in each distance bin.
         
    Optional arguments:
      * title:
         Include a title on the plot
      * colorbar:
         Include a colorbar on the plot
    
    Returns:
      None
    """

    region_size = model_dict['region_size']
    lag_length = model_dict['lag_length']
    
    print '---> Plotting correlation maps for '+str(region_size)+'x'+str(region_size)+' sub-regions'
    corr_con_levs=[0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95]
    
    # Plot correlation maps at each lag
    for lag in xrange(lag_length):
        plot_name='asop_coherence.'+model_dict['name']
        if 'grid_type' in model_dict:
            plot_name=plot_name+'_'+model_dict['grid_type']
        if 'time_type' in model_dict:
            plot_name=plot_name+'_'+model_dict['time_type']
        plot_name=plot_name+'_precip_'+str(region_size)+'x'+str(region_size)+'maps_lag'+str(lag)+'.ps'
        cfp.setvars(file=plot_name,text_fontsize=18,axis_label_fontsize=18)
        cfp.gopen(figsize=[8,8])
        cfp.gset(xmin=0,xmax=region_size,ymin=0,ymax=region_size)
        cfp.levs(manual=np.array(corr_con_levs))
        cfp.cscale(cmap='parula',reverse=1,ncols=len(corr_con_levs)+1,white=0)
        cfp.axes(xticks=np.arange(region_size-1)+0.5,yticks=np.arange(region_size)+0.5,
                 xticklabels=np.arange(region_size)-region_size//2,yticklabels=np.arange(region_size)-region_size//2,
                 xlabel='$\Delta$x (approximately '+str(model_dict['dx'])+' km at the equator)',
                 ylabel='$\Delta$y (approximately '+str(model_dict['dy'])+' km)')

        if title == True:
            title_string = 'Correlation map for '+model_dict['legend_name']
            if 'time_desc' in model_dict:
                title_string = title_string + ' - ' +model_dict['time_desc']
            if 'grid_desc' in model_dict:
                title_string = title_string + ' - ' +model_dict['grid_desc']
            title_string = title_string + ' - Lag '+str(lag)
        else:
            title_string = ''
        if colorbar == True:
            cfp.con(f=corr_map[lag,:,:],x=np.arange(region_size)+0.5,y=np.arange(region_size)+0.5,blockfill=1,
                    lines=False,line_labels=False,ptype=0,colorbar=1,colorbar_title='Correlation with (0,0) at lag 0, mean over all sub-regions',
                    title=title_string)
        else:
            cfp.con(f=corr_map[lag,:,:],x=np.arange(region_size)+0.5,y=np.arange(region_size)+0.5,blockfill=1,
                    lines=False,line_labels=False,ptype=0,colorbar=0,title=title_string)
        for region_x in xrange(region_size):
            for region_y in xrange(region_size):
                if corr_map[lag,region_y,region_x] > 0.5:
                    cfp.plotvars.plot.text(region_x+0.5,region_y+0.5,str(corr_map[lag,region_y,region_x])[0:4],
                                           horizontalalignment='center',color='white',fontsize=20)
                elif corr_map[lag,region_y,region_x] < 0.0:
                    cfp.plotvars.plot.text(region_x+0.5,region_y+0.5,str(corr_map[lag,region_y,region_x])[0:5],
                                           horizontalalignment='center',color='black',fontsize=20)
                else:
                    cfp.plotvars.plot.text(region_x+0.5,region_y+0.5,str(corr_map[lag,region_y,region_x])[0:4],
                                           horizontalalignment='center',color='black',fontsize=20)
        cfp.gclose()

    # Plot correlation vs. distance diagram
    print '---> Plotting lag vs. distance diagram'
    plot_name='asop_coherence.'+model_dict['name']
    if 'grid_type' in model_dict:
        plot_name=plot_name+'_'+model_dict['grid_type']
    if 'time_type' in model_dict:
        plot_name=plot_name+'_'+model_dict['time_type']
    plot_name=plot_name+'_precip_'+str(region_size)+'x'+str(region_size)+'_lag'+str(lag_length)+'.ps'
    cfp.setvars(file=plot_name,text_fontsize=20,axis_label_fontsize=18)
    cfp.gopen(figsize=[9,8])
    ticklabels=['Centre','0.5']
    max_dist=0
    for dist in xrange(2,region_size):
        if npts[dist] > 0 :
            ticklabels.append(str(dist-0.5))
            max_dist=dist
    ticklabels.append(str(max_dist+0.5))
    ticklabels.append(str(max_dist+1.5))
    lag_vs_distance=np.insert(lag_vs_distance,0,autocorr,1)
    cfp.gset(xmin=0,xmax=max_dist+1,ymin=0,ymax=lag_length)
    cfp.levs(manual=np.array(corr_con_levs))
    cfp.cscale(cmap='parula',reverse=1,ncols=len(corr_con_levs)+1,white=0)
    xtickvals=np.arange(max_dist+1)+2.0
    xtickvals=np.insert(xtickvals,0,[0.5,1.0])
    cfp.axes(xticks=xtickvals,yticks=np.arange(lag_length)+0.5,xticklabels=ticklabels,yticklabels=np.arange(lag_length),
             xlabel='$\Delta$x bins ($\Delta$x approximately '+str(model_dict['dx'])+' km at the equator)',ylabel='Lag')

    if title == True:
        title_string = 'Correlation map for '+model_dict['legend_name']
        if 'time_desc' in model_dict:
            title_string = title_string + ' - ' +model_dict['time_desc']
        if 'grid_desc' in model_dict:
            title_string = title_string + ' - ' +model_dict['grid_desc']
    else:
        title_string = ''
    if colorbar == True:
        cfp.con(f=lag_vs_distance[:,0:max_dist+2],x=np.arange(max_dist+2)+0.5,y=np.arange(lag_length)+0.5,blockfill=1,
                lines=False,line_labels=False,ptype=0,colorbar_title='Correlation with centre at lag=0, mean over all sub-regions',
                title=title_string)

    else:
        cfp.con(f=lag_vs_distance[:,0:max_dist+2],x=np.arange(max_dist+2)+0.5,y=np.arange(lag_length)+0.5,blockfill=1,
                lines=False,line_labels=False,ptype=0,colorbar=0,title=title_string)

    for dist in xrange(max_dist+2):
        for lag in xrange(lag_length):
            if lag_vs_distance[lag,dist] == -999:
                print '-999'
            #                cfp.plotvars.plot.text(dist+0.5,lag+0.5,'XXX',horizontalalignment='center',color='black',fontsize=20,verticalalignment='center')
            elif lag_vs_distance[lag,dist] > 0.5:
                cfp.plotvars.plot.text(dist+0.5,lag+0.5,str(lag_vs_distance[lag,dist])[0:4],
                                       horizontalalignment='center',color='white',fontsize=20)
            elif lag_vs_distance[lag,dist] < 0.0:
                cfp.plotvars.plot.text(dist+0.5,lag+0.5,str(lag_vs_distance[lag,dist])[0:5],
                                       horizontalalignment='center',color='black',fontsize=20)
            else:
                cfp.plotvars.plot.text(dist+0.5,lag+0.5,str(lag_vs_distance[lag,dist])[0:4],
                                       horizontalalignment='center',color='black',fontsize=20)
    cfp.gclose()

def compute_equalarea_corr(precip,model_dict):

    """
    Computes spatial correlations as a function of physical distance (in km).  Note that unlike
    compute_equalgrid_corr, this routine does *not* compute lagged correlations; it computes only
    instantaneous correlations.
    
    Method:
    As in compute_equalgrid_corr, the analysis domain is divided into square sub-regions that are
    box_size in length.  Precipitation at each point in the sub-region is correlated against the
    central point in the domain.  These correlations are binned by the physical distance from the
    central point, using bins delta_x wide starting from 0.5*delta_x.  Correlations are averaged
    within the bin and across all sub-regions.
    
    Limitations:

    * The physical distance in the x direction is computed using the latitude of the point at the centre
    of the box.  This is not entirely accurate, but likely accurate enough for these purposes
    (comparing correlations over 10s or 100s of km).

    * The number of distance bins is limited to max_box_distance, defined in the "parameters" function above.

    Arguments:
        * precip
            An iris cube of precipitation to analyse
        * model_dict
            The dictionary for this dataset
            
    Returns:
        * distance_correlations (max_box_distance):
            Correlations as a function of physical distance, binned by physical distance from the
            central point, using bins delta_x wide.
        * distance_ranges (3, max_box_distance):
            The minimum, median and maximum distance away from the central points, considering
            all points in that distance bin
        * distance_max (max_box_distance):
            The number of bins that contain valid points.  Can be used to subscript 
            "distance_correlations" and "distance_ranges" (i.e., distance_correlations[0:distance_max]).
    """
    
    
    print '---> Computing correlations for '+str(model_dict['box_size'])+'x'+str(model_dict['box_size'])+' km sub-boxes.'
    nlon=len(precip.coord('longitude').points)
    latitude = precip.coord('latitude').points
    nlat=len(latitude)
    max_box_distance,max_boxes,max_timesteps = parameters()
    box_length_x = np.int(model_dict['box_size']//model_dict['dx'])
    box_length_y = np.int(model_dict['box_size']//model_dict['dy'])
    distance_x = np.zeros(box_length_x)
    distance_y = np.zeros(box_length_y)
    nboxes=0
    autocorr=0
    npts=np.zeros(max_box_distance,dtype=np.int32)
    distance_lists = np.zeros((max_box_distance,max_box_distance*max_boxes))
    distance_ranges = np.zeros((3,max_box_distance))
    distance_correlations = np.zeros((max_box_distance))
    
    print '----> Info: Sub-boxes are '+str(box_length_x)+'x'+str(box_length_y)+' gridboxes in this model (nx x ny).'
    for box_xstart in xrange(0,nlon-box_length_x,box_length_x):
        box_xcentre = box_xstart+box_length_x//2
        for box_ystart in xrange(0,nlat-box_length_y,box_length_y):
            box_ycentre = box_ystart+box_length_y//2
            central_precip = precip[:,box_ycentre,box_xcentre]
            for box_x in xrange(box_length_x):
                distance_x[box_x]=np.abs(box_xstart+box_x-box_xcentre)*np.cos(latitude[box_xcentre]*np.pi/180.)
            for box_y in xrange(box_length_y):
                distance_y[box_y]=np.abs(box_ystart+box_y-box_ycentre)*model_dict['dy']/float(model_dict['dx'])
            for box_x in xrange(box_length_x):
                for box_y in xrange(box_length_y):
                    distance=np.int(round(np.sqrt(distance_x[box_x]*distance_x[box_x]+distance_y[box_y]*distance_y[box_y])))
                    km_distance=np.sqrt(distance_x[box_x]*distance_x[box_x]+distance_y[box_y]*distance_y[box_y])*model_dict['dx']
                    distance_lists[distance,npts[distance]]=km_distance
                    remote_precip=precip[:,box_y+box_ystart,box_x+box_xstart]
                    if (box_x + box_xstart == box_xcentre) and (box_y + box_ystart == box_ycentre):
                        autocorr=autocorr+1
                    else:
                        distance_correlations[distance]=distance_correlations[distance]+np.corrcoef([central_precip.data,remote_precip.data])[1,0]
                    npts[distance]=npts[distance]+1
            nboxes = nboxes+1
            if nboxes >= max_boxes :
                raise Exception('ERROR: Number of sub-boxes ('+str(nboxes)+') exceeds maximum number of sub-boxes ('+str(max_boxes)+') exceeded.  Increase size of sub-boxes (-b option) or increase parameter max_boxes in code.')
    for dist in xrange(max_box_distance):
        if npts[dist] > 0 :
            distance_correlations[dist]=distance_correlations[dist]/npts[dist]
            distance_ranges[0,dist]=np.amin(distance_lists[dist,0:npts[dist]])
            distance_ranges[1,dist]=np.median(distance_lists[dist,0:npts[dist]])
            distance_ranges[2,dist]=np.amax(distance_lists[dist,0:npts[dist]])
            distance_max=dist
        else:
            distance_correlations[dist]=-999
    print '---> Info: There are '+str(nboxes)+' sub-boxes in your input data.'
    return(distance_correlations,distance_ranges,distance_max)

def compute_autocorr(precip,model_dict):
    
    """
    Compute the lagged auto-correlation of precipitatio
    across all points in the analysis domain.
    
    Arguments:
    * precip:
        An iris cube of precipitation to analyse
    * model_dict:
        The dictionary of information about this dataset
    
    Returns:
    * time_correlations (max_timesteps):
        Composite lagged auto-correlations across all points
    
    
    * time_max (max_timesteps):
        The maximum valid lag in the time_correlation array (can be
        used as a subscript).
    """
    
    print '---> Computing auto-correlations'
    nlon=len(precip.coord('longitude').points)
    nlat=len(precip.coord('latitude').points)
    autocorr_length = model_dict['autocorr_length']
    max_box_distance,max_boxes,max_timesteps = parameters()
    # +1 to account for lag-zero correlation
    autocorr_nt = np.int(autocorr_length//(model_dict['dt']))+1
    time_max = autocorr_nt
    time_correlations = np.zeros(max_timesteps)
    print '----> Info: Computing auto-correlations for '+str(autocorr_nt)+' lags.'
    if autocorr_nt > max_timesteps:
        raise Exception('Error: Number of lags for auto-correlation exceeds maximum ('+str(max_timesteps)+').  Increase parameter max_timesteps in code or reduce autocorrelation length.')
    
    for lon in xrange(nlon):
        for lat in xrange(nlat):
            for lag in xrange(autocorr_nt):
                this_precip = precip.data[:,lat,lon]
                time_correlations[lag] = time_correlations[lag] + np.corrcoef(this_precip,np.roll(this_precip,lag,0))[1,0]

    time_correlations = time_correlations / (nlon*nlat)
    return(time_correlations,time_max)

def compute_spacetime_summary(precip,ndivs):
    
    """
    Computes summary metrics of spatial and temporal coherence,
    as in Klingaman et al. (2017).
        
    Method:
    Precipitation data are binned into "ndivs" divisions (at each gridpoint).
        
    The temporal coherence metric measures the relative frequency of
    persistent upper- and lower-division precipitation to the relative
    frequency of intermittent precipitation (upper-division then
    lower-division, or lower-division then upper-division) on consecutive
    timesteps at the same gridpoint.
        
    The spatial coherence metric measures the relative frequency of
    persistent upper- and lower-quartile precipitation to the relative
    frequency of intermittent precipitation, using data at neighboring
    gridpoints.
        
    In Klingaman et al. (2017), the divisions are quartiles, but this
    can be adjusted with the "ndivs" argument (see below).
        
    Positive values of either measure indicate coherent precipitation.
    Negative values of either measure indicate intermittent precipitation.
        
    The function prints the individual upper- and lower-quartile metrics,
    as well as the combined metric (the ratio).  These values are only
    printed to the screen; they are not plotted.
        
    Arguments:
    * precip:
        An iris cube of precipitation to analyse
    * ndivs:
        The number of divisions for the spatial and coherence metrics.
        Example: ndivs=4 computes the metrics based on upper-quartile and
        lower-quartile precipitation.
            
    Returns:
    * space_inter:
        The combined spatial coherence metric.
    * time_inter:
        The combined temporal coherence metric.
    """
    
    print '---> Computing summary statistics for spatial and temporal intermittency'
    nlon=precip.shape[2]
    nlat=precip.shape[1]
    
    lower_thresh = np.empty((nlat,nlon))
    upper_thresh = np.empty((nlat,nlon))
    for lon in xrange(nlon):
        for lat in xrange(nlat):
            this_precip = precip.data[:,lat,lon]
            this_precip = this_precip[np.where(this_precip > 1)]
            nt = np.size(this_precip)
            if nt > ndivs:
                precip_sorted = np.sort(this_precip)
                lower_thresh[lat,lon] = precip_sorted[np.int(np.floor(nt/ndivs))]
                upper_thresh[lat,lon] = precip_sorted[np.int(np.floor(nt*(ndivs-1)/float(ndivs)))]
            else:
                lower_thresh[lat,lon] = 0
                upper_thresh[lat,lon] = 0

    onon_count=0
    onoff_count=0
    offon_count=0
    offoff_count=0
    non=0
    noff=0
    for lon in xrange(nlon):
        for lat in xrange(nlat):
            this_precip = precip.data[:,lat,lon]
            for t in xrange(nt-1):
                if this_precip[t] > 1:
                    if (this_precip[t] > upper_thresh[lat,lon]):
                        non=non+1
                        if (this_precip[t+1] < lower_thresh[lat,lon]):
                            onoff_count=onoff_count+1
                        elif (this_precip[t+1] > upper_thresh[lat,lon]):
                            onon_count=onon_count+1
                    elif (this_precip[t] < lower_thresh[lat,lon]):
                        noff=noff+1
                        if (this_precip[t+1] < lower_thresh[lat,lon]):
                            offoff_count=offoff_count+1
                        elif (this_precip[t+1] > upper_thresh[lat,lon]):
                            offon_count=offon_count+1

    onon_count = onon_count/float(non)
    offoff_count = offoff_count/float(noff)
    onoff_count = onoff_count/float(non)
    offon_count = offon_count/float(noff)
    time_inter = 0.5*((onon_count+offoff_count)-(onoff_count+offon_count))
    print '-----> Info: Temporal intermittency measure p(upper|upper): ',onon_count
    print '-----> Info: Temporal intermittency measure p(lower|lower): ',offoff_count
    print '-----> Info: Temporal intermittency measure p(upper|lower): ',offon_count
    print '-----> Info: Temporal intermittency measure p(lower|upper): ',onoff_count
    print '----> Info: Combined temporal intermittency measure: ',time_inter
    
    onon_count=0
    onoff_count=0
    offon_count=0
    offoff_count=0
    non=0
    noff=0
    for lon in xrange(1,nlon-1,3):
        for lat in xrange(1,nlat-1,3):
            for t in xrange(nt-1):
                if (precip.data[t,lat,lon] > upper_thresh[lat,lon]):
                    onon_count = onon_count + np.sum(np.where(precip.data[t,lat-1:lat+2,lon-1:lon+2] > upper_thresh[lat,lon],1,0))
                    onoff_count = onoff_count + np.sum(np.where(precip.data[t,lat-1:lat+2,lon-1:lon+2] < lower_thresh[lat,lon],1,0))
                    non=non+1
                if (precip.data[t,lat,lon] < lower_thresh[lat,lon] and precip.data[t,lat,lon] > 1):
                    offoff_count = offoff_count + np.sum(np.where(precip.data[t,lat-1:lat+2,lon-1:lon+2] < lower_thresh[lat,lon],1,0))
                    offon_count = offon_count + np.sum(np.where(precip.data[t,lat-1:lat+2,lon-1:lon+2] > upper_thresh[lat,lon],1,0))
                    noff=noff+1

    onon_count = onon_count / float(non*8.0)
    onoff_count = onoff_count / float(non*8.0)
    offoff_count = offoff_count / float(noff*8.0)
    offon_count = offon_count / float(noff*8.0)
    space_inter = 0.5*((onon_count+offoff_count)-(onoff_count+offon_count))
    print '-----> Info: Spatial intermittency measure p(upper|upper): ',onon_count
    print '-----> Info: Spatial intermittency measure p(lower|lower): ',offoff_count
    print '-----> Info: Spatial intermittency measure p(upper|lower): ',offon_count
    print '-----> Info: Spatial intermittency measure p(lower|upper): ',onoff_count
    print '----> Info: Combined spatial intermittency measure: ',space_inter
    
    return (space_inter,time_inter)

def plot_equalarea_corr(distance_correlations,distance_ranges,distance_max,model_dict=None,colors=None,legend_names=None,set_desc=None):

    """
    Plots correlations as a function of physical distance from one or several datasets,
    using correlation data from compute_equalarea_corr.  The output is a line graph.
    See Fig. 3a in Klingaman et al. (2017) for an example.
    
    Note that the correlation at the central point (0 km) is not plotted, as this is 1.0 by definition.
       
    Arguments:
    * distance_correlations (n_datasets,max_box_distance) or (max_box_distance):
        Composite correlations as a function of physical distance, averaged over
        all sub-regions, as output from compute_equalarea_corr.  If a 2D array, then
        the routine assumes that the input contains > 1 sets of composite correlations
        from multiple datasets.
    * distance_ranges (n_datasets,3,max_box_distance) or (3,max_distance) :
        For each bin of physical distance, the minimum, median and maximum distance
        of points that fall into that bin, as output from compute_equalarea_corr.  If
        a 3D array, then the routine assumes that the input contains > 1 sets of
        range values from multiple datasets.
    * distance_max (n_datasets) or scalar:
        The furthest distance bin for which the data in distance_corelations and
        distance_ranges is valid, as output from compute_equalarea_corr.  If a 1D
        array, then the routine assumes that the input contains > 1 set of values
        from multiple datasets.
            
    Arguments that may be required (see below):
    * model_dict:
        The dictionary containing information about this dataset.  Required only if plotting
        data for one dataset.
    * colors (n_datasets):
        A list of line colors for each dataset.  Required only if plotting data for more than
        one dataset.
    * legend_names (n_datasets):
        A list of legend names for each dataset.  Required only if plotting data for more than
        one dataset.
    * set_desc:
        A string containing a description for this set of datasets.  Used in output plot filename.
        Required only if plotting data for more than one dataset.
    """

    print '--> Plotting correlations vs. distance for all models.'
    max_box_distance,max_boxes,max_timesteps = parameters()
    
    if distance_correlations.ndim == 1:
        if model_dict == None:
            raise Exception('You are plotting correlations for only one dataset, but you have not specified a dataset dictionary with the model_dict option to plot_equalarea_corr.')
        else:
            colors=model_dict['color']
            legend_names=model_dict['legend_names']
            set_desc=model_dict['name']
            nmodels=1
    elif distance_correlations.ndim == 2:
        nmodels=distance_correlations.shape[0]
        if colors == None:
            raise Exception('You are plotting correlations for more than one dataset, but you have not specified a list of plot colors with the colors option to plot_equalarea_corr.')
        if legend_names == None:
            raise Exception('You are plotting correlations for more than one dataset, but you have not specified a list of legend names with the legend_names option to plot_equalarea_corr.')
        if set_desc == None:
            raise Exception('You are plotting correlations for more than one dataset, but you have not specified a description for this dataset with the set_desc option to plot_equalarea_corr.')

    cfp.setvars(file='asop_coherence.'+set_desc+'_precip_spatial_correlations.ps',text_fontsize=20,axis_label_fontsize=20,legend_text_size=18)
    cfp.gopen(figsize=[10,9])

    if nmodels > 1:
        dmax=np.amax(distance_ranges[:,2,:])
        xmax=dmax*1.05
        cfp.gset(xmin=0,xmax=xmax,ymin=-0.1,ymax=1.0)
        for model in xrange(nmodels):
            xpts=distance_ranges[model,1,1:distance_max[model]].flatten()
            ypts=distance_correlations[model,1:distance_max[model]].flatten()
            cfp.lineplot(x=xpts,y=ypts,linestyle=':',marker='o',color=colors[model],markersize=8,label=legend_names[model],
                         xticks=np.arange(0,dmax+1,dmax/10),yticks=np.arange(12)*0.1-0.1)
            for dist in xrange(1,distance_max[model]):
                xpts=[distance_ranges[model,0,dist],distance_ranges[model,2,dist]]
                ypts=[distance_correlations[model,dist],distance_correlations[model,dist]]
                cfp.plotvars.plot.plot(xpts,ypts,linewidth=2,color=colors[model])
            cfp.lineplot(x=xpts,y=ypts,linestyle='None',color=colors[model],markersize=8,legend_location='upper right',
                         xticks=np.arange(0,dmax+1,dmax/10),yticks=np.arange(12)*0.1-0.1)
    else:
        dmax=np.amax(distance_ranges[2,:])
        xmax=dmax*1.05
        cfp.gset(xmin=0,xmax=xmax,ymin=-0.1,ymax=1.0)
        xpts=distance_ranges[1,1:distance_max].flatten()
        ypts=distance_correlations[1:distance_max].flatten()
        cfp.lineplot(x=xpts,y=ypts,linestyle=':',marker='o',color=colors,markersize=8,label=legend_names,legend_location='upper right',
                     xticks=np.arange(0,dmax+1,dmax/10),yticks=np.arange(12)*0.1-0.1)
        for dist in xrange(1,distance_max):
            xpts=[distance_ranges[0,dist],distance_ranges[2,dist]]
            ypts=[distance_correlations[dist],distance_correlations[dist]]
            cfp.plotvars.plot.plot(xpts,ypts,linewidth=2,color=colors)

    cfp.plotvars.plot.plot([0,xmax],[0,0],linestyle=':',color='black')
#    cfp.plotvars.plot.set_xticks(np.arange(0,xmax,max_box_size/10))
    cfp.plotvars.plot.set_xticklabels(np.arange(0,dmax+1,dmax/10,dtype=np.int),fontsize=18)
#    cfp.plotvars.plot.set_yticks(np.arange(12)*0.1-0.1)
    cfp.plotvars.plot.set_yticklabels(np.arange(12)*0.1-0.1,fontsize=18)
    cfp.plotvars.plot.set_xlabel('Distance from central gridpoint (km)',fontsize=20)
    cfp.plotvars.plot.set_ylabel('Lag=0 correlation (mean of sub-regions)',fontsize=20)
    cfp.gclose()

def plot_autocorr(time_correlations,time_max,dt=None,model_dict=None,colors=None,legend_names=None,set_desc=None,legend=True):

    """
    Plots correlations as a function of physical time from one or several datasets,
    using lagged auto-correlation data from compute_autocorr.  The output is a line graph.
    See Fig. 3b in Klingaman et al. (2017) for an example.

    Note that the lag-0 correlation is not plotted, as this is 1.0 by definition.

    Arguments:
    * time_correlations (n_datasets,max_timesteps) or (max_timesteps):
        Composite correlations as a function of physical time, averaged over all points
        in the analysis region, as output from compute_autocorr.  If a 2D array, then
        the routine assumes that the input contains > 1 sets of composite correlations
        from multiple datasets.
    * time_max (n_datasets) or scalar:
        The longest lag for which the data is time_correlations is valid, as output from
        compute_autocorr.  If a 1D array, then the routine assumes that the input contains
        > 1 sets of values from multiple datasets.
        
    Arguments that may be required (see below):
    * model_dict:
        The dictionary containing information about this dataset.  Required only if plotting
        data for one dataset.
    * dt (n_datasets):
        An array containing the temporal sampling frequency for each input dataset.  Required 
        only if plotting data for more than one dataset.
    * colors (n_datasets):
        A list of line colors for each dataset.  Required only if plotting data for more than
        one dataset.
    * legend_names (n_datasets):
        A list of legend names for each dataset.  Required only if plotting data for more than
        one dataset.
    * set_desc:
        A string containing a description for this set of datasets.  Used in output plot filename.
        Required only if plotting data for more than one dataset.
        
    Optional arguments:
    * legend:
        If set to True, include a legend on the graph.  Default is True.
    """

    print '--> Plotting correlations vs. time for all models.'
    
    if time_correlations.ndim == 1:
        if model_dict == None:
            raise Exception('You are plotting correlations for only one dataset, but you have not specified a dataset dictionary with the model_dict option to plot_autocorr.')
        else:
            colors=model_dict['color']
            if legend == True:
                legend_names=model_dict['legend_names']
            set_desc=model_dict['name']
            dt = model_dict['dt']
            nmodels=1
    elif time_correlations.ndim == 2:
        nmodels=time_correlations.shape[0]
        if colors == None:
            raise Exception('You are plotting correlations for more than one dataset, but you have not specified a list of plot colors with the colors option to plot_autocorr.')
        if legend_names == None and legend == True:
            raise Exception('You are plotting correlations for more than one dataset, but you have not specified a list of legend names with the legend_names option to plot_autocorr.')
        if set_desc == None:
            raise Exception('You are plotting correlations for more than one dataset, but you have not specified a description for this dataset with the set_desc option to plot_autocorr.')
        if dt == None:
            raise Exception('You are plotting correlations for more than one dataset, but you have not specified an array of temporal sampling intervals with the dt option to plot_autocorr.')

    cfp.setvars(file='asop_coherence.'+set_desc+'_precip_temporal_correlations.ps',text_fontsize=20,axis_label_fontsize=20,legend_text_size=18)
    cfp.gopen(figsize=[10,9])

    dt_min = dt/60
    tmax=np.amax((time_max-1)*dt_min)
    xmax=tmax+np.amax(dt_min)*0.5
    cfp.gset(xmin=0,xmax=xmax,ymin=-0.1,ymax=1.0)

    if nmodels > 1:
        for model in xrange(nmodels):
            xpts=(np.arange(time_max[model]))*dt_min[model]
            ypts=time_correlations[model,0:time_max[model]]
            print xpts,ypts
            cfp.lineplot(x=xpts[1:],y=ypts[1:],linestyle=':',marker='o',color=colors[model],markersize=8,label=legend_names[model],
                         xticks=np.arange(11)*tmax//10,yticks=np.arange(12)*0.1-0.1)
    else:
        xpts=(np.arange(time_max))*dt_min
        ypts=time_correlations[0:time_max]
        cfp.lineplot(x=xpts[1:],y=ypts[1:],linestyle=':',marker='o',color=colors,markersize=8,label=legend_names,
                     xticks=np.arange(11)*tmax//10,yticks=np.arange(12)*0.1-0.1)

    if legend == True:
        cfp.lineplot(x=xpts,y=ypts,linestyle='None',color=colors[model],markersize=8,legend_location='lower right',xticks=np.arange(11)*tmax//10,yticks=np.arange(12)*0.1-0.1)
    else:
        cfp.lineplot(x=xpts,y=ypts,linestyle='None',color=colors[model],markersize=8,xticks=np.arange(11)*tmax//10,yticks=np.arange(12)*0.1-0.1)

    cfp.plotvars.plot.plot([0,xmax],[0,0],linestyle=':',color='black')
    #    cfp.plotvars.plot.set_xticks(np.arange(11)*xmax//10)
    #    cfp.plotvars.plot.set_yticks(np.arange(12)*0.1-0.1)
    cfp.plotvars.plot.set_xticklabels(np.arange(11)*tmax//10,fontsize=18)
    cfp.plotvars.plot.set_yticklabels(np.arange(12)*0.1-0.1,fontsize=18)
    cfp.plotvars.plot.set_xlabel('Time (minutes)',fontsize=20)
    cfp.plotvars.plot.set_ylabel('Auto-correlation (mean of all points)',fontsize=20)
    cfp.gclose()


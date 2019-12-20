"""
Analysing Scales of Precipitation diagnostics
Precipitation-Duration 2-D Histogram package
These functions allow the user to compute and plot diagnostics of the
histogram for preciptiaiton duration in observations and/or
models.  These functions should work for any resolution or timescale.
These functions assume that the input precipitation dataset is in netCDF.
Written by Amulya Chevuturi
a.chevuturi@reading.ac.uk
(C) The author 2020
"""

from __future__ import division
import iris
import cf, cfplot as cfp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


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
    try:
        precip.coord('latitude').guess_bounds()
    except:
        pass
    try:
        precip.coord('longitude').guess_bounds()
    except:
        pass
    return(precip)
    
def compute_twod_histogram(pr,pbins,cbins):

    """
    Computes 2D histograms of precipitation from an input iris cube of precipitation.
    The 2D histogram is the histogram of precipitation for differnt duration bins.
    It tells you frequency of a precipitaiton bin lasting for a particular duration bin
    Arguments:
    * precip :
       An numpy array of precipitation (TIMExLATxLON)
    * pbins:
       A python list of the precipitation bins for which to compute the histogram,
       in the same units as "precip" above.
    *  cbins:
       A python list of the duration bins for which to compute the histogram,
       in the same units as time in the "precip" above.
    Returns:
    * twod_hist:
       The two-dimensional histogram, as a numpy array, normalised by the size of
       the dataset (TIMExLATXLON) to compare across datasets.        
    """

    print '---> Computing 2D histogram'
    histogram = np.zeros((pr[0,:,0].size*pr[0,0,:].size,len(pbins)-1, len(cbins)-1))
    p = pr.reshape((pr[:,0,0].size, pr[0,:,0].size*pr[0,0,:].size))
    for l in range(pr[0,:,0].size*pr[0,0,:].size):
        for i in range(len(pbins)-1):
	    if i == 0:
	        condition = (p[:,l]>=pbins[i]) & (p[:,l]<=pbins[i+1])
	    else:
	        condition = (p[:,l]>pbins[i]) & (p[:,l]<=pbins[i+1])
            count = np.diff(np.where(np.concatenate(([condition[0]],condition[:-1] != condition[1:],[True])))[0])[::2]
	    histogram[l,i,:] = np.histogram(count, cbins)[0]
    hist_twod = (np.nansum(histogram, axis=0))/(pr.size)

    return hist_twod


def plot_twod_histogram(hist_twod,model_dict,pbins,cbins,title=True,colorbar=True):

    """
    Creates a PostScript plot of the 2D histograms calculated in compute_duration_histogram.
    Arguments:
    * hist_twod:
       The 2D histogram computed in compute_histogram.
    * model_dict:
       The dictioary containing details of this model.
    * pbins:
       A python list of the precipitation bins for which to compute the histogram,
       in the same units as "precip".
    *  cbins:
       A python list of the duration bins for which to compute the histogram,
       in the same units as time in the "precip".
    Optional arguments:
    * title:
       A logical to control whether the title is printed at the top of the plot (useful for creating multi-panel plots).
    * colorbar:
       A logical to control whether the colorbar is printed at the bottom of the plot (useful for creating multi-panel plots).
    """

    print '---> Plotting 2D histogram'
    y = np.arange(len(pbins))
    yt = list(pbins); yt[0] = '<'+str(yt[1]); yt[-1] = '>'+str(yt[-2])
    x = np.arange(len(cbins))
    xt = list(cbins); xt[-1] = '>'+str(xt[-2])
    cmap = plt.get_cmap("viridis_r")
    hist_con_levs=[1e-5,2e-5,4e-5,7e-5,1e-4,2e-4,4e-4,7e-4,1e-3,2e-3,4e-3,7e-3,1e-2,2e-2,4e-2,7e-2,1e-1]
    norm = mpl.colors.BoundaryNorm(hist_con_levs,ncolors=cmap.N,clip=True)
    
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    im = plt.pcolormesh(hist_twod, cmap=cmap, norm=norm)
    plt.xticks(x, xt, rotation=90)
    plt.xlabel('Duration (Days)')
    plt.yticks(y,yt)
    plt.ylabel('Precipitation Bins (mm day$^{-1}$)')
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_aspect((x1-x0)/(y1-y0))
    if title == True:
        plt.title(model_dict['legend_name']+' '+model_dict['time_desc']+' '+model_dict['grid_desc'])
    if colorbar == True:
        cbar = fig.add_axes([0.88, 0.15, 0.03, 0.77])
        fig.colorbar(im, cax=cbar, orientation='vertical', extend='both', format='%.0e')
    fig.tight_layout()
    #plt.show()
    fig.savefig('asop_2D-hist_'+model_dict['legend_name']+'_'+model_dict['time_desc']+'_'+model_dict['grid_desc']+'.ps')
    plt.close(fig)
    

def compute_oned_histogram(pr,tbins,cbins):

    """
    Computes 1D histograms of precipitation from an input iris cube of precipitation.
    The 1D histogram is the histogram of precipitation for differnt duration bins.
    It tells you frequency of a precipitaiton bin lasting for a particular duration bin
    Arguments:
    * precip :
       An numpy array of precipitation (TIMExLATxLON)
    * tbins:
       A python list of ONLY two precipitation bins for which to compute the histogram,
       in the same units as "precip" above. Two precipitation bins are upper and lower
       THRESHOLDS of precipitation bins needed. 
    *  cbins:
       A python list of the duration bins for which to compute the histogram,
       in the same units as time in the "precip" above.
    Returns:
    * oned_hist:
       The two-dimensional histogram, as a numpy array, normalised by the size of
       the dataset (TIMExLATXLON) to compare across datasets.        
    """

    print '---> Computing 1D histogram'
    histogram = np.zeros((pr[0,:,0].size*pr[0,0,:].size,len(cbins)-1))
    p = pr.reshape((pr[:,0,0].size, pr[0,:,0].size*pr[0,0,:].size))
    for l in range(pr[0,:,0].size*pr[0,0,:].size):
        condition = (p[:,l]>=tbins[0]) & (p[:,l]<=tbins[1])
	count = np.diff(np.where(np.concatenate(([condition[0]],condition[:-1] != condition[1:],[True])))[0])[::2]
        histogram[l,:] = np.histogram(count, cbins)[0]
    hist_oned = (np.nansum(histogram, axis=0))/(pr.size)
    return hist_oned

def plot_oned_histogram(hist_oned,model_dict,tbins,cbins,colors,legends,title=True,legend=True):

    """
    Creates a PostScript plot of the 2D histograms calculated in compute_duration_histogram.
    Arguments:
    * hist_oned:
       The 2D histogram computed in compute_histogram.
    * model_dict:
       The dictioary containing details of this model.
    * pbins:
       A python list of the precipitation bins for which to compute the histogram,
       in the same units as "precip".
    *  cbins:
       A python list of the duration bins for which to compute the histogram,
       in the same units as time in the "precip".
    Optional arguments:
    * title:
       A logical to control whether the title is printed at the top of the plot (useful for creating multi-panel plots).
    * colorbar:
       A logical to control whether the colorbar is printed at the bottom of the plot (useful for creating multi-panel plots).
    """

    print '---> Plotting 1D histogram'
    x = np.arange(len(cbins))
    xt = list(cbins); xt[-1] = '>'+str(xt[-2])
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    #plt.axhline(0, color='0.5')
    for i in range(len(hist_oned)):
        plt.plot(x[:-1]+0.5,hist_oned[i], color=colors[i], label=legends[i], marker='.', ls=':', markersize=10, lw=1)
    plt.xticks(x, xt, rotation=90)
    plt.xlabel('Duration (Days)')
    #plt.yticks(np.arange(0,0.1,0.02))
    plt.ylabel('Normalized Frequency')
    #plt.ylim(-0.001,0.08)
    if title == True:
      plt.title('1D Histogram for Precipitation between '+str(tbins[0])+'-'+str(tbins[1])+'mm/day')
    if legend == True:
      plt.legend(loc=1, prop={'size': 15})
    fig.tight_layout()
    #plt.show()
    fig.savefig('asop_1D-hist_'+model_dict['time_desc']+'_'+model_dict['grid_desc']+'.ps')
    plt.close(fig)
    

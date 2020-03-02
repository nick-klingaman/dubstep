"""
Functions to read in and plot precip rain rate histograms as 1D 
line graphs.

Should work for any timescale. Assumes values are in mm/day.

NOTE:
Designed for SMALL or UNIFORM regions, or single gridpoints, as takes spatial
averages of the histograms calculated at each point. This is in order to retain
as much information as possible about the TIME variation of intensity, rather than
the spatial variation.

Written by Gill Martin

(C) British Crown Copyright 2016
"""

import os.path
import sys
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np

from make_hist_maps import extract_region
from make_hist_maps import calc_rain_contr


def guessbounds(cube):
    """
    Guess bounds if cube doesn't have any
    """

    lon_coords = [coord for coord in cube.coords() if "longitude" in coord.name()]
    lat_coords = [coord for coord in cube.coords() if "latitude" in coord.name()]
    for c in [lon_coords[0].name(),lat_coords[0].name()]:
	if not cube.coord(c).has_bounds():
	    cube.coord(c).guess_bounds()
	

def calc_1d_contr(cube,region=None):
    """
    Function to output 1d histogram of contributions
    
    Args:
    * cube:
       histogram cube of contributions from each bin
    * region:
       region to calculate histogram over. NOTE that histograms for each gridbox
       in the region are averaged together (the histogram is NOT re-calculated for
       the region from the original data). Therefore this should be a small and/or 
       uniform (similar histograms in each gridbox) region.
    """

# Extract a region if necessary, or plot whole region

    if region:
        cube_exreg = extract_region(cube,region)
    else:
        cube_exreg = cube
    
    guessbounds(cube_exreg)
    grid_areas = iris.analysis.cartography.area_weights(cube_exreg)
    hist1d=cube_exreg.collapsed(['latitude','longitude'], iris.analysis.MEAN,weights=grid_areas)

    return hist1d
    
    
def plot_1dhist(plotname,region,filenames,runtitles,plottitle,timescale=None,filenames_obs=None,runtitles_obs=None,frac=None,col_offset=None,log=None):

    """
    Function to plot 1d histograms. Can overplot several (e.g. at different timescale or resolution, 
    or plot differences between one or more pairs.
    Also allows histograms from a different type of dataset (here it's from observations, assuming the
    default is models) to be overplotted in the same colours but with dashed lines, for comparison.
    NOTE that plotting differences excludes overplotting alternative datasets.
    
    Args:
    * plotname:
       filename for resulting plot.
    * region:
       [W,S,E,N] limits of region to average over
    * filenames:
       list of filenames for files of histograms (output from make_hist_maps)
       OR list of PAIRS OF filenames for differencing (b-a) in which case filenames_obs etc should be OMITTED.
    * runtitles:
       list of runtitles for legend (must have title for all filenames)
       OR list of PAIRS OF runtitles for differencing, matching the filename pairs (b-a).
    * plottitle:
       main title for plot

    Optional arguments:
    * timescale:
       If set, adds the timescale to the plot subtitle
    * filenames_obs:
       list of obs filenames for files of histograms (output from make_hist_maps)
    * runtitles_obs:
       list of obs runtitles (must have title for all obs filenames)
    * frac:
       If frac is not None, calculates fractional contribution, rather than absolute
    * log:
       If log is not None, uses log scale for X axis only
    * col_offset:
       If set, offsets start of colour set by whatever this is set to
    """

# Check that a 4-element region has been entered, if required

    if region is None or len(region) != 4:
        raise Exception('Please enter a 4-element region!'+str(region))

# Check that all filenames have a runtitle

    if len(runtitles) != len(filenames):
        raise Exception('Not all run filenames have a runtitle!')
    if isinstance(filenames[0], list):
        if (not isinstance(runtitles[0], list)) or (len(runtitles[0]) != len(filenames[0])):
            raise Exception('Not all run filenames have a runtitle!')

    if filenames_obs is not None:
        if runtitles_obs is None:
	    raise Exception('Please supply runtitles for obs filenames!')
	elif len(runtitles_obs) != len(filenames_obs):
            raise Exception('Not all obs filenames have a runtitle!')

# Go through filenames reading in cube_of_hist from file and calculating regional averages of the
# histograms of contributions to average rainfall

    nice_cmap = plt.get_cmap('brewer_Paired_12')
    num_colors=12
    color=iter(nice_cmap(np.linspace(0,1,num_colors)))
    if col_offset is not None:
        for j in range(0,col_offset):
	    c=next(color)    
    
    for i, filename in enumerate(filenames):
        c=next(color)
	
        if isinstance(filename, str):
            cube_of_hist=iris.load_cube(filename)
            avg_rain_bins, avg_rain_bins_frac = calc_rain_contr(cube_of_hist)
            if frac is not None:
                hist1d=calc_1d_contr(avg_rain_bins_frac,region)
	    else:
                hist1d=calc_1d_contr(avg_rain_bins,region)
	
        else:
            cube_of_hist_a=iris.load_cube(filename[0])
            cube_of_hist_b=iris.load_cube(filename[1])
            avg_rain_bins_a, avg_rain_bins_frac_a = calc_rain_contr(cube_of_hist_a)    
            avg_rain_bins_b, avg_rain_bins_frac_b = calc_rain_contr(cube_of_hist_b)
            if frac is not None:
                hist1da=calc_1d_contr(avg_rain_bins_frac_a,region)
                hist1db=calc_1d_contr(avg_rain_bins_frac_b,region)
	    else:
                hist1da=calc_1d_contr(avg_rain_bins_a,region)
                hist1db=calc_1d_contr(avg_rain_bins_b,region)
	
            hist1d=hist1db-hist1da

        if type(runtitles[i]) is list:
            iplt.plot(hist1d[1:90],label='{} minus {}'.format(runtitles[i][1], runtitles[i][0]),color=c,linewidth=2.5)
        else:
            iplt.plot(hist1d[1:90],label=runtitles[i],color=c,linewidth=2.5)

# Now add plot titles and legend

    a0=str(int(round(region[0])))
    a2=str(int(round(region[2])))
    a1=str(int(round(region[1])))
    a3=str(int(round(region[3])))
    
    if frac is not None:
        if type(filenames[0]) is list:
	    method = 'DIFFERENCE in Frac'
	else:
            method = 'Fractional'
	units = ''
    else:
        if isinstance(filenames[0], list):
	    method = 'DIFFERENCE in Avg'
	else:
            method = 'Average'
	units = 'mm/day'
	    
    if timescale is None:
        timescale = ''
	
    plt.title('{} precip contribution from {} events \n {} \n {} - {}E  {} - {}N'.format(method, timescale, plottitle, a0, a2, a1, a3), fontsize=12)
    plt.ylabel('{} Precip contribution {}'.format(method, units), fontsize=10)

# Now add dashed lines for another type of data (e.g. obs) if supplied 
# ONLY if we haven't plotted differences in the above code (indicated by single, not paired, filenames)

    if isinstance(filenames[0], str) and filenames_obs is not None:
        nice_cmap = plt.get_cmap('brewer_Greys_09')
        color=iter(nice_cmap([8,6,4,2]))
        linestyle=iter(['dashed','dotted','solid','dashed'])
        for i, (v, w) in enumerate(zip(filenames_obs,runtitles_obs)):
            c=next(color)
            c2=next(linestyle)
            filename=v
	    runtitleo=w
            cube_of_hist=iris.load_cube(filename)
            avg_rain_bins, avg_rain_bins_frac = calc_rain_contr(cube_of_hist)

            if frac is not None:
                hist1d=calc_1d_contr(avg_rain_bins_frac,region)
            else:
                hist1d=calc_1d_contr(avg_rain_bins,region)
		
            iplt.plot(hist1d[1:90],label=runtitleo,color=c,linewidth=2.5,linestyle=c2)

# Label axes and add legends
    
    plt.xlabel('Precip bin mm/day', fontsize=10)
    if log is not None:
        plt.xscale('log')
        plt.xlim((0.1,1000))
        plt.legend(ncol=2, fontsize=9, loc='upper left')
    else:
        plt.xlim((0.1,400))
        plt.legend(ncol=2, fontsize=9, loc='upper right')

# Save the figure to the supplied plot filename

    plt.savefig(plotname)

    plt.clf()



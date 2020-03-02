"""
Functions to read in and plot precip rain rate histogram maps
Should work for any timescale. Assumes values are in kg/m2/day.

Written by Gill Martin

(C) British Crown Copyright 2016
"""

import os.path
import sys
import numpy as np

import iris
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import from_levels_and_colors
from matplotlib.colors import NoNorm

from make_hist_maps import extract_region


def _make_panel(plot_position, ppn_cube, ppn_bounds, cmap, norm):
    """
    Make a single panel for the plot of precipitation histograms by extracting and averaging
    histogram over bins as defined by ppn_bounds.
    Assumes histogram cube was obtained from precipitation data in kg/m2/day
    """
    
    ax = plt.subplot(plot_position[0],plot_position[1],plot_position[2], projection=ccrs.PlateCarree())
    ax.coastlines()
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = False
    gl.xlabels_bottom = False
    gl.ylabels_right = False
      
    if plot_position[2] >= plot_position[0] * (plot_position[1]) - 1: 
        #True for bottom row of panels 
        gl.xlabels_bottom = True
	 
    gl.xlabel_style = gl.ylabel_style = {'size': 6}   
	# the x and y label styles are now tied to each other  
	# (both are pointers to a single dictionary) 
	 
    if plot_position[2] <= 2: 
        #True for top row of panels 
        title = '%s events %s to %s mm/day' % ( ppn_cube.attributes['timescale'], ppn_bounds[0], int(ppn_bounds[1]))
    elif plot_position[2] >= plot_position[0] * (plot_position[1]) - 1:
        #True for bottom row of panels 
        title = '%s events > %s mm/day' % ( ppn_cube.attributes['timescale'], int(ppn_bounds[0]))
    else:
        title = '%s events %s to %s mm/day' % ( ppn_cube.attributes['timescale'], int(ppn_bounds[0]), int(ppn_bounds[1]))
    
    plt.title(title , fontsize=8)
    
    ppn_con = iris.Constraint(precipitation_flux = lambda cell: ppn_bounds[0] <= cell < ppn_bounds[1])
	 
    field_to_plot = ppn_cube.extract(ppn_con).collapsed('precipitation_flux', iris.analysis.SUM)
    cf = iplt.pcolormesh(field_to_plot, cmap=cmap, norm=norm)   
    
    # return objects that you may want to work on, for example you could include the axes or gridlines objects   
    
    return cf 

     
def plot_rain_contr(cube_a,cube_b,plotname,runtitle,timescale_a,timescale_b,all_ppn_bounds=None,region=None,frac=0):
    """
    Function to plot set of histogram maps as rainfall contributions
    on two different timescales side by side.
    
    Inputs can be ACTUAL (assumed) or FRACTIONAL (if frac=1) contributions
    
    Args:
    * cube_a, cube_b:
       histogram cubes of rainfall contributions from each bin on timescale_a and timescale_b
    * plotname:
       filename for plot
    * runtitle:
       Name for dataset being plotted
    * timescale_a, timescale_b:
       Strings detailing time frequency of data from which histogram cubes were calculated
       e.g. 'Daily'
    * all_ppn_bounds:
       (optional) defines range of bin limits to lump together for plotting, otherwise is set
       to 4 default groups: [(0.005, 10.), (10., 50.), (50., 100.), (100., 3000.)] kg/m2/day
       Note that plots will work best for no more than 5 or 6 groups.
    * region:
       (optional) defines region to plot, otherwise just plots everywhere there are data.
    * frac:
       (optional) If frac=1, expects FRACTIONAL contributions have been input
                  otherwise assumes ACTUAL contributions have been input
    """

# Check that the correct fields have been input

    test_a=cube_a.collapsed('precipitation_flux', iris.analysis.SUM)
    test_b=cube_b.collapsed('precipitation_flux', iris.analysis.SUM)

    if frac != 0:
        if round(np.max(test_a.data),12) != 1.0 or round(np.max(test_b.data),12) != 1.0:
            raise Exception('One or more input cube(s) is not a fractional histogram')
    else:
        if round(np.max(test_a.data),12) == 1.0 or round(np.max(test_b.data),12) == 1.0:
            raise Exception('One or more input cube(s) is a fractional histogram')


# Extract region if required

    if region:
        avg_rain_bins_a = extract_region(cube_a,region)
        if avg_rain_bins_a is None:
	    raise Exception('No data points in first dataset exist within region!')

        avg_rain_bins_b = extract_region(cube_b,region)
        if avg_rain_bins_b is None:
	    raise Exception('No data points in second dataset exist within region!')
    else:
        avg_rain_bins_a = cube_a
	avg_rain_bins_b = cube_b
    
# Put timescale and runtitle data into cube attributes so that they can be extracted later.
   
    avg_rain_bins_a.attributes['timescale'] = timescale_a
    avg_rain_bins_a.attributes['runtitle'] = runtitle
    avg_rain_bins_b.attributes['timescale'] = timescale_b
    avg_rain_bins_b.attributes['runtitle'] = runtitle		

# Set up contour levels and colours
    
    if frac != 0:
        contour_levels = [0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    else:
	contour_levels = [0.25,0.5,0.75,1,1.5,2,2.5,3,4,6,8,10,12]

# Make a range of colours from dark blue through near-white to red by concatenating
# part of two Cynthia Brewer color maps.

    colors = np.concatenate([plt.get_cmap('brewer_RdBu_11'  )(range(1, 10)), 
             plt.get_cmap('brewer_YlOrBr_09')(range(5, 0, -1))])  
    
    cmap, norm = from_levels_and_colors(contour_levels, colors, extend='both')

# Plot contribution amounts, lumped into bins as defined in all_ppn_bounds
				
    if not all_ppn_bounds:
        all_ppn_bounds = [(0.005, 10.), (10., 50.), (50., 100.), (100., 3000.)]

    cubes = [avg_rain_bins_a,avg_rain_bins_b]

    plot_shape = [len(all_ppn_bounds) , len(cubes) ]
    fig = plt.figure(figsize=(10,6), dpi=300)
    
    plot_number = 0
    for ppn_bounds in all_ppn_bounds:
        for cube in cubes:
            plot_number += 1
            cf = _make_panel(tuple(plot_shape + [plot_number]), cube, ppn_bounds, cmap, norm)


# Make an axes to put the shared colorbar in

    colorbar_axes = plt.gcf().add_axes([0.25, 0.05, 0.5, 0.02])
    colorbar = plt.colorbar(cf, colorbar_axes, orientation='horizontal')
    colorbar.ax.tick_params(labelsize=8)
    colorbar.ax.xaxis.set_label_position('top')

    if frac != 0:
        colorbar.set_label('Fractional contribution from events',fontsize=10)
    else:
        colorbar.set_label('Actual contribution from events (mm/day)',fontsize=10)
   
    plt.suptitle('Histogram of total rainfall events (mm/day)' + '\n' + runtitle, fontsize=10)

    plt.savefig(plotname)

    plt.clf()


# (C) British Crown Copyright 2016
"""
Functions to create histograms of precipitation rate from CMIP5 data
Output can be saved as netCDF files using iris.save
Written by Gill Martin
(C) British Crown Copyright 2016
"""

import os.path
import sys
import numpy as np

import iris


def extend_cube_with_dimcoord(cube,dimcoord,dimcoord_index=0):
    """
    Function to create a new cube with an extra dimension as specified
    
    Args:
    * cube:
        input data cube
    * dimcoord:
        dimension coordinates to be added
    * dimcoord_index:
        locates where extra dimension is added
    """
    
    coords = list(cube.dim_coords)
    coords.insert(dimcoord_index,dimcoord)
    
    dim_coords_and_dims = [(c,i) for i,c in enumerate(coords)]
    
    newcube_shape = tuple(c.shape[0] for c in coords)
    
    newcube = iris.cube.Cube(np.zeros(newcube_shape), 
                             dim_coords_and_dims = dim_coords_and_dims,
                             standard_name = cube.standard_name,
                             long_name = cube.long_name,
                             var_name = cube.var_name,
                             attributes = cube.attributes,
                             units = cube.units)
        
    return newcube
    
                            
def hist_count(data, pbin_min, pbin_max):
    """
    Function to calculate where values fall between pbin_min and pbin_max.
    
    Args:
    * data (cube):
       raw data to compare against pbin limits.
       
    * pbin_min (float):
       lower limit of pbin
    
    * pbin_max (float):
       upper limit of pbin
           
    """
    
    if data.coord('time') is None:
        raise Exception('cube does not have a time dimension:'+data)
	
# Count occurrence of values within each pbin, with no overlaps at pbin edges

    return data.collapsed('time', iris.analysis.COUNT, function=lambda values: (pbin_min <= values) & (values < pbin_max))


def read_data_cube(filename):
    """
    Function to read in data from filename
    
    Args:
    * filename:
       name of netCDF file containing single cube of time-varying values from one variable

    """
    
    data_cube = iris.load_cube(filename)
		
    return data_cube 
    
    
def extract_region(cube, region):
    """
    Function to extract data from a region
    Args:
    * cube:
       input cube
    * region:
       region to create histograms within.
       List of [west longitude, south latitude, east longitude, north latitude] in degrees.
    """
    import cartopy.crs as ccrs

    if len(region) != 4:
	raise ValueError('region does not have 4 elements:'+str(region))
    else:
        lon_coords = [coord for coord in cube.coords() if "longitude" in coord.name()]
        lat_coords = [coord for coord in cube.coords() if "latitude" in coord.name()]

	region_use = region
	
        # Check coordinate system and transform region values accordingly

        if cube.coord(lon_coords[0].name()).coord_system:
	    coord_sys = cube.coord(lon_coords[0].name()).coord_system.as_cartopy_crs()
            if type(coord_sys) == ccrs.RotatedGeodetic:
	        req_coords = ccrs.Geodetic()
                region_use[0],region_use[1] = rot_pole.transform_point(region[0], region[1], req_coords)
                region_use[2],region_use[3] = rot_pole.transform_point(region[2], region[3], req_coords)
	        
	long_limits = iris.coords.CoordExtent(lon_coords[0].name(), region_use[0], region_use[2])
	lat_limits = iris.coords.CoordExtent(lat_coords[0].name(), region_use[1], region_use[3])
   
	reg_cube = cube.intersection(long_limits, lat_limits)
	if reg_cube is None:
	    raise ValueError('No data points exist within region!'+str(region))

    return reg_cube


def make_hist_ppn(ppn_cube,region=None):
    """
    Function to create set of histograms for precipitation data input from filename
    NOTE that this function explicitly assumes precipitation data are input.
    Could be altered to work for any variable provided bin definitions are adjusted.
    
    Args:
    * ppn_cube:
       input timeseries of ppn data
    * region:
       region to create histograms within.
       List of [west longitude, south latitude, east longitude, north latitude] in degrees.
       If not specified, use all available gridpoints.

    Internal variables:
    * pbin2, pbin:
       define bin limits.
       The pbin2 distribution function is defined in Klingaman et al. (2016) doi:10.5194/gmd-2016-161
       pbin[0] = 0.0 is lowest bin limit. 
       pbin[1] to pbin[99] = pbin2 range from 0.005 to 2209 kg/m2/day. 
       Counts of values into element [i] of the output histogram are made such that pbin[i] <= values < pbin[i+1]
       Counts of values >= 2209 kg/m2/day are included in the 100th element.  
       This means that histogram will sum to the number of intervals in the sample.

    Returned:
    * cube_of_hist
       cube containing histogram of counts in each of 99 bins 
    """
    
# Ensure precipitation data are in correct units for pbin definitions

    if (ppn_cube.units == '') or (ppn_cube.units is None):
        raise ValueError('Data cube has no units!')    
    else:
        if (str(ppn_cube.units)[0:1] == 'm'):
            ppn_cube.convert_units('mm day-1')
        elif (str(ppn_cube.units)[0:2] == 'kg') or (str(ppn_cube.units)[0:1] == 'g'):
            ppn_cube.convert_units('kg m-2 day-1')
        else:
            raise ValueError('Unrecognised units:',str(ppn_cube.units))
	
# Extract region if required

    if region:
        reg_ppn_cube = extract_region(ppn_cube, region)
    else:
    	reg_ppn_cube = ppn_cube

# Set up 100 bin limits for precipitation data in kg/m2/day, making lowest limit zero (i.e. pbin[0] = 0.0) so that 
# histogram will sum to total no. of intervals in sample.

    pbin2=np.exp(np.log(0.005)+np.sqrt(np.linspace(0,98,99)*((np.square(np.log(120.)-np.log(0.005)))/59.)))
    pbin=np.zeros(100)
    pbin[1:]=pbin2[0:99]

# Make cube of counts for the first bin as template for the other 99

    ecount=hist_count(reg_ppn_cube, pbin[0], pbin[1])
    ecount.var_name='Precipitation'
    ecount.long_name='No of events'

# Use extend_cube_with_dimcoord function to create a new cube with an extra dimension of 100 bins
# And insert the data from ecount into the first data slice.
 
    pbincoord = iris.coords.DimCoord(pbin, standard_name='precipitation_flux', long_name='precipitation_bin', units=ppn_cube.units)
    cube_of_hist = extend_cube_with_dimcoord(ecount, pbincoord)
    cube_of_hist.data[0] = ecount.data   
   
# Loop through the other 98 pbins, replacing the relevant data slices in cube_of_hist.
# Then include the count of any values >= pbin[99] in the 100th bin.

    for x in range(1,99):
        ecount1=hist_count(reg_ppn_cube, pbin[x], pbin[x+1])
        cube_of_hist.data[x] = ecount1.data

    ecount1=reg_ppn_cube.collapsed('time', iris.analysis.COUNT, function=lambda values: pbin[99] <= values)	
    cube_of_hist.data[99] = ecount1.data
	
    return cube_of_hist

   
def calc_rain_contr(cube_of_hist):
    """
    Function to calculate histogram maps as total and fractional precipitation contributions
    for histogram cube of counts in each bin.
    
    Args:
    * cube_of_hist:
       histogram cube of counts
       
    Calls: 
    * Function make_cube to create merged cube of histogram
    """

# Set up 100 bin mid-points for precipitation data in kg/m2/day, based on function used to define bin limits.
# Note that these will be slightly to the right of the geometric centre of the bin.
 
# Define first bin midpoint as 0.0 since this bin includes all events where there was no rain. This effectively
# limits the calculation of rainfall contributions to events >0.005 kg/m2/day.

    bin_mid2=np.exp(np.log(0.005)+np.sqrt(np.linspace(0.5,98.5,99)*((np.square(np.log(120.)-np.log(0.005)))/59.)))
    bin_mid=np.zeros(100)
    bin_mid[1:]=bin_mid2[0:99]

# Reshape bin_mid so that can multiply it with the correct dimension in cube_of_hist

    bin_mid3=bin_mid.reshape(100,1,1)

# Average rain in gridbox over data sample = sum over bins(no. events in bin x avg. ppn rate in bin) / (no. of intervals in sample)
# Note that cube_of_hist is designed such that the sum over the whole histogram = number of intervals in sample.

    tot_intervals = cube_of_hist.collapsed('precipitation_flux', iris.analysis.SUM)
    
    tot_rain_bins = cube_of_hist * bin_mid3
    avg_rain_box = tot_rain_bins.collapsed('precipitation_flux', iris.analysis.SUM) / tot_intervals

# Contribution to average rainfall in gridbox over data sample by events in bin = 
#             (no. of events in bin x avg. ppn rate in bin)/(no. of intervals in sample)

    avg_rain_bins = tot_rain_bins / tot_intervals

# Fractional contribution from bin = contribution from bin / average rain in gridbox

    avg_rain_bins_frac = avg_rain_bins / avg_rain_box

    return(avg_rain_bins,avg_rain_bins_frac)


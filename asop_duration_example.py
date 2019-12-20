import asop_duration as asop
import numpy as np

"""
    Example use of ASoP Duration package to compute
    and plot diagnostics of the histogram for
    preciptiaiton duration in a dataset.
    
    This example uses TRMM 3B42v7A and CMORPH v1.0
    data at a common grid of N96. 
    
    Written by Amulya Chevuturi
    a.chevuturi@reading.ac.uk
    (C) The author 2020
    
"""

def get_dictionary(key):
    
    """
    The Coherence package relies on a dataset dictionary,
    for which currently the user must specify most values.
        
    This function shows how to build a dictionary.  The dictionary to
    be returned to the main program is selected through the use of
    a "key", here either "TRMM" or "CMORPH"
        
    Arguments:
        * key:
            A string used to select the correct dataset
            
    Returns:
        * asop_dict:
            A dictionary containing required and optional
            parameters for the ASoP Coherence package.
    
        Required dictionary keys:
        
        infile       - Path to input file (must be readable by Iris)
        name         - Name for model in plot files (no spaces)
        legend_name  - Name for model in legends and titles on plots (spaces allowed)
        dt           - Timestep of input data (in seconds)
        dx           - Longitudinal grid spacing at equator (in km)
        dy           - Latitudinal grid spacing (in km)
        constraint   - standard_name of data to read from netCDF file (e.g., precipitation flux)
        scale_factor - Multiplier necessary to convert precipitation to units of mm/day
        region       - Region of data to read [minlat, maxlat, minlon, maxlon]
        box_size     - Length of sub-regions (square boxes) to consider for correlation analysis
                        as a function of physical distance
                        (in km, recommended value is > 6.5*dx).
        color        - Name of color to use on line graphs (must be recognised by matplotlib).
        region_size  - Length of sub-regions (square boxes) to consider for correlation analysis
                        as a function of model gridpoints (including lag correlations, see below)
                        (in units of native gridpoints, odd integers strongly recommended).
        lag_length   - Maximum lag to consider for correlation analysis as a function of model
                        gridpoints, for constructing distance vs. lag correlation diagrams.
                        Correlations for lags from 0 (coincidence) until lag_length will be computed.
        autocorr_length - The maximum autocorrelation to analyse (in seconds).
        
        Optional dictionary keys, which are useful mainly if analysing the same
        dataset on more than one grid / temporal sampling interval:
        
        grid_type    - A string describing the grid, used in output filenames
                        (recommend no spaces).
        time_type    - A string describing the temporal sampling, used in output filenames
                        (recommend no spaces).
        grid_desc    - A string describing the grid, used in plot titles (can contain spaces).
        time_desc    - A string describing the temporal sampling, used in plot titles
                        (can contain spaces).
    """
    
    asop_dict = {}
    if key == 'TRMM':
        asop_dict['infile']       = '/gws/nopw/j04/klingaman/amulya/data/obs/TRMM/TRMM_3B42_V7_0.25deg-DAILY_DJF_1998-2018_brazil_n48.nc'
        asop_dict['name']         = 'TRMM'
        asop_dict['dt']           = 10800
        asop_dict['dx']           = 27
        asop_dict['dy']           = 27
        asop_dict['constraint']   = 'precipitation_flux'
        asop_dict['scale_factor'] = 1.0
        asop_dict['legend_name']  = 'TRMM'
        asop_dict['region']       = [-10,5,288,310]
        asop_dict['box_size']     = 1680
        asop_dict['color']        = 'red'
        asop_dict['region_size']  = 7
        asop_dict['lag_length']   = 6
        asop_dict['grid_type']    = 'N48'
        asop_dict['time_type']    = 'Day'
        asop_dict['grid_desc']    = 'N48'
        asop_dict['time_desc']    = 'Day'
        asop_dict['autocorr_length'] = 60*60*24
    elif key == 'CMORPH':
        asop_dict['infile']       = '/gws/nopw/j04/klingaman/amulya/data/obs/CMORPH/CMORPH_V1.0_0.25deg-DAILY_DJF_1998-2017_brazil_n48.nc'
        asop_dict['name']         = 'CMORPH'
        asop_dict['dt']           = 10800
        asop_dict['dx']           = 27
        asop_dict['dy']           = 27
        asop_dict['constraint']   = 'precipitation_flux'
        asop_dict['scale_factor'] = 1.0
        asop_dict['legend_name']  = 'CMORPH'
        asop_dict['region']       = [-10,5,288,310]
        asop_dict['box_size']     = 1680
        asop_dict['color']        = 'blue'
        asop_dict['region_size']  = 7
        asop_dict['lag_length']   = 6
        asop_dict['grid_type']    = 'N48'
        asop_dict['time_type']    = 'Day'
        asop_dict['grid_desc']    = 'N48'
        asop_dict['time_desc']    = 'Day'
        asop_dict['autocorr_length'] = 60*60*24

    return(asop_dict)

if __name__ == '__main__':

    datasets = ('TRMM','CMORPH')
    n_datasets = len(datasets)

    # Allocate memory for multi-model fields
    max_box_distance,max_timesteps,max_boxes = asop.parameters()
    all_distance_correlations = np.zeros((n_datasets,max_box_distance))
    all_distance_ranges = np.zeros((n_datasets,3,max_box_distance))
    all_distance_max = np.zeros((n_datasets),dtype=np.int)
    all_time_correlations = np.zeros((n_datasets,max_timesteps))
    all_time_max = np.zeros((n_datasets),dtype=np.int)
    all_dt = np.zeros((n_datasets),dtype=np.int)
    all_colors = []
    all_legend_names = []
    hist_1d = [0]*n_datasets

    for i in xrange(n_datasets):
    
        print '--> '+datasets[i]
        asop_dict = get_dictionary(datasets[i])
    
        # Read precipitation data
        precip = asop.read_precip(asop_dict).data
	#convert the iris cube to a numpy array
	#precip = prp.data

        # Define the precipitation bins as prbins
	# Define the duration length bins as cobins
        # Note that on plots, the first and last edges will be
        # replaced by < and > signs, respectively.
        prbins=[0,1,2,4,6,9,12,16,20,25,30,40,60,90,130,180,2e20]
        cobins=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,2e20]
	thbins=[0,0.1]

        # Compute 2D precipitation-duration histograms
        hist_2d = asop.compute_twod_histogram(precip,prbins,cobins)

        # Plot 2D histograms
        asop.plot_twod_histogram(hist_2d,asop_dict,prbins,cobins,title=True,colorbar=True)
	
        # Compute 2D precipitation-duration histogram for a threshold bin
	hist_1d[i] = asop.compute_oned_histogram(precip,thbins,cobins)

        # Save color and legend information
        all_colors.append(asop_dict['color'])
        all_legend_names.append(asop_dict['legend_name'])
	
    # Plot 1D histogram
    asop.plot_oned_histogram(hist_1d,asop_dict,thbins,cobins,all_colors,all_legend_names,title=True,legend=True)

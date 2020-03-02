def ASoP_tests():

    import os.path
    import sys
    import numpy as np
    
    import iris
    
    from ASoP1_spectral import make_hist_maps
    from ASoP1_spectral import plot_hist_maps
    from ASoP1_spectral import plot_hist1d
    
    THIS_DIR = os.path.dirname(os.path.abspath(__file__))
    my_data_path = os.path.join(THIS_DIR, 'reference/test_precipitation_dataset.nc')

# Read in test dataset provided, load data and create histograms
 
    filename=my_data_path
    ppn=make_hist_maps.read_data_cube(filename)
    histmap=make_hist_maps.make_hist_ppn(ppn)

# Save the histograms in another netCDF file. This should compare exactly with the
# test_precipitation_hist.nc file provided.

    my_histmap_file = os.path.join(THIS_DIR, 'mytest_histmap.nc')
    iris.save(histmap, my_histmap_file)

# Plot set of maps of contributions to rainfall over the period
# NOTE: Designed to compared two datasets (e.g. on different timescales) side by side.
#       Since only one dataset is used here, plots the SAME data side by side.

    filename_a=my_histmap_file
    filename_b=filename_a
    
    histcube_a=make_hist_maps.read_data_cube(filename_a)
    histcube_b=make_hist_maps.read_data_cube(filename_b)

    runtitle='Artificial rainfall dataset'

# Calculate actual and fractional contributions to rainfall over period

    avg_rain_bins_a,avg_rain_bins_frac_a=make_hist_maps.calc_rain_contr(histcube_a)
    avg_rain_bins_b,avg_rain_bins_frac_b=make_hist_maps.calc_rain_contr(histcube_b)
    
# Check that the fractional contributions sum to 1.0 at all points

    test_a=avg_rain_bins_frac_a.collapsed('precipitation_flux', iris.analysis.SUM)
    test_b=avg_rain_bins_frac_b.collapsed('precipitation_flux', iris.analysis.SUM)

    if round(np.max(test_a.data),12) != 1.0 or round(np.max(test_b.data),12) != 1.0:
        raise Exception('One or more fractional histograms does not sum to 1.0')
   
#------------
# Plotting
#------------

# First plot actual contribution to total rainfall in period. This should compare exactly with the
# test_precipitation_histmap_act.png file provided.

    plotname=os.path.join(THIS_DIR, 'mytest_histmap_act.png')
    plot_hist_maps.plot_rain_contr(avg_rain_bins_a,avg_rain_bins_b,
                                       plotname,runtitle,'hourly','hourly')
    
# Now plot fractional contribution to total rainfall in period. This should compare exactly with the
# test_precipitation_histmap_frac.png file provided.
    
    plotname=os.path.join(THIS_DIR, 'mytest_histmap_frac.png')
    plot_hist_maps.plot_rain_contr(avg_rain_bins_frac_a,avg_rain_bins_frac_b,
                                       plotname,runtitle,'hourly','hourly',frac=1)
    
    
# Plot 1-d histograms of fractional and actual contributions for a small region. 
# These should compare exactly with the test_hist1d_frac(act).png files provided. 

    region=[215.0,-5.0,225.0,5.0]
    filenames=[my_histmap_file]
    runtitles=['Artificial rainfall dataset']
    plottitle='Artificial rainfall dataset for testing'
    timesc='hourly'

    plotname=os.path.join(THIS_DIR, 'mytest_hist1d_frac.png')
    plot_hist1d.plot_1dhist(plotname,region,filenames,runtitles,plottitle,
                                       timesc,frac=1,col_offset=2,log=1)

    plotname=os.path.join(THIS_DIR, 'mytest_hist1d_act.png')
    plot_hist1d.plot_1dhist(plotname,region,filenames,runtitles,plottitle,
                                       timesc,col_offset=2,log=1)
	       
				       
				       

# (C) British Crown Copyright 2016, Met Office.
"""
Main function to show how to use the code package to process
two model precipitation datasets and compare them with each
other and with two observational datasets.

These could be from different data sources or different timescales
or spatial resolutions of the same data source.

Author: Gill Martin
"""

import iris
import argparse
   
import make_hist_maps
import plot_hist_maps
import plot_hist1d

def main(data_filename1,data_filename2,hist_filename1,hist_filename2,obs_hist_filename1,obs_hist_filename2):
    """ 
    Read in data and create histogram cubes, save these to netcdf files.
    Then plot histogram maps and some regional 1d histograms
    
    Arguments:
    * data_filename1,data_filename2
        input filenames of time-varying precipitation data
      
    * hist_filename1,hist_filename2
        chosen filenames for output histogram cubes (netcdf files)

    * obs_hist_filename1,obs_hist_filename2
        filenames for previously-calculated observed histogram cubes (netcdf files)
    """ 

    # Make and save histogram cubes
    making_histogram_files(data_filename1,data_filename2,hist_filename1,hist_filename2)

    # Plot histogram maps

    runtitle='My dataset'
    plotname_root='compare_ppn1_ppn2'
    plot_histogram_maps(hist_filename1,hist_filename2,runtitle,plotname_root)

    # Set up list of filenames and runtitles for 1d histogram plots
    # Please edit plotname_root, plottitle, runtitles and timescale, as well as region, 
    # as required.
    
    filenames=[]
    runtitles=[]
    
    filenames.append(hist_filename1)
    runtitles.append('My first dataset')
    filenames.append(hist_filename2)
    runtitles.append('My second dataset')
    
    filenames_obs=[]
    runtitles_obs=[]

    filenames_obs.append(obs_hist_filename1)
    runtitles_obs.append('First obs dataset')
    filenames_obs.append(obs_hist_filename2)
    runtitles_obs.append('Second obs dataset')    

    timescale='Timescale'
    plottitle='My two datasets'
    myregion=[70.0,-5.0,80.0,10.0]
    plotname_root='compare_as_1dhistograms'

    # Plot 1d histograms of model data with obs overplotted
    plot_1d_histograms(filenames,runtitles,filenames_obs,runtitles_obs,timescale,
                                         myregion,plottitle,plotname_root)
    
    # Set up list of filenames and runtitles for 1d histogram DIFFERENCE plots
    # Please edit plotname_root, plottitle, runtitles and timescale, as well as region, 
    # as required.
    
    filenames=[]
    runtitles=[]
    
    filenames.append([hist_filename1,hist_filename2])
    runtitles.append(['My first dataset','My second dataset'])
    
    timescale='Timescale'
    plottitle='Differences between my two datasets'
    myregion=[70.0,-5.0,80.0,10.0]
    plotname_root='compare_as_1dhist_differences'

    # Plot differences between 1d histograms from 1 model datasets
    plot_1d_histogram_diffs(filenames,runtitles,timescale,
                                         myregion,plottitle,plotname_root)

    # Print message confirming all completed OK
    
    print 'Processing completed OK!'
    
    return
    

def making_histogram_files(data_filename1,data_filename2,hist_filename1,hist_filename2):
    """ 
    Read in data and create histogram cubes, save these to netcdf files.
    
    Arguments:
    * data_filename1,data_filename2
        input filenames of time-varying precipitation data
      
    * hist_filename1,hist_filename2
        chosen filenames for output histogram cubes (netcdf files)
    """ 

    ppndata1=make_hist_maps.read_data_cube(data_filename1)
    ppn_hist_cube1=make_hist_maps.make_hist_ppn(ppndata1)
    iris.save(ppn_hist_cube1, hist_filename1)

    ppndata2=make_hist_maps.read_data_cube(data_filename2)
    ppn_hist_cube2=make_hist_maps.make_hist_ppn(ppndata2)
    iris.save(ppn_hist_cube2, hist_filename2)

    return

def plot_histogram_maps(hist_filename1,hist_filename2,runtitle,plotname_root):
    """ 
    Plot histogram maps
    """
    
    ppn_hist_cube1=make_hist_maps.read_data_cube(hist_filename1)
    ppn_hist_cube2=make_hist_maps.read_data_cube(hist_filename2)
    
    avg_rain_bins_a,avg_rain_bins_frac_a=make_hist_maps.calc_rain_contr(ppn_hist_cube1)
    avg_rain_bins_b,avg_rain_bins_frac_b=make_hist_maps.calc_rain_contr(ppn_hist_cube2)
    
    # (optional) Define how you want to lump the bins together (below is the default)
    
    all_ppn_bounds = [(0.005, 10.), (10., 50.), (50., 100.), (100., 3000.)]

    # Plot as actual contributions for specific region, e.g. 60 to 160E,10S to 10N
    
    plotname='{}_actual_contributions.png'.format(plotname_root)
    plot_hist_maps.plot_rain_contr(avg_rain_bins_a,avg_rain_bins_b,plotname,
                  runtitle,'Timescale 1','Timescale 2',all_ppn_bounds,region=[60.0,-10.0,160.0,10.0])    
    
    # Plot as fractional contributions
    
    plotname='{}_fractional_contributions.png'.format(plotname_root)
    plot_hist_maps.plot_rain_contr(avg_rain_bins_frac_a,avg_rain_bins_frac_b,plotname,
                  runtitle,'Timescale 1','Timescale 2',all_ppn_bounds,region=[60.0,-10.0,160.0,10.0],
		  frac=1)    

    return
  

def plot_1d_histograms(filenames,runtitles,filenames_obs,runtitles_obs,timescale,
                                         myregion,plottitle,plotname_root):
    """
    Plot 1d histograms for a small region.
   
    This example uses histogram cubes pre-calculated from two different model datasets 
    on the same timescale, and compares with those from two observational datasets.
   
    NOTE that the region and the timescale will appear automatically in the plot title
    """
  
    plotname='{}_actual.png'.format(plotname_root)
    plot_hist1d.plot_1dhist(plotname,myregion,filenames,runtitles,plottitle,timescale,
                         filenames_obs,runtitles_obs,log=1)

    plotname='{}_fractional.png'.format(plotname_root)
    plot_hist1d.plot_1dhist(plotname,myregion,filenames,runtitles,plottitle,timescale,
                         filenames_obs,runtitles_obs,frac=1,log=1)

    return


def plot_1d_histogram_diffs(filenames,runtitles,timescale,
                                         myregion,plottitle,plotname_root):
    """
    Plot 1d histograms for a small region.
   
    This example uses histogram cubes pre-calculated from two different model datasets 
    on the same timescale, and compares with those from two observational datasets.
   
    NOTE that the region and the timescale will appear automatically in the plot title
    """
  
    plotname='{}_actual.png'.format(plotname_root)
    plot_hist1d.plot_1dhist(plotname,myregion,filenames,runtitles,plottitle,timescale,
                         log=1)

    plotname='{}_fractional.png'.format(plotname_root)
    plot_hist1d.plot_1dhist(plotname,myregion,filenames,runtitles,plottitle,timescale,
                         frac=1,log=1)

    return
  

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process two model ' 
        'precipitation datasets and compare them with each other and with two '
        'observational datasets.')
    parser.add_argument('data_one', help='first time-varying '
        'precipitation data file')
    parser.add_argument('data_two', help='second time-varying '
        'precipitation data file')
    parser.add_argument('hist_one', help='chosen filename for first '
        'output histogram cube')
    parser.add_argument('hist_two', help='chosen filename for second '
        'output histogram cube')
    parser.add_argument('obs_one', help='filename for first '
        'previously-calculated observed histogram cube')
    parser.add_argument('obs_two', help='filename for second '
        'previously-calculated observed histogram cube')
    args = parser.parse_args()

    main(args.data_one,args.data_two,args.hist_one,args.hist_two,args.obs_one,args.obs_two)


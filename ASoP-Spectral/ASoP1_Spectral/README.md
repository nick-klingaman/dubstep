**(C) British Crown Copyright 2016, Met Office**

=======
Modules
=======

* make_hist_maps.py:
    * extend_cube_with_dimcoord: create a new cube with an extra dimension as specified
    * hist_count: Function to calculate where values fall between bin_min and bin_max
    * read_data_cube: Function to read in data from filename
    * extract_region: Function to extract chosen region
    * make_hist_ppn: Function to create set of histograms for precipitation data input 
    * calc_rain_contr: calculate total and fractional precipitation contributions from histogram counts

* plot_hist_maps.py:
    * _make_panel: make a single panel for the plot
    * plot_rain_contr: Function to plot set of histogram maps as rainfall contributions on two different timescales side by side.

* plot_hist1d.py:
    * guessbounds: Function to guess bounds for gridboxes
    * calc_1d_contr: Function to output 1d histogram of contributions
    * plot_1dhist: Function to plot 1d histograms over given region. Can also plot DIFFERENCES between 1d histograms in pairs of files, if input as list of pairs.

* ASoP1_spectral_main.py: 
    * an additional set of modules included to illustrate how the modules above might be used. 

      This can be called from the command line thus:
        
          $ python2.7 ASoP1_spectral/ASoP1_spectral_main.py one.nc two.nc one_op.nc two_op.nc one_hist.nc two_hist.nc

      where:
        * one.nc, two.nc are the input data files
        * one_op.nc, two_op.nc are the chosen filenames for output files
        * one_hist.nc, two_hist.nc are files of pre-calculated histograms from two control (e.g. observation) datasets.

      A help message is provided:
        
          $ python2.7 ASoP1_spectral/ASoP1_spectral_main.py --help                                                
          usage: ASoP1_spectral_main.py [-h]
                                  data_one data_two hist_one hist_two obs_one
                                  obs_two
        
          Process two model precipitation datasets and compare them with each other and
          with two observational datasets.
        
          positional arguments:
            data_one    first time-varying precipitation data file
            data_two    second time-varying precipitation data file
            hist_one    chosen filename for first output histogram cube
            hist_two    chosen filename for second output histogram cube
            obs_one     filename for first previously-calculated observed histogram cube
            obs_two     filename for second previously-calculated observed histogram
                        cube
        
          optional arguments:
            -h, --help  show this help message and exit



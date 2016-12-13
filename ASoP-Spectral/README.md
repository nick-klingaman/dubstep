# ASoP1-Spectral
**Code to create histogram maps of rainfall intensities on any given timescale.** 

**Part of the ASoP (Analysing Scales of Precipitation) v1.0 package.**

This set of Python modules is sufficient to generate and plot histograms of precipitation amounts on any input timescale and horizontal grid.

The code uses Python2.7 and requires installation of the Iris (see http://scitools.org.uk/iris/) and Matplotlib (see  http://matplotlib.org/) code libraries.

The code uses either kg m-2 day-1 or mm day-1 as the unit for precipitation throughout. Input datasets on other timescales are converted to kg m-2 day-1 or mm day-1 using Iris method ``convert_units``, provided their input units are specified and are in the same form. 

Datasets of accumulated rainfall over specified time periods must include the accumulation time as part of the unit otherwise convert_units will fail. For example, 3-hourly accumulations must be labelled in mm 3h-1 or kg m-2 3h-1 and not just mm or kg m-2.

NOTE that plots are labelled with "mm/day" although it is recognised that this assumes water density = 1 kg m-3.

The method uses 100 bins whose limits are defined thus:

        bin2=np.exp(np.log(0.005)+np.sqrt(np.linspace(0,98,99)*((np.square(np.log(120.)-np.log(0.005)))/59.)))
        bin=np.zeros(100)
        bin[1:]=bin2[0:99]

This samples rainfall intensities in the range 0.005 to 2209 mm/day. Counts of values into element ``[i]`` of the output histogram are made such that ``pbin[i] <= values < pbin[i+1]``. Any values larger(smaller) than that are all lumped into the 100th(first) bin.

Note that this implies that the histogram of events will sum to total no. of intervals in sample.

Histogram counts can be converted to actual and fractional rainfall contributions using:

* Average rain in box over sample = sum over bins(no. events x avg. ppn rate in bin) / (no of intervals in sample)
* Contribution of average rainfall by events in bin (in mm/day) = (no. of events x avg. ppn rate in bin in mm/day)/(no of intervals in sample)
* Fractional contribution = Actual contribution / Average rain in box over sample

By calculating these contributions at many grid-boxes in a region, we can produce maps of the contributions of various precipitation intensity bins to the total precipitation at each grid-box. 

Regional averages of the spectra can also be produced for direct comparison between datasets. This is done by averaging together the histograms calculated at each point in the region, NOT by creating a new histogram using all the points within the region (doing the latter would remove some of the effects of grid-scale temporal variability (as it may be raining at one point when it is not at another), which is what this analysis method was designed to investigate).

Regional averages are included for easier comparison between timescales and datasets, but since they do introduce spatial averaging, they are best done for relatively small (or spatially-consistent) regions only.

======
Author
======

Gill Martin, Met Office, FitzRoy Road, Exeter UK

**(c) British Crown Copyright 2016, Met Office**

**Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.**

You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

========
Citation
========

Users who apply the software resulting in presentations or papers are kindly asked to cite:

* *Klingaman, N. P., Martin, G. M., and Moise, A.: ASoP (v1.0): A set of methods for analyzing scales of precipitation in general circulation models, Geosci. Model Dev., http://www.geosci-model-dev-discuss.net/gmd-2016-161/*

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


=====
Tests
=====

The set of tests is designed to check that the code produces identical results ot the sample files included in the package. This can be used to check whether any changes made have altered the functionality. The test package: 

* Uses test input dataset (a random distribution of "rainfall" data). 
* Generates histogram counts and saves as netCDF file. Should be identical to sample histogram dataset supplied.
* Checks that the fractional contributions sum to 1.0 at all points

Files supplied (in reference subdirectory):

* test_precipitation_dataset.nc : input artifical rainfall dataset (mm/day)
* test_precipitation_hist.nc : output file of histogram counts
* test_precipitation_histmap_act.png : map of actual rainfall contributions
* test_precipitation_histmap_frac.png : map of fractional rainfall contributions
* test_hist1d_act.png : 1D histogram of actual rainfall contributions for sub-region (215E - 225E, 5S - 5N)
* test_hist1d_frac.png : 1D histogram of fractional rainfall contributions for sub-region (215E - 225E, 5S - 5N)

Files generated by tests:

* mytest_histmap.nc
* mytest_histmap_act.png
* mytest_histmap_frac.png
* mytest_hist1d_act.png
* mytest_hist1d_frac.png



# ASoP-Coherence
**Code to diagnose spatial and temporal coherence of precipitation intensities on any given timescale.**

**Part of the ASoP (Aanalysing Scales of Precipitation) v1.0 package.**

This set of Python modules generates a variety of diagnostics and metrics of spatial and temporal coherence (or intermittency) of precipitation amounts on any input timescale and horizontal grid.

## Methods

The code creates the following diagnostics:
* One-dimensional histograms of precipitation, using user-specified bins.
* Two-dimensional histograms of precipitation on successive temporal intervals at the same gridpoint, using the same bins as the one-dimensional histograms.
* Correlations in space and time as functions of the native (input) dataset horizontal grid and temporal sampling.
* Correlations in space and time as functions of physical distance (in km) and physical time (in minutes), which are more useful for comparing correlations between datasets or between the same dataset at different resolutions.
* Summary metrics of spatial and temporal coherence.

These diagnostics are described in greater detail in the peer-reviewed journal article cited below (Klingaman et al., 2017).

## Dependencies
The code uses Python2.7 and requires installation of the following Python packages:
* Iris: http://scitools.org.uk/iris
* Matplotlib: http://matplotlib.org
* Numpy: http://www.numpy.org
* cf-plot: http://ajheaps.github.io/cf-plot

## Author

Nicholas Klingaman, National Centre for Atmospheric Science and Department of Meteorology, University of Reading, Reading, UK

**(c) The Author, 2017**

**Licensed under the Apache License, Version 2.0 (the "License"); you may not use these files except in compliance with the License.**

You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and limitations under the License.

## Citation
Users who apply the software resulting in presentations or papers are kindly asked to cite:

* *Klingaman, N. P., Martin, G. M., and Moise, A.: ASoP (v1.0): A set of methods for analyzing scales of precipitation in general circulation models, Geosci. Model Dev.,10, 57-83, doi:10.5194/gmd-10-57-2017, 2017*

## Modules

There are two files: asop_coherence.py contains all required functions to produce the diagnostics; asop_coherence_example.py shows an example application of these functions to TRMM and CMORPH rainfall data.

* asop_coherence.py: ASoP-Coherence functions
  * read_precip: Function to read in data using iris
  * compute_histogram: Computes 1D and 2D histograms from precipitation data
  * plot_histogram: Plots 2D and 1D histograms
  * compute_equalgrid_corr: Computes spatial correlations, including lagged correlations, as a function of native gridpoints, for a single dataset
  * plot_equalgrid_corr: Plots spatial correlations, including lagged correlations, as a function of native gridpoints, for a single dataset
  * compute_equalarea_corr: Computes spatial correlations as a function of physical distance (only coincident correlations, not lagged), for a single dataset
  * plot_equalarea_corr: Plots spatial correlations as a function of physical distance, for either single or multiple datasets
  * compute_autocorr: Computes lagged autocorrelations as a function of physical time, for a single dataset
  * plot_autocorr: Plots lagged autocorrelations as a function of physical time, for either single or multiple datasets
 
* asop_coherence_example.py: ASoP-Coherence example
  * get_dictionary: Demonstrates how to build a dictionary of dataset information needed by the ASoP-Coherence functions
  * __main__: Main procedure that calls each of the ASoP-Coherence functions in turn.

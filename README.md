# asop_duration
Author: Amulya Chevuturi (a.chevuturi@reading.ac.uk)

Adapted from Analysing Scales of Precipitaiton (ASoP) diagnostics https://github.com/nick-klingaman/ASoP
Reference: Klingaman, N. P., Martin, G. M., and Moise, A.: ASoP (v1.0): a set of methods for analyzing scales of precipitation in general circulation models, Geosci. Model Dev., 10, 57–83, https://doi.org/10.5194/gmd-10-57-2017, 2017.

2D and 1D histograms for Precipitation Intensity vs. Time Duration:
To analyse the differences in rainfall variability, we use a new set of diagnostics (“duration diagnostics”) which 
compare the precipitation duration characteristics across model resolutions and observation. Using these diagnostics, 
we measure the probability of the duration of discrete precipitation bins and probability of length of dry spells. 
These diagnostics can be applied to data across resolutions and at any time step. In the first metric, we construct 
a two-dimensional (2-D) histogram of discrete precipitation bins against discrete duration bins, where duration is 
the length of time for which the precipitation intensity (within a particular discrete bin) occurs continuously. We 
calculate the 2-D histogram across all grid points within the analysis domain, and is normalized by the time length and 
grid size of the dataset, to be able to compare across datasets. The second metric produces a 1-D histogram for the 
frequency of length/duration of dry spells, where dry spell is defined as all the time steps having precipitation below 
0.1 mm/day. This metric calculates the frequency of occurrence of each duration bin of dry spell for all grid points within 
the analysis domain, which is normalized by the size of the dataset, to be able to compare across datasets. 

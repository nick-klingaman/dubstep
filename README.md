# asop_duration
2D and 1D histograms for Precipitation vs. Duration

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

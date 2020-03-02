from __future__ import division
import numpy as np
from netCDF4 import Dataset
from scipy.stats import linregress as ols
#from scipy.stats.stats import pearsonr

###########################################################
#Metrics descripbed below
###########################################################



###########################################################
'''
 Metric:
   Zeng's Gamma
 Description:
   This function (cal_gam) calculates Zeng's Gamma for input array(s).

   Zeng's gam is the correlation between evaporation and precipitation
   scaled by the standard deviation of the evaporation and normalized by the
   standard deviation of precipitation to keep the index dimensionless.

 Requirements:
   Takes precipitation (mm/day) and evapotranspiration (mm/day) as input variables. Data
   should be formatted as numpy array at each grid point, constrained to the same time period.

 Output:
   One single value for Gamma. 

 Reference:
   Zeng, X. et al., 2010. Comparison of Land-Precipitation Coupling Strength
   Using Observations and Models. Journal of Hydrometeorology, 11, 979-994.
'''

def cal_gam(pr,et):
  slope, intercept, r, p, std_err = ols(pr,et)
  if p < 0.01:
    gam = (r)*(np.nanstd(et)/np.nanstd(pr))
  else:
    gam = np.nan
  return(gam)
  
'''
OLDER FUNCTION
  r, p = pearsonr(pr,et)
  gam = (r)*(np.nanstd(e)/np.nanstd(p))
  return(gam)
'''
###########################################################



#############################################################################
'''
 Metric:
   Terrestrial Coupling Index (TCI)

 Description:
    This script calculates the TCI for all grid cells in input array(s).

    The terrestrial coupling index determines the coupling of soil moisture
    and surface fluxes. It is calculated as the slope of the sol moisture-
    surface flux relationship, weighted by the standard deviation in soil
    moisture to determine the degree to which soil moisture changes drive
    surface flux variability.
    
    Can also calculate Temperature Evapotranpiration Metric (TET) by using
    met.cal_tci(t2ano, etano, ols_out='r',weighting=False)

 Requirements:
    Takes soil moisture and evapotranspiration as input variables, but could
    be applied to other variables. Data should be formatted as numpy arrays,
    constrained to the same time period.

 Output:
    One single value for the TCI.
    
 References:
    Dirmeyer, P. A. 2011. The terrestrial segment of soil moisture and climate
    coupling. Geophysical Research Letters, 38, L16702.

    Guo, Z. et al., 2006. Glace: The Global Land and Atmosphere Coupling
    Experiment. Part II: Analysis. Journal of Hydrometeorology, 7, 611-625.
'''

def cal_tci(s, e, ols_out='slope',weighting=True):
  slope, intercept, r, p, std_err = ols(s,e)
  # If significant then save value otherwise NaN
  if p < 0.01:
    if ols_out == 'slope':
       tc = slope
    elif ols_out == 'r':
       tc = r
  else:
    tc = -999.0
  # Weight by variability of denominator (see Dirmeyer et al., 2011). This emphasises places where actual impact is large
  if weighting is True:
    if (tc != -999.0):
      tci = tc * np.nanstd(s)
    else:
      tci = tc * np.nan
    return(tci)

  if weighting is False:
    if (tc != -999.0):
      tci = tc * 1
    else:
      tci = tc * np.nan
    return(tci)

'''
OLDER FUNCTION
  r, p = pearsonr(s,e)
  if p < 0.05:
    tci = r*(np.nanstd(s))
  else:
    tci = np.nan

def cal_tet(t2,et):
  #slope, intercept, r, p, std_err = ols(t2,et)
  r, p = pearsonr(t2,et)
  if p < 0.01:
    tet = r.copy()
  else:
    tet = np.nan
  return(tet)
'''
###########################################################




###########################################################
'''
 Metric:
 Two-legged Metric
 
 Description:
    This script is the diagnostics and plotting script for the global
    two-legged metric.
    
    The two-legged coupling metric is used to trace energy or moisture 
    feedback pathways from the surface to the atmosphere in a mechanistic 
    way. The metric is based on having a physical understanding of the factors
    that control interactions between the land and the atmosphere, and can 
    be used to identify areas where land-atmosphere coupling is particularly strong.
    The feedback pathway is broken down into two stages: the surface leg,
    which measures the strength of regression between a surface state variable (S)
    and a surface flux variable (F), and the atmospheric leg, which measures the
    regression relationship between the flux and an atmospheric variable (A).

 Requirements:
    Takes three variables relating to the land surface, land-atmosphere flux
    and the atmosphere. Data should be constrained to the same time period.
    This function uses the above defined TCI function (cal_tci)
    and it needs to be defined before it. 
    
 Output:
    Three values -
    first is the surface leg (TCI between first and second input arrays),
    second is the atmosphere leg (TCI between second and third input arrays),
    third is the total feed back (mulitple of first and second output values).
  
 References:
    Dirmeyer, P. A. 2011. The terrestrial segment of soil moisture and climate
    coupling. Geophysical Research Letters, 38, L16702.

    Dirmeyer, P. A. et al., 2014. Intensified land surface control on boundary 
    layer growth in a changing climate. Geophysical Research Letters, 41, 1290-1294.

    Guo, Z. et al., 2006. Glace: The Global Land and Atmosphere Coupling
    Experiment. Part II: Analysis. Journal of Hydrometeorology, 7, 611-625.
'''

def cal_leg(s, e, p, ols_out='slope',weighting=True):
  sur = cal_tci(s, e, ols_out=ols_out,weighting=weighting)
  atm = cal_tci(e, p, ols_out=ols_out,weighting=weighting)
  leg = sur*atm
  return(sur,atm,leg)

###########################################################

def rmm_index(SY,EY,MON):

  #Reads in MJO Index vales as the phase (ph) and amplitude (am) for SY (Start Year) till EY (End Year) for the MON (Numpy array of Months)
  #from the a txt file called rmm.txt downloaded from http://www.bom.gov.au/climate/mjo/graphics/rmm.74toRealtime.txt
  #MON can br define the months of a season uptp any months. Should be written in ascending order or the months and should be in the array form
  #For example for NDJFM season, your MON value should be written as: MON = np.array([1,2,3,11,12])
  
  #Following is an example of how to use this definition in your code for deliating the months of NDJFM in the years between 1999-2010:
  
  #import read_mjo_index as rmm
  #SY = 1999
  #EY = 2010
  #MON = np.array([1,2,3,11,12])
  #am,ph = rmm.rmm_index(SY,EY,MON)

  import numpy as np
  import pandas as pd
  
  df = pd.read_table('/gws/nopw/j04/ncas_climate_vol1/users/amulya/data/obs/mjo/rmm.txt', delim_whitespace=True, skiprows=2, usecols=range(7), header=None) 

  amplt = df[df.columns[6]].values
  phase = df[df.columns[5]].values
  yrs   = df[df.columns[0]].values
  mons  = df[df.columns[1]].values

  amplt = amplt[(yrs>=SY)&(yrs<=EY)]
  phase = phase[(yrs>=SY)&(yrs<=EY)]
  mons  = mons[(yrs>=SY)&(yrs<=EY)]
  yrs   = yrs[(yrs>=SY)&(yrs<=EY)]
  
  a=[]; p =[]
  for i in range(SY,EY+1):
    a.append([0]*len(MON))
    p.append([0]*len(MON))

  for y in range(SY,EY+1):
    for m in range(len(MON)):
      a[y-SY][m] = amplt[(yrs==y)&(mons==MON[m])]
      p[y-SY][m] = phase[(yrs==y)&(mons==MON[m])]

  return a,p



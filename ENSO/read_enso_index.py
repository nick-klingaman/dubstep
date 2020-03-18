def on_index(SY,EY,MON):

  #Reads in ONI vales for SY (Start Year) till EY (End Year) for the MON (Numpy array of Months)
  #from the a txt file called oni.txt downloaded from https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/detrend.nino34.ascii.txt
  #MON can define the months of a season uptp only 6 months. Should be written in ascending order or the months and should be in the array form
  #For example for NDJFM season, your MON value should be written as: MON = np.array([1,2,3,11,12])
  
  #Following is an example of how to use this definition in your code for deliating the months of NDJFM in the years between 1999-2010:
  
  #import read_enso_index as nino
  #SY = 1999
  #EY = 2010
  #MON = np.array([1,2,3,11,12])
  #oni = nino.on_index(SY,EY,MON)

  #To just do the seasons instead of the months, you would have to remove the first (1999) J-F-M and the last (2010) N-D months
  #But this cannot be done by this function, you will have to add the index yourself in your script which is something like this:
  #oni = oni[3:-2]
     
  import numpy as np
  import pandas as pd
  
  df = pd.read_table('oni.txt', delim_whitespace=True, header=0) 
  anom = df['ANOM'].values
  yrs  = df['YR'].values
  mons  = df['MON'].values
  anom = anom[(yrs>=SY)&(yrs<=EY)]
  mons  = mons[(yrs>=SY)&(yrs<=EY)]

  if MON.size==1:
    index = anom[(mons==MON[0])]
  if MON.size==2:
    index = anom[(mons==MON[0])|(mons==MON[1])]
  if MON.size==3:
    index = anom[(mons==MON[0])|(mons==MON[1])|(mons==MON[2])]
  if MON.size==4:
    index = anom[(mons==MON[0])|(mons==MON[1])|(mons==MON[2])|(mons==MON[3])]
  if MON.size==5:
    index = anom[(mons==MON[0])|(mons==MON[1])|(mons==MON[2])|(mons==MON[3])|(mons==MON[4])]
  if MON.size==6:
    index = anom[(mons==MON[0])|(mons==MON[1])|(mons==MON[2])|(mons==MON[3])|(mons==MON[4])|(mons==MON[5])]
  if MON.size==7:
    index = anom[(mons==MON[0])|(mons==MON[1])|(mons==MON[2])|(mons==MON[3])|(mons==MON[4])|(mons==MON[5])|(mons==MON[6])]
  if MON.size==8:
    index = anom[(mons==MON[0])|(mons==MON[1])|(mons==MON[2])|(mons==MON[3])|(mons==MON[4])|(mons==MON[5])|(mons==MON[6])|(mons==MON[7])]
  if MON.size==9:
    index = anom[(mons==MON[0])|(mons==MON[1])|(mons==MON[2])|(mons==MON[3])|(mons==MON[4])|(mons==MON[5])|(mons==MON[6])|(mons==MON[7])|(mons==MON[8])]
  if MON.size==10:
    index = anom[(mons==MON[0])|(mons==MON[1])|(mons==MON[2])|(mons==MON[3])|(mons==MON[4])|(mons==MON[5])|(mons==MON[6])|(mons==MON[7])|(mons==MON[8])|(mons==MON[9])]
  if MON.size==11:
    index = anom[(mons==MON[0])|(mons==MON[1])|(mons==MON[2])|(mons==MON[3])|(mons==MON[4])|(mons==MON[5])|(mons==MON[6])|(mons==MON[7])|(mons==MON[8])|(mons==MON[9])|(mons==MON[10])]
  if MON.size==12:
    index = anom[(mons==MON[0])|(mons==MON[1])|(mons==MON[2])|(mons==MON[3])|(mons==MON[4])|(mons==MON[5])|(mons==MON[6])|(mons==MON[7])|(mons==MON[8])|(mons==MON[9])|(mons==MON[10])|(mons==MON[11])]

  return index



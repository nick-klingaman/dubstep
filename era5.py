from __future__ import division
#import cf, cfplot as cfp
import numpy as np
from netCDF4 import Dataset
from scipy.stats.stats import pearsonr
from scipy.stats import linregress as ols
import metrics as met

LAT=47
LON=47
TIME=55
WEEKS=5
weeki = np.array([0,8,16,24])
SY = 1999
EY = 2010
MON = np.array([1,2,3,11,12])

yr = [str(x) for x in np.arange(SY,EY+1)]
mon = [str(x).zfill(2) for x in MON]

yrmon=[]
for i in range(len(yr)):
    for j in range(len(mon)):
        yrmon.append("".join((yr[i],mon[j])))

yrmon = yrmon[3:-2]

pr = np.empty((TIME,len(weeki),WEEKS,LON,LAT)).squeeze()
et = np.empty((TIME,len(weeki),WEEKS,LON,LAT)).squeeze()
sm = np.empty((TIME,len(weeki),WEEKS,LON,LAT)).squeeze()
t2 = np.empty((TIME,len(weeki),WEEKS,LON,LAT)).squeeze()
for i in range(len(yrmon)):
    #pr[i,:,:,:,:]=Dataset('/gws/nopw/j04/ncas_climate_vol1/users/myoung02/datasets/DUBSTEP_data/SAmerica_CHIRPS_weekly_'+yrmon[i]+'.nc', 'r').variables['week_precip'][:]
    pr[i,:,:,:,:]=Dataset('/gws/nopw/j04/klingaman/amulya/data/surface_s2s_skill/precip/dubstep_output/ERA5/SA_ERA5_precip_weekly_'+yrmon[i]+'.nc', 'r').variables['week_precip'][:]
    et[i,:,:,:,:]=Dataset('/gws/nopw/j04/klingaman/amulya/data/surface_s2s_skill/slhf/dubstep_output/ERA5/SA_ERA5_ET_weekly_'+yrmon[i]+'.nc', 'r').variables['week_ET'][:]
    sm[i,:,:,:,:]=Dataset('/gws/nopw/j04/klingaman/amulya/data/surface_s2s_skill/soilm/dubstep_output/ERA5/SA_ERA5_sm20_weekly_'+yrmon[i]+'.nc', 'r').variables['week_sm20'][:]
    t2[i,:,:,:,:]=Dataset('/gws/nopw/j04/klingaman/amulya/data/surface_s2s_skill/t2m/dubstep_output/ERA5/SA_ERA5_t2m_weekly_'+yrmon[i]+'.nc', 'r').variables['week_t2m'][:]

pr[pr<0] = 0.0
sm[sm<0] = 0.0

#pr = pr.reshape((TIME,len(weeki),WEEKS,LON,LAT))
prano = pr - np.nanmean(pr, axis=(0,1), keepdims=True)
prano = prano.reshape((TIME*len(weeki),WEEKS,LON,LAT))
#et = et.reshape((TIME,len(weeki),WEEKS,LON,LAT))
etano = et - np.nanmean(et, axis=(0,1), keepdims=True)
etano = etano.reshape((TIME*len(weeki),WEEKS,LON,LAT))
#sm = sm.reshape((TIME,len(weeki),WEEKS,LON,LAT))
smano = sm - np.nanmean(sm, axis=(0,1), keepdims=True)
smano = smano.reshape((TIME*len(weeki),WEEKS,LON,LAT))
#t2 = t2.reshape((TIME,len(weeki),WEEKS,LON,LAT))
t2ano = t2 - np.nanmean(t2, axis=(0,1), keepdims=True)
t2ano = t2ano.reshape((TIME*len(weeki),WEEKS,LON,LAT))

gam = np.empty((WEEKS,LON,LAT)).squeeze()
tci = np.empty((WEEKS,LON,LAT)).squeeze()
tet = np.empty((WEEKS,LON,LAT)).squeeze()
pleg = np.empty((3,WEEKS,LON,LAT)).squeeze()
tleg = np.empty((3,WEEKS,LON,LAT)).squeeze()
for w in range(WEEKS):
  for x in range(LON):
    for y in range(LAT):
      gam[w,x,y] = met.cal_gam(prano[:,w,x,y], etano[:,w,x,y])
      tci[w,x,y] = met.cal_tci(smano[:,w,x,y], etano[:,w,x,y], ols_out='slope',weighting=True)
      tet[w,x,y] = met.cal_tci(t2ano[:,w,x,y], etano[:,w,x,y], ols_out='r',weighting=False)
      pleg[:,w,x,y] = met.cal_leg(smano[:,w,x,y], etano[:,w,x,y], prano[:,w,x,y], ols_out='slope',weighting=True)
      tleg[:,w,x,y] = met.cal_leg(smano[:,w,x,y], etano[:,w,x,y], t2ano[:,w,x,y], ols_out='slope',weighting=True)
  print w
    
#gam = np.nanmean(gam, axis=1)
#tci = np.nanmean(tci, axis=1)
#leg = np.nanmean(leg, axis=1)
np.savez('../output/era5', gam=gam, tci=tci, tet=tet, pleg=pleg, tleg=tleg)






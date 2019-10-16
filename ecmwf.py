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

ENS=11
LAG=3
ALL=ENS*LAG

yr = [str(x) for x in np.arange(SY,EY+1)]
mon = [str(x).zfill(2) for x in MON]

yrmon=[]
for i in range(len(yr)):
    for j in range(len(mon)):
        yrmon.append("".join((yr[i],mon[j])))

yrmon = yrmon[3:-2]

pr = np.empty((TIME,len(weeki),WEEKS,LAG,ENS,LON,LAT)).squeeze()
et = np.empty((TIME,len(weeki),WEEKS,LAG,ENS,LON,LAT)).squeeze()
sm = np.empty((TIME,len(weeki),WEEKS,LAG,ENS,LON,LAT)).squeeze()
t2 = np.empty((TIME,len(weeki),WEEKS,LAG,ENS,LON,LAT)).squeeze()
for i in range(len(yrmon)):
    pr[i,:,:,:,:,:]=Dataset('/gws/nopw/j04/ncas_climate_vol1/users/myoung02/datasets/DUBSTEP_data/SAmerica_ECMWF_ensemble_weekly_hindcasts_lead1-5_'+yrmon[i]+'.nc', 'r').variables['week_precip'][:]
    et[i,:,:,:,:,:]=Dataset('/gws/nopw/j04/klingaman/amulya/data/surface_s2s_skill/slhf/dubstep_output/ECMWF/SA_ECMWF_ensemble_ET_weekly_hindcasts_lead1-5_'+yrmon[i]+'.nc', 'r').variables['week_ET'][:]
    sm[i,:,:,:,:,:]=Dataset('/gws/nopw/j04/klingaman/amulya/data/surface_s2s_skill/soilm/dubstep_output/ECMWF/SA_ECMWF_ensemble_sm20_weekly_hindcasts_lead1-5_'+yrmon[i]+'.nc', 'r').variables['week_sm20'][:]
    t2[i,:,:,:,:,:]=Dataset('/gws/nopw/j04/klingaman/amulya/data/surface_s2s_skill/t2m/dubstep_output/ECMWF/SA_ECMWF_ensemble_t2m_weekly_hindcasts_lead1-5_'+yrmon[i]+'.nc', 'r').variables['week_t2m'][:]

pr[pr<0] = 0.0
sm[sm<0] = 0.0

pr = pr.reshape((TIME,len(weeki),WEEKS,ALL,LON,LAT))
prano = pr - np.nanmean(pr, axis=(0,1), keepdims=True)
prano = prano.reshape((TIME*len(weeki),WEEKS,ALL,LON,LAT))
et = et.reshape((TIME,len(weeki),WEEKS,ALL,LON,LAT))
etano = et - np.nanmean(et, axis=(0,1), keepdims=True)
etano = etano.reshape((TIME*len(weeki),WEEKS,ALL,LON,LAT))
sm = sm.reshape((TIME,len(weeki),WEEKS,ALL,LON,LAT))
smano = sm - np.nanmean(sm, axis=(0,1), keepdims=True)
smano = smano.reshape((TIME*len(weeki),WEEKS,ALL,LON,LAT))
t2 = t2.reshape((TIME,len(weeki),WEEKS,ALL,LON,LAT))
t2ano = t2 - np.nanmean(t2, axis=(0,1), keepdims=True)
t2ano = t2ano.reshape((TIME*len(weeki),WEEKS,ALL,LON,LAT))

gam = np.empty((WEEKS,ALL,LON,LAT)).squeeze()
tci = np.empty((WEEKS,ALL,LON,LAT)).squeeze()
tet = np.empty((WEEKS,ALL,LON,LAT)).squeeze()
pleg = np.empty((3,WEEKS,ALL,LON,LAT)).squeeze()
tleg = np.empty((3,WEEKS,ALL,LON,LAT)).squeeze()
for w in range(WEEKS):
  for e in range(ALL):
    for x in range(LON):
      for y in range(LAT):
        gam[w,e,x,y] = met.cal_gam(prano[:,w,e,x,y], etano[:,w,e,x,y])
	tci[w,e,x,y] = met.cal_tci(smano[:,w,e,x,y], etano[:,w,e,x,y], ols_out='slope',weighting=True)
	tet[w,e,x,y] = met.cal_tci(t2ano[:,w,e,x,y], etano[:,w,e,x,y], ols_out='r',weighting=False)
	pleg[:,w,e,x,y] = met.cal_leg(smano[:,w,e,x,y], etano[:,w,e,x,y], prano[:,w,e,x,y], ols_out='slope',weighting=True)
	tleg[:,w,e,x,y] = met.cal_leg(smano[:,w,e,x,y], etano[:,w,e,x,y], t2ano[:,w,e,x,y], ols_out='slope',weighting=True)
  print w
    
gam = np.nanmean(gam, axis=1)
tci = np.nanmean(tci, axis=1)
tet = np.nanmean(tet, axis=1)
pleg = np.nanmean(pleg, axis=2)
tleg = np.nanmean(tleg, axis=2)

np.savez('../output/ecmwf', gam=gam, tci=tci, tet=tet, pleg=pleg, tleg=tleg)





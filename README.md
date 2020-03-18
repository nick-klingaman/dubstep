# dubstep
Diagnostics from the CSSP Brazil DUBSTEP project on sub-seasonal rainfall

This repository incorporates several pieces of code, which are maintained at their own repositories.  This repository is updated irregularly from those repositories, so we recommend checking individual repositories for the latest versions.

## "s2s_skill" - Sub-seasonal lagged ensembles for weekly rainfall
Code to compute weekly rainfall from S2S forecasts, which can be adapted to any other single-level field of interest, and to compute lagged ensembles from ECMWF and NCEP forecasts (relative to the less-frequent UKMO dates) is contained in "s2s_skill".  The code can also compute climatologies and anomalies from these, as well as forecast metrics such as the correlation coefficient (against observations) and the Brier Skill Score.

## "conditional_skill" - Conditional skill analysis
Code to compute prediction performance for sub-seasonal forecasts, conditional upon the phase of large-scale variability, such as the Madden-Julian Oscillation or the El Nino-Southern Oscillation.  The code is maintained here: https://github.com/achevuturi/conditional_skill_analysis.  Original code by Amulya Chevuturi (University of Reading).

## "WAM-2layers" - Water Accounting Model
The Water Accounting Model, which can trace evaporation to precipitation or precipitation to evaporation, for a specified region, in reanalysis or model data.  This code is maintained here: https://github.com/ruudvdent/WAM2layersPython.  Original code by Ruud van der Ent (Technical University of Delft), modified by Liang Guo (University of Reading).

## "asop" - Analysis of Scales of Precipitation
Code to compute intensity spectra of precipitation, as well as spatial and temporal coherence of rainfall features, in observations and model output.  This code is maintained here: https://github.com/nick-klingaman/ASoP.  Original code by Gill Martin (Met Office), Nick Klingaman (University of Reading) and Aurel Moise (Bureau of Meteorology).

## "asop_duration" - Wet-spell and dry-spell duration
Adaptation of ASoP diagnostics to compute a distribution of wet-spell or dry-spell duration, for given precipitation thresholds.  This code is maintained here: https://github.com/achevuturi/asop_duration.  Original code by Amulya Chevuturi (University of Reading).

## "coupling_strength" - Land-atmosphere coupling strengths
Code to compute sub-seasonal land-atmosphere coupling strengths in reanalysis or S2S models, from weekly averaged precipitation, evaporation, soil moisture and 2m air temperature.  Includes instantaneous and lagged correlations.  This code is maintained here: https://github.com/achevuturi/coupling_strength_scripts.  Original code by Amulya Chevuturi (University of Reading).

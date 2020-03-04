# dubstep
Diagnostics from the CSSP Brazil DUBSTEP project on sub-seasonal rainfall

Description of scripts (M. Young 04/03/2020
):

1) prepare_chirps.py
Subsets and regrids (to 1.5d) CHIRPS into daily files for S2S analysis over Brazil

2) grab_data.py
Individual functions to read in CHIRPS and the various S2S hindcast datasets (UKMO,NCEP,ECMWF and BAM). Called in following scripts 2), 3), 4), 5)

3) save_weekly_forecasts_bam.py
- saves weekly means of BAM hindcasts (6-hrly) into individual netcdf files for each month/year. BAM has 2 start dates per month.
- saves the equivalent weekly mean CHIRPS data for evaluation of BAM at the same dates.

4) save_weekly_forecasts_ukmo_ncep.py
- saves weekly means of UKMO hindcasts (4 starts per month) and NCEP hindcasts. For NCEP, saves hindcasts of a given week from the previous 7 start dates (i.e. 7 lags) relative to each given UKMO start date.
- saves the equivalent weekly mean CHIRPS data for evaluation of UKMO, NCEP and ECMWF at the same dates.

5) save_weekly_forecasts_ecmf.py
- saves weekly means of ECMWF hindcasts relative to UKMO weekly forecast start dates. Saves hindcasts from the three previous ECMWF start dates (i.e. 3 lags) closest to the given UKMO start date.

6) save_weekly_hindcasts_probabilities.py
- saves tercile probabilities for BAM, UKMO, NCEP and ECMWF hindcasts using the both the individual model ensemble members and the lagged ensembles (BAM = 11 members; UKMO = 7 members; NCEP = 4 members * 7 lags; ECMWF = 11 members * 3 lags)
- saves equivalent observed tercile category outcomes (1 or 0) from the CHIRPS observations

6) read_hindcasts_for_analysis.py
- reads in weekly mean hindcast and CHIRPS data saves in scripts 2), 3), 4), 5) and masks the data for the period of interest (e.g. hindcasts of precipiation for forecast weeks falling within a given season, e.g. NDJFM)


7) analyse_skill_extseason.py
- Calls 'read_hindcasts_for_analysis.py' to read the data.
- Computes anomaly correlations, RMSE, bias, Brier Skill Score,
- Computes statistical significance of anomaly correlations by applying a two-tailed t-test using the adjusted sample size due to autocorrelation (lag 1)
- saves the resulting gridded statistics in '.npy' files

8) plot_skill_extseason.py
- Plots gridded skill scores produced from 'analyse_skill_extseason.py' as spatial maps and regional average skill vs lead time.

9) date_str.py
Functions to sort out date strings.

# conditional_skill_analysis
DUBSTEP Activity on Conditional Skill Analysis

####################################################################################

Scripts for analysing skill in S2S dataset based on large scale modes of variability (ENSO and MJO) 

Author: Amulya Chevuturi (a.chevuturi@reading.ac.uk)

Adapted from work of Matt Young and Nick Klingaman (NCAS, University of Reading, UK)

Reference: https://github.com/nick-klingaman/dubstep

Each folder (ENSO or MJO) contains scripts to calculate conditional analysis on the respective large-scale modes of variability

Script read_enso_index.py or read_mjo_index.py reads in the index from the txt file that can be downloaded online; oni.txt or rmm.txt respectivey.

date_str.py provides the date and time structure and grab_data.py reads the data in a particular format for the read_hindcast_analysis.py script. Script read_hindcast_analysis.py reads in all the data for analysis using the analyse_skill_extseason.py. analyse_skill_extseason.py writes out the composites based on phases of ENSO or MJO. 

Bootstrapping for calculating significance is done using enso_boot.py (mjo_boot.py) script. Scripts enso_ano.py, enso_ccdiff.py and enso_ccreg.py (mjo_ano.py, mjo_ccdiff.py and mjo_ccreg.py) are plotting scripts. 

####################################################################################


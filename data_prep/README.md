# README

This folder contains codes used to reproduce the clean input data.

## get_data_capacity.R
Run this script with new data from hospital association to update hospital data stream. Detailed steps:
+ Save spreadsheet from hospital association in the format of `COVID Data Extract (m-dd-yy).csv` in the data folder (only locally will be ignored from git commit).
+ Update the line `data_stream <- "4-29-20"` in the get_data_capacity.R file to the same `m-dd-yy` string in the filename.
+ Source the script. Some built-in checking will be performed to make sure data is in the same format as before.
+ Three files will be updated in the `data/` folder
    + ct_current_hosp.csv
    + ct_cum_hosp.csv
    + ct_hosp_cap.csv 
`
## get_data_population.R
This script downloads ACS estimates and reproduces `data/CT_demog_ACS2018_5year_estimates.csv`.

## plot_map.R
This script tests shapefile and generates adj matrix.

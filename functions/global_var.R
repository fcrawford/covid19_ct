########################################################################
## This script holds all global variables, functions and dependencies
########################################################################
if("is.calibration" %in% ls() == FALSE){
	is.calibration <- FALSE
}

library(deSolve)
library(yaml)
library(lubridate)
library(plyr)
library(reshape2)
library(dplyr)
library(foreign)
library(mvtnorm)
library(readr)

if(!is.calibration){
	# things that not needed on cluster
	library(rgdal)
	library(SUMMER)
	library(ggplot2)
	library(mapproj)
	library(tidyr)
	library(patchwork)
	library(imputeTS)	
   library(forecast)
	# Spatial polygons
	CT_MAP <- readOGR("../map/wgs84/countyct_37800_0000_2010_s100_census_1_shp_wgs84.shp", verbose=FALSE)
}

########################
# Global variables of set dates: do not change
# starting day
day0 = ymd("2020-03-01")

# actual dates: do not change

lockdown_start_date = ymd("2020-03-17") # lockdown effect starts with school closure 
schools_reopen_date = ymd("2020-09-01") # school reopening

INT_START_DATES = as.list(c(lockdown_start_date = lockdown_start_date, schools_reopen_date = schools_reopen_date))


########################
source("../functions/model.R")
source("../functions/run_ct_model.R")
source("../functions/get_data.R")
source("../functions/intervention_functions.R")
source("../functions/truncated_distributions.R")
source("../functions/plot_functions.R")

#################################################### 
# Other global variables that should not be changed
####################################################
global_dat <- get_ct_data(day0=day0)

# State and county observed data
DAT_CT_STATE <- global_dat$dat_ct_state
DAT_CT_COUNTY <- global_dat$dat_ct_county
# Estimated hospitalizations from congregate settings
HOSP_CONG <- global_dat$hosp_cong
# incidence measures: positive testing (community) and CLI ED visits
INCIDENCE <- global_dat$incidence
# County capacity functions, ordered
COUNTY_CAPACITIES <- global_dat$county_capacities
# Mobility 
MOB <- global_dat$mob
# daily PCR testing
TESTING <- global_dat$testing
# relative change in hospital death hazard
DEATH_HAZ <- global_dat$smooth_hdeath_haz
# Spatial adj  matrix
CT_ADJ <- global_dat$adj
# Population count vector in the same order as CT_ADJ
CT_POPULATIONS <- global_dat$populations
CT_POP = sum(CT_POPULATIONS) # total CT population
# County names in the correct order
CT_NAMES <- colnames(CT_ADJ)
# Number of regions
nregions <- length(CT_NAMES)


# distribution of initial numbers exposed by county: used to calculate state0 on day0
E_INIT_COUNTY <- c(0.4,    # 1. "Fairfield" 
                   0.003,  # 2. "New London" 
                   0,      # 3. "Litchfield" 
                   0.001,  # 4. "Windham"  
                   0.001,  # 5. "Tolland" 
                   0.35,   # 6. "Hartford"
                   0.005,  # 7. "Middlesex" 
                   0.24)   # 8. "New Haven"

# initial numbers exposed in each county from v1
# E_init_state0 <- c(7.5,    # 1. "Fairfield" 
#                  0.13,    # 2. "New London" 
#                  0,       # 3. "Litchfield" 
#                  0.04,    # 4. "Windham"  
#                  0.05,    # 5. "Tolland" 
#                  6.5,     # 6. "Hartford"
#                  0.25,    # 7. "Middlesex" 
#                  4.5)     # 8. "New Haven"





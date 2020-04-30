########################################################################
## This script holds all global variables, functions and dependencies
########################################################################

library(deSolve)
library(yaml)
library(lubridate)
library(plyr)
library(reshape2)
library(dplyr)
library(rgdal)
library(SUMMER)
library(ggplot2)
library(mapproj)
library(lubridate)

########################
# Global variables of set dates: do not change
# starting day
day0 = ymd("2020-03-01")
# actual dates: do not change
state_schools_close = dmy("13/03/2020")
state_lockdown_order = dmy("20/03/2020") # order date
state_lockdown_start = dmy("23/03/2020") # actual start date

########################
source("../functions/model.R")
source("../functions/run_ct_model.R")
source("../functions/get_data.R")
source("../functions/intervention_functions.R")
source("../functions/truncated_distributions.R")
source("../functions/plot_functions.R")


########################
# Load data objects
global_dat <- get_ct_data(day0=day0)

DAT_CT_STATE <- global_dat$dat_ct_state
DAT_CT_COUNTY <- global_dat$dat_ct_county
COUNTY_CAPACITIES <- global_dat$county_capacities
CT_MAP <- global_dat$CTmap
CT_ADJ <- global_dat$adj
CT_POPULATIONS <- global_dat$populations
CT_NAMES <- colnames(CT_ADJ)


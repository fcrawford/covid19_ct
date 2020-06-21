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
library(tidyr)
library(patchwork)

########################
# Global variables of set dates: do not change
# starting day
day0 = ymd("2020-03-01")

# actual dates: do not change
state_schools_close = dmy("17/03/2020")
#state_lockdown_order = dmy("20/03/2020") # order date
state_lockdown_start = dmy("23/03/2020") # actual lockdown start date
state_phase1_start = dmy("20/05/2020") # start of phase 1 lockdown release
state_phase2_start = dmy("17/06/2020") # start of phase 2 lockdown release

INT_START_DATES = as.list(c(state_schools_close=state_schools_close, 
                      state_lockdown_start=state_lockdown_start, 
                      state_phase1_start=state_phase1_start, 
                      state_phase2_start=state_phase2_start))



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

# State and county observed data, for plotting purposes only
DAT_CT_STATE <- global_dat$dat_ct_state
DAT_CT_COUNTY <- global_dat$dat_ct_county
# County capacity functions, ordered
COUNTY_CAPACITIES <- global_dat$county_capacities
# Mobility 
MOB <- global_dat$mob
# daily PCR testing
TESTING <- global_dat$testing
# Spatial polygons
CT_MAP <- global_dat$CTmap
# Spatial adj  matrix
CT_ADJ <- global_dat$adj
# Population count vector in the same order as CT_ADJ
CT_POPULATIONS <- global_dat$populations
# County names in the correct order
CT_NAMES <- colnames(CT_ADJ)


# initial numbers exposed in each county: used to generate state0
E_init_state0 <- c(7.5,    # 1. "Fairfield" 
                  0.13,    # 2. "New London" 
                  0,       # 3. "Litchfield" 
                  0.04,    # 4. "Windham"  
                  0.05,    # 5. "Tolland" 
                  6.5,     # 6. "Hartford"
                  0.25,    # 7. "Middlesex" 
                  4.5)     # 8. "New Haven"

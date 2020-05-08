##########################################################
## Script to produce output for the dashboard
##########################################################
setwd("functions/")
source("../functions/global_var.R")

mydaymax              = ymd("2020-09-01") 
myschools_reopen_date = ymd("2021-09-01") 
mylockdown_end_date   = ymd("2020-06-01") 
mytesting_on_date     = mydaymax + 1
mydistancing_on_date  = mylockdown_end_date + 1 # distancing on at end of lockdown
mydistancing_stepdown_dates = seq(ymd("2020-07-01"), ymd("2021-04-01"), length.out=10)

nsim = 5

#########################################################
# First scenario
#########################################################
# get initial conditions
mystate0 = get_state0("../data/ct_init_a07.csv")
# get parameter values 
myparams = yaml.load_file("../parameters/params_a07.yml")  
# get sample from posterior
myposterior = read.csv("../data/posterior_a07.csv", stringsAsFactors=FALSE) 

myparams$testing_effect_A = 0.2
myparams$testing_effect_Im = 0.5

####################################
myparams$distancing_effect = 0.4
res1 = get_sir_results(daymax=mydaymax,
                      lockdown_end_date=mylockdown_end_date,
                      schools_reopen_date=myschools_reopen_date,
                      testing_on_date=mytesting_on_date,
                      distancing_on_date=mydistancing_on_date, 
                      distancing_stepdown_dates=mydistancing_stepdown_dates,
                      nsim=nsim,
                      params = myparams,
                      state0 = mystate0,
                      posterior = myposterior,
                      draw_rparams = FALSE )


#########################################################
# Second scenario
#########################################################
myparams$distancing_effect = 0.6
res2 = get_sir_results(daymax=mydaymax,
                      lockdown_end_date=mylockdown_end_date,
                      schools_reopen_date=myschools_reopen_date,
                      testing_on_date=mytesting_on_date,
                      distancing_on_date=mydistancing_on_date, 
                      distancing_stepdown_dates=mydistancing_stepdown_dates,
                      nsim=nsim,
                      params = myparams,
                      state0 = mystate0,
                      posterior = myposterior,
                      draw_rparams = FALSE )


#########################################################
# Third scenario
#########################################################
myparams$distancing_effect = 0.8
res3 = get_sir_results(daymax=mydaymax,
                      lockdown_end_date=mylockdown_end_date,
                      schools_reopen_date=myschools_reopen_date,
                      testing_on_date=mytesting_on_date,
                      distancing_on_date=mydistancing_on_date, 
                      distancing_stepdown_dates=mydistancing_stepdown_dates,
                      nsim=nsim,
                      params = myparams,
                      state0 = mystate0,
                      posterior = myposterior,
                      draw_rparams = FALSE )



#########################################################
# Output
#########################################################
# Organize output into list, names used as 'Scenario' variable in output
res <- list(`Scenario 1`   = res1,
			`Scenario 2` = res2, 
			`Scenario 3`    = res3)
# The get_dashboard_output function:
#	1. Returns the data frame
#	2. Show a plot for visual check
#	3. Saves CSV file
summary <- get_dashboard_output(data=res, plot=TRUE, 
								filename = paste0("../output/model output ", Sys.Date(), ".csv"))


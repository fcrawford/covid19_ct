##
## The working directory of this file is functions/ folder.
##

source("../functions/global_var.R")
source("../functions/run_ct_model.R")
source("../functions/plot_functions.R")

mydaymax              = ymd("2020-11-01") 

###########################
# intervention end dates
mylockdown_end_date = ymd("2021-09-01") # NEVER: taken care by random effects 
myschools_reopen_end_date = ymd("2021-09-01") # NEVER: taken care by random effects

# number of simulations 
# scenario to simulate: school_effect = 0.2 --> 20% increase
nsim = 2500
school_effect <- 0.2
####################################


## choose your scenario and get state0, params and a sample from posterior

# base scenario (low asymptomatic): q_A = 0.36, mean(q_Is) = 0.064
mystate0 = get_state0("../data/ct_init.csv")
myparams = yaml.load_file("../parameters/params.yml")  
myposterior = read.csv("../posterior/posterior_params.csv", stringsAsFactors=FALSE) 
myposterior.re = read.csv("../posterior/posterior_re.csv", stringsAsFactors=FALSE) 


## set value of school reopening effect
myparams$school_reopening_effect = log(1 + school_effect)

###########################
# intervention end dates and ramping

mytmax = as.numeric(difftime(mydaymax, day0, units="days"))
mydayseq = seq(day0, mydaymax, by="day")

myint_off_dates = as.list(c(lockdown_end_date = mylockdown_end_date, schools_reopen_end_date = myschools_reopen_end_date))

myint_ramping_times = as.list(c(14, 2))


distancingfun_list = get_distancing_fun_list(mydayseq, INT_START_DATES, myint_off_dates, myint_ramping_times)
mobilityfun = get_mobility_fun(mydayseq, MOB)
testingfun = get_testing_fun(mydayseq, TESTING)

# combine in the list of interventions 
INTERVENTIONS = list(distancing_list=distancingfun_list, mobility=mobilityfun, testing=testingfun) 


ptm = proc.time()
res1 = get_sir_results( daymax=mydaymax,
                        int_off_dates = myint_off_dates,
                        nsim=nsim,
                        params = myparams,
                        state0 = mystate0,
                        posterior = myposterior,
                        posterior.re = myposterior.re,
                        draw_rparams = FALSE )
print(proc.time() - ptm)



#### Plots ####

# w = 900
# h = 300
if(school_effect == 0){
    str = "\nNo change in contact rates starting September 1"
}else{
    str = paste0("\nContact rates increase by ", round(school_effect*100), "% starting September 1")
}


## plots: set 1 ##
pdf(paste0("../figures/test_scenario_", round(school_effect*100), ".pdf"), width=11.2, height=7.2)

par(mfrow=c(2, 2), mar=c(4,4,1,1))
##   Row 1
connecticut_summary_intx = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Overall community contact pattern in Connecticut", str),
                                             title=str,
                                             region_name=NULL, 
                                             which.plot="contact_pattern",
                                             ylab="",
                                             ymax=NULL)

connecticut_R_eff = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Community R effective in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="R_eff",
                                             ymax=2, ylab="Effective reproductive number")
abline(h=1, lty="dotted", col="gray")

##   Row 2
connecticut_summary_hosp1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("State-wide hospitalizations in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rHsum",
                                             add.cong = TRUE)

connecticut_summary_deaths1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("State-wide cumulative deaths in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rD",
                                             add.cong = TRUE)



par(mfrow=c(2, 2), mar=c(4,4,1,1))


## Row 3
connecticut_curentI = plot_ct_region(data=res1$summary, 
                                     end_day=mydaymax, 
                                     title.override = paste0("Current community infections in Connecticut", str),
                                     title=str,
                                     region_name="Connecticut", 
                                     which.plot="currentI",
                                     ymax=NULL, ylab="Current infections")

connecticut_summary_dailyI = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Daily new infections in community in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="dailyI",
                                             ymax=NULL)


##   Row 4
connecticut_cumulative_incid = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Cumulative incidence in community in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="alive_cum_incid_prop",
                                             ymax=NULL)

connecticut_cum_modI = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Cumulative infections in community in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="cum_modI",
                                             ymax=NULL)



dev.off()

# be aware of the date format
# get_snapshot(data=res1$summary, date="05/27/2020")
# get_snapshot(data=res1$summary, date="05/27/2020", where="New Haven")


source("global_var.R")

mydaymax              = ymd("2020-08-31") 

###########################
# intervention end dates
myschools_reopen_date = ymd("2020-09-01") # assuming schools reopen in September
mylockdown_end_date   = ymd("2021-09-01") # NEVER: replaced by phased lockdown lifting

# phased lockdown release end dates: never, but can be changed if release is rolled back
myph1_release_end_date = ymd("2021-09-01") 
myph2_release_end_date = ymd("2021-09-01")

myint_off_dates = as.list(c(schools_reopen_date = myschools_reopen_date, 
                      lockdown_end_date = mylockdown_end_date, 
                      ph1_release_end_date = myph1_release_end_date, 
                      ph2_release_end_date = myph2_release_end_date))



# number of simulations 
nsim = 1000

####################################

str = "\ntest"


## choose your scenario and get state0, params and a sample from posterior

# base scenario (low asymptomatic): q_A = 0.36, mean(q_Is) = 0.064
mystate0 = get_state0("../data/ct_init.csv")
myparams = yaml.load_file("../parameters/params.yml")  
myposterior = read.csv("../data/posterior_a036.csv", stringsAsFactors=FALSE) 
myposterior = post.save

####################################

ptm = proc.time()
# run simulation sampling parameters from joint posterior 
res1 = get_sir_results( daymax=mydaymax,
                        int_off_dates = myint_off_dates,
                        nsim=nsim,
                        params = myparams,
                        state0 = mystate0,
                        posterior = myposterior,
                        draw_rparams = FALSE )
proc.time() - ptm

####################################

w = 900
h = 300


par(mfrow=c(3,3), mar=c(4,4,1,1))

# plot_interventions(res1$raw_results, mydaymax, subtitle=NULL)

connecticut_summary_intx = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name=NULL, 
                                             which.plot="intervention_pattern",
                                             ylab="",
                                             ymax=NULL)


connecticut_summary_dailyI = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="dailyI",
                                             ymax=NULL)

connecticut_cum_modI = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="cum_modI",
                                             ymax=NULL)


connecticut_summary_hosp1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rHsum",
                                             ymax=NULL)

connecticut_summary_deaths1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rD",
                                             ymax=NULL)

connecticut_cumulative_incid_prop = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="alive_cum_incid_prop",
                                             ymax=NULL)




connecticut_cumulative_incid = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="alive_cum_incid_num",
                                             ymax=NULL)

connecticut_R_eff = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="R_eff",
                                             ymax=2, ylab="Effective reproductive number")
abline(h=1, lty="dotted", col="gray")


connecticut_curentI = plot_ct_region(data=res1$summary, 
                                     end_day=mydaymax, 
                                     title=str,
                                     region_name="Connecticut", 
                                     which.plot="currentI",
                                     ymax=NULL, ylab="Current infections")


# be aware of the date format
get_snapshot(data=res1$summary, date="05/27/2020")
get_snapshot(data=res1$summary, date="05/27/2020", where="New Haven")


source("global_var.R")

mydaymax              = ymd("2020-09-01") 
myschools_reopen_date = ymd("2021-09-01") 


nsim = 100

####################################

str = "\nwith release on 6/1 and phased reductions in distancing at businesses"
mylockdown_end_date   = ymd("2020-05-20") 
mytesting_on_date  = ymd("2020-05-20")
mydistancing_on_date  = mylockdown_end_date+1 # distancing on at end of lockdown
mydistancing_stepdown_dates = seq(ymd("2020-06-01"), ymd("2020-10-01"), length.out=10)


## choose your scenario

# get initial conditions
mystate0 = get_state0("../data/ct_init_a05.csv")

# get parameter values 
myparams = yaml.load_file("../parameters/params_a05.yml")  


#myparams$testing_effect_A = 0.2
#myparams$testing_effect_Im = 0.5
#myparams$distancing_effect = 0.6

####################################

res1 = get_sir_results(daymax=mydaymax,
                      lockdown_end_date=mylockdown_end_date,
                      schools_reopen_date=myschools_reopen_date,
                      testing_on_date=mytesting_on_date,
                      distancing_on_date=mydistancing_on_date, 
                      distancing_stepdown_dates=mydistancing_stepdown_dates,
                      nsim=nsim,
                      params = myparams,
                      state0 = mystate0)


####################################

w = 900
h = 300


par(mfrow=c(4,2), mar=c(4,4,1,1))



plot_interventions(res1$raw_results, mydaymax, subtitle=NULL)



connecticut_summary_dailyI = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="dailyI",
                                             ymax=NULL)



connecticut_summary_hosp1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="cum_modI",
                                             ymax=NULL)



ymax = 8e3
connecticut_summary_hosp1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rHsum",
                                             ymax=ymax)

connecticut_summary_deaths1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rD",
                                             ymax=NULL)

connecticut_cumulative_incid = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="alive_cum_incid_num",
                                             ymax=NULL)

connecticut_cumulative_incid_prop = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="alive_cum_incid_prop",
                                             ymax=NULL)

connecticut_R_eff = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="R_eff",
                                             ymax=2, ylab="Effective reproductive number")



abline(h=1, lty="dotted", col="gray")




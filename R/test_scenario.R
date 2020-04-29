

source("run_ct_model.R")

mydaymax              = ymd("2020-09-01") 
myschools_reopen_date = ymd("2021-09-01") 


nsim = 100

####################################

str = "\nwith release on 6/1 and phased reductions in distancing at businesses"
mylockdown_end_date   = ymd("2020-06-01") 
mytesting_on_date  = mydaymax
mydistancing_on_date  = mylockdown_end_date+1 # distancing on at end of lockdown
mydistancing_stepdown_dates = seq(ymd("2020-07-01"), ymd("2021-04-01"), length.out=10)

params_init$testing_effect_A = 0.2
params_init$testing_effect_Im = 0.5
params_init$distancing_effect = 0.6


####################################

res1 = get_sir_results(daymax=mydaymax,
                      lockdown_end_date=mylockdown_end_date,
                      schools_reopen_date=myschools_reopen_date,
                      testing_on_date=mytesting_on_date,
                      distancing_on_date=mydistancing_on_date, 
                      distancing_stepdown_dates=mydistancing_stepdown_dates,
                      nsim=nsim)


####################################

w = 900
h = 300

par(mfrow=c(5,1), mar=c(4,4,1,1))


ymax = 2e4
connecticut_summary_dailyI = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="dailyI",
                                             ymax=ymax)



ymax = 1e6
connecticut_summary_hosp1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="cum_modI",
                                             ymax=ymax)



ymax = 1e4
connecticut_summary_hosp1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rHsum",
                                             ymax=ymax)
ymax = 1e4
connecticut_summary_deaths1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rD",
                                             ymax=ymax)


plot_interventions(res1$raw_results, mydaymax, titles="")





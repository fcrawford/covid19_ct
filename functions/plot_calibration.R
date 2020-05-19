source("../functions/global_var.R")
mystate0 = get_state0("../data/ct_init.csv")
myparams = yaml.load_file("../parameters/params.yml")  
myposterior = read.csv("../data/posterior_a036.csv", stringsAsFactors=FALSE) 

mydaymax              = ymd("2020-09-01") 
myschools_reopen_date = ymd("3000-01-01")
mytesting_on_date     = ymd("2020-05-15")
mylockdown_end_date   = ymd("2020-06-01") 
mydistancing_stepdown_dates = seq(ymd("2020-06-15"), mydaymax, length.out=10)
nsim = 300
myseed = 12345
res1 = get_sir_results(daymax=mydaymax,
                      lockdown_end_date=mylockdown_end_date,
                      schools_reopen_date=myschools_reopen_date,
                      testing_on_date=mytesting_on_date,
                      distancing_on_date=mylockdown_end_date+1,
                      distancing_stepdown_dates=mydistancing_stepdown_dates,
                      params=myparams,
                      state0=mystate0,
                      nsim=nsim,
                      seed=myseed,
                      posterior = myposterior,
                      draw_rparams = FALSE )

pdf("../figures/calibration.pdf", width = 8, height = 3.2)
par(mfrow = c(1, 2))
out <-    plot_ct_region(data=res1$summary, region_name="Connecticut", which.plot="rD", title.override="Cumulative deaths", goodness=TRUE)
out <-    plot_ct_region(data=res1$summary, region_name="Connecticut", which.plot="rHsum", title.override="Required hospitalizations", goodness=TRUE)
dev.off()


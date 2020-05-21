source("global_var.R")

mydaymax              = ymd("2020-09-01") 
myschools_reopen_date = ymd("2021-09-01") 


nsim = 3000
# scenario of releasing contacts 
if(!"scenario_name" %in% ls()) scenario_name = c("slow", "fast")[2]
# # scenario of asymptomatic proportion
if(!"case_name" %in% ls())     case_name = c("low", "medium", "high")[2]

####################################

mylockdown_end_date   = ymd("2020-05-19") 
mytesting_on_date  = ymd("2020-05-20")

mydistancing_on_date  = mylockdown_end_date+1 # distancing on at end of lockdown


if(scenario_name == "fast") {
  str = "\nwith release on 5/20 and biweekly increments in contact"
  mydistancing_stepdown_dates = seq(mydistancing_on_date+1, mydistancing_on_date+1+days(14*10), by=14)
} else {
  str = "\nwith release on 5/20 and monthly increments in contact"
  mydistancing_stepdown_dates = seq(mydistancing_on_date+1, mydistancing_on_date+1+days(30*10), by=30)
}


## choose your scenario and get state0, params and a sample from posterior

# base scenario (low asymptomatic): q_A = 0.36, mean(q_Is) = 0.064
if(case_name == "low"){
    mystate0 = get_state0("../data/ct_init.csv")
    myparams = yaml.load_file("../parameters/params.yml")  
    myposterior = read.csv("../data/posterior_a036.csv", stringsAsFactors=FALSE)     
}

# alternative scenario 1 (medium asymptomatic): q_A = 0.5, mean(q_Is) = 0.05
if(case_name == "medium"){
    mystate0 = get_state0("../data/ct_init_a05.csv")
    myparams = yaml.load_file("../parameters/params_a05.yml")  
    myposterior = read.csv("../data/posterior_a05.csv", stringsAsFactors=FALSE) 
}
# alternative scenario 2 (high asymptomatic): q_A = 0.7, mean(q_Is) = 0.03
if(case_name == "high"){
    mystate0 = get_state0("../data/ct_init_a07.csv")
    myparams = yaml.load_file("../parameters/params_a07.yml")  
    myposterior = read.csv("../data/posterior_a07.csv", stringsAsFactors=FALSE)     
}


# set testing effects and distancing effect for nsim = 1
myparams$testing_effect_A = 0.2
myparams$testing_effect_Im = 0.5

myparams$distancing_effect = myparams$lockdown_effect

# set testing effects and distancing effect if drawing from posterior
myposterior$testing_effect_A = rtruncdist(nrow(myposterior), mean=0.2, sd=0.03, lower=0.1, upper=0.3)
myposterior$testing_effect_Im = rtruncdist(nrow(myposterior), mean=0.5, sd=0.04, lower=0.35, upper=0.65)

myposterior$distancing_effect = myparams$distancing_effect


####################################
if(file.exists(paste("../output/", scenario_name, "_", case_name, ".RData", sep=""))){
    load(paste("../output/", scenario_name, "_", case_name, ".RData", sep=""))
}else{
    ptm = proc.time()
    # run simulation sampling parameters from joint posterior 
    res1 = get_sir_results( daymax=mydaymax,
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
    print(proc.time() - ptm)
    save(res1, file = paste("../output/", scenario_name, "_", case_name, ".RData", sep=""))
}


####################################

w = 900
h = 300

pdf(paste(scenario_name, "_", case_name, ".pdf", sep=""), width=12, height=6)

par(mfrow=c(2,2), mar=c(4,4,1,1))

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
                                             which.plot="rHsum",
                                             ymax=NULL)

connecticut_summary_deaths1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rD",
                                             ymax=NULL)

dev.off()

#########################################


w = 900
h = 150

pdf(paste("cumincidence_", scenario_name, "_", case_name, ".pdf", sep=""), width=6, height=3)

par(mfrow=c(1,1), mar=c(4,4,1,1))

connecticut_cumulative_incid_prop = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="alive_cum_incid_prop",
                                             ymax=1, 
                                             ylab="Proportion")


dev.off()

#########################################

pdf(paste(scenario_name, "_", case_name, "_full.pdf", sep=""), width=6*2, height=3*3)

par(mfrow=c(3,2), mar=c(4,4,1,1))


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

connecticut_cumulative_incid = plot_ct_region(data=res1$summary, 
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
                                             ymax=NULL, ylab="Effective reproductive number")
abline(h=1, lty="dotted", col="gray")
dev.off()

#########################################

pdf(paste0("calibration_", scenario_name, "_", case_name, ".pdf"), width=12, height=4)
par(mfrow = c(1, 2))
out <-    plot_ct_region(data=res1$summary, region_name="Connecticut", which.plot="rHsum", title.override="Hospitalizations", goodness=TRUE)
out <-    plot_ct_region(data=res1$summary, region_name="Connecticut", which.plot="rD", title.override="Cumulative deaths", goodness=TRUE)
dev.off()







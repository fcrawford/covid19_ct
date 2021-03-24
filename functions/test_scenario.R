##
## The working directory of this file is functions/ folder.
##

source("../functions/global_var.R")
source("../functions/run_ct_model.R")
source("../functions/plot_functions.R")

if (FALSE){
   source("global_var.R")
   source("run_ct_model.R")
   source("plot_functions.R")
}


## simulate model projections through this date
mydaymax              = ymd("2021-04-30") 

##########################################
# intervention end dates for testing hypothetical interventions in the future
intv1_end_date = ymd("2022-12-31") # in the future


# number of simulations 
nsim = 1000
####################################


## get state0, params and a sample from posterior
mystate0 = get_state0("../data/ct_init.csv") # only used in nsim==1
myparams = yaml.load_file("../parameters/params.yml")  
if(FALSE){
   myposterior = read.csv("../posterior/posterior_params.csv", stringsAsFactors=FALSE) 
   myposterior.re = read.csv("../posterior/posterior_re.csv", stringsAsFactors=FALSE) 
}
myposterior = post.save
myposterior.re = post.save.re


# changes in parameters for this simulation
myparams$alpha_Is = 1/12

myparams$time_num = 0

myparams$H_lag = 5
myparams$D_lag = 5


###########################
# intervention end dates and ramping
mytmax = as.numeric(difftime(mydaymax, day0, units="days"))
mydayseq = seq(day0, mydaymax, by="day")

int_off_dates = as.list(c(intv1_end_date = intv1_end_date))
int_ramping_times = as.list(c(2))

distancingfun_list = get_distancing_fun_list(mydayseq, INT_START_DATES, int_off_dates, int_ramping_times)
mobilityfun = get_mobility_fun(mydayseq, MOB)
testingfun = get_testing_fun(mydayseq, TESTING)

# random effects
me_start = 9
random_effect_at <- seq(me_start, (length(dayseq)-28), by = 21) # random effects start on this day, and the first one is 0
random_effect <- rep(0, length(random_effect_at)) # these effects are calibrated to data, 0 is starting value

random_effect_fun = get_random_effect_fun(mydayseq, random_effect_at, random_effect)

# combine in the list of interventions 
INTERVENTIONS = list(distancing_list=distancingfun_list, random_effect = random_effect_fun, mobility=mobilityfun, testing=testingfun) 

ptm = proc.time()
res1 = get_sir_results( daymax=mydaymax,
                        int_off_dates = int_off_dates,
                        nsim=nsim,
                        params = myparams,
                        state0 = mystate0,
                        posterior = myposterior,
                        posterior.re = myposterior.re,
                        draw_rparams = FALSE, 
                        CI = 0.95)
print(proc.time() - ptm)



#### Plots ####

# w = 900
# h = 300
#if(school_effect == 0){
    str = " "
    # str = "\nunder observed contact rate increase in September"
#}else{
#    str = paste0("\nContact rates increase by ", round(school_effect*100), "% starting September 15")
#}

### only include important plots ###

png(paste0("../figures/proj.png"), width=2400, height=1200, res=150)

par(mfrow=c(2, 2), mar=c(4,4,1,1))
##   Row 1
connecticut_summary_hosp1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("State-wide hospitalizations in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rHsum",
                                             avg='median',
                                             ymax=3000,
                                             add.cong = TRUE)

connecticut_summary_deaths1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("State-wide cumulative deaths in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rD",
                                             avg='median',
                                             add.cong = TRUE)

##   Row 2
connecticut_summary_intx = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Overall community contact pattern in Connecticut", str),
                                             title=str,
                                             region_name=NULL, 
                                             which.plot="contact_pattern",
                                             ymax=0.5,
                                             ylab="")

connecticut_cumulative_incid = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Cumulative incidence in community in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="alive_cum_incid_prop",
                                             avg = 'median',
                                             ymax=0.5, ylab="Proportion")

dev.off()





##### single plot ######
png(paste0("../figures/proj_hosp.png"), width=1200, height=600, res=150)

par(mfrow=c(1, 1), mar=c(4,4,1,1))

connecticut_summary_hosp1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("State-wide hospitalizations in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rHsum",
                                             avg='mean',
                                             ymax=3000,
                                             add.cong = TRUE)
dev.off()




#### R_eff #####
connecticut_R_eff = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Community R effective in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="R_eff",
                                             avg="mean",
                                             ymax=2, ylab="Effective reproductive number")
abline(h=1, lty="dotted", col="gray")








## other plots: set 1 ##
pdf(paste0("../figures/proj_school_effect_", round(school_effect*100), ".pdf"), width=11.2, height=7.2)

png(paste0("../figures/proj_p1.png"), width=1800, height=1200, res=150)

par(mfrow=c(2, 2), mar=c(4,4,1,1))
##   Row 1
connecticut_summary_intx = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Overall community contact pattern in Connecticut", str),
                                             title=str,
                                             region_name=NULL, 
                                             which.plot="contact_pattern",
                                             ylab="")

connecticut_R_eff = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Community R effective in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="R_eff",
                                             avg="median",
                                             ymax=2, ylab="Effective reproductive number")
abline(h=1, lty="dotted", col="gray")

##   Row 2
connecticut_summary_hosp1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("State-wide hospitalizations in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rHsum",
                                             avg='median',
                                             ymax=3000,
                                             add.cong = TRUE)

connecticut_summary_deaths1 = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("State-wide cumulative deaths in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="rD",
                                             add.cong = TRUE)

dev.off()



png(paste0("../figures/proj_p2.png"), width=1800, height=1200, res=140)
par(mfrow=c(2, 2), mar=c(4,4,1,1))

## Row 3
connecticut_curentI = plot_ct_region(data=res1$summary, 
                                     end_day=mydaymax, 
                                     title.override = paste0("Current community infections in Connecticut", str),
                                     title=str,
                                     region_name="Connecticut", 
                                     which.plot="currentI",
                                     avg='median',
                                     ymax=300000, ylab="Current infections")

connecticut_summary_dailyI = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Daily new infections in community in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="dailyI",
                                             avg='median',
                                             ymax=30000)


##   Row 4
connecticut_cumulative_incid = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Cumulative incidence in community in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="alive_cum_incid_prop",
                                             avg = 'median',
                                             ymax=0.5)

connecticut_cum_modI = plot_ct_region(data=res1$summary, 
                                             end_day=mydaymax, 
                                             title.override = paste0("Cumulative infections in community in Connecticut", str),
                                             title=str,
                                             region_name="Connecticut", 
                                             which.plot="cum_modI",
                                             avg = 'median',
                                             ymax=1500000)

dev.off()

# be aware of the date format
# get_snapshot(data=res1$summary, date="05/27/2020")
# get_snapshot(data=res1$summary, date="05/27/2020", where="New Haven")


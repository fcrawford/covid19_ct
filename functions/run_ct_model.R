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
#####################

source("../functions/model.R")
source("../functions/get_data.R")
source("../functions/intervention_functions.R")
source("../functions/truncated_distributions.R")
source("../functions/plot_functions.R")

########################
# Dates

# starting day
day0 = ymd("2020-03-01")

# actual dates: do not change
state_schools_close = dmy("13/03/2020")
state_lockdown_order = dmy("20/03/2020") # order date
state_lockdown_start = dmy("23/03/2020") # actual start date


#############################
# hypothetical end date of lockdown

# lockdown_end_date=dmy("01/06/2020") 


########################

##### get reported data from CT #####

data <- get_ct_data(day0=day0)
dat_ct_state <- data$dat_ct_state
dat_ct_county <- data$dat_ct_county


############ regional population and adjacency #############

# load region adjacency matrix 
adj = read.csv("../map/CT_adj_matrix.csv", stringsAsFactors=FALSE)
rownames(adj) = adj$X
adj = adj[,-1]
adj = as.matrix(adj)

# county populations
nregions = nrow(adj)

pop = read.csv('../data/ct_population.csv', stringsAsFactors=FALSE)
populations = list()
populations[1:nregions] = pop$population
names(populations) = pop$county

region_names = pop$county


#######################

# set county level capacities 

county_capacities = list()
for(nm in region_names) {
  county_cap = filter(dat_ct_capacity, County==nm)
  county_days = as.numeric(ymd(county_cap$Date) - day0)
  county_cap_fun = approxfun(county_days, county_cap$Capacity, rule=2)
  county_capacities[[nm]] = county_cap_fun
}


# the ordering of counties is the standard throughout the code
# make sure adj and initial conditions have the same ordering! 

########################
# draw random params

rparams = function(params) {
  params_tmp = params
  # sample new param values
  params_tmp$beta_pre = rtruncdist(1, mean=(params$beta_pre*0.9975), sd=params$sd_beta_pre, lower=params$lower_beta_pre, upper=params$upper_beta_pre)
  params_tmp$q_Im = rtruncdist(1, mean=(params$q_Im), sd=params$sd_q_Im, lower=params$lower_q_Im, upper=params$upper_q_Im)
  params_tmp$gamma_H = rtruncdist(1, mean=params$gamma_H, sd=params$sd_gamma_H, lower=params$lower_gamma_H, upper=params$upper_gamma_H)
  params_tmp$m_H = rtruncdist(1, mean=params$m_H, sd=params$sd_m_H, lower=params$lower_m_H, upper=params$upper_m_H)
  params_tmp$m_Hbar_mult = rtruncdist(1, mean=params$m_Hbar_mult, sd=params$sd_m_Hbar_mult, lower=params$lower_m_Hbar_mult, upper=params$upper_m_Hbar_mult)
  params_tmp$lockdown_effect = rtruncdist(1, mean=params$lockdown_effect, sd=params$sd_lockdown_effect, lower=params$lower_lockdown_effect, upper=params$upper_lockdown_effect)
  #params_tmp$q_A = rtruncdist(1, mean=(params$q_A), sd=params$sd_q_A, lower=params$lower_q_A, upper=params$upper_q_A)
  # params_tmp$delta = rtruncdist(1, mean=params$delta, sd=params$sd_delta, lower=params$lower_delta, upper=params$upper_delta)
  return(params_tmp)
}



########################
# function: get state0 from initial confitions
get_state0 = function(init_file_csv){
init <- read.csv(init_file_csv, stringsAsFactors=FALSE) 

      E_init = init$E
      I_s_init = init$Is
      I_m_init = init$Im
      A_init = init$A
      H_init = init$H
      Hbar_init = rep(0,nregions)
      D_init = init$D
      R_init = init$R
      S_init = as.numeric(populations) - (E_init + I_s_init + I_m_init + A_init + H_init + Hbar_init + D_init + R_init)

# this is state0 to be passed to get_sir_results 
state0 = c(S=S_init, E=E_init, I_s=I_s_init, I_m=I_m_init, A=A_init, H=H_init, Hbar=Hbar_init, D=D_init, R=R_init)
return(state0)
}


########################
#### set default parameters and initial conditions
# these can be changed when the model run is called

# get initial conditions
#mystate0 = get_state0("../data/ct_init.csv")

# default parameters
#myparams = yaml.load_file("../parameters/params.yml")  



#######################
# Run the model 


get_sir_results = function(daymax=ymd("2020-09-01"), 
                           lockdown_end_date, 
                           schools_reopen_date,
                           testing_on_date,
                           distancing_on_date,
                           distancing_stepdown_dates,
                           nsim=1,
                           params = myparams,
                           state0 = mystate0,
                           seed = NULL) {
  # parameters, set seed if given
  if(!is.null(seed)) set.seed(seed)
  pars <- list()
  if(nsim == 1){
     pars[[1]] =  myparams 
  } else { 
    for(i in 1:nsim) pars[[i]] <- rparams(myparams)
  } 
   

  dayseq = seq(day0, daymax, by="day")
  tmax = as.numeric(difftime(daymax, day0, units="days"))

  lockfun = get_state_lockdown_fun(dayseq, offdate=lockdown_end_date)
  schoolsfun = get_school_in_session_fun(dayseq, schools_reopen_date=schools_reopen_date)
  testingfun = get_testing_on_fun(dayseq, testing_on_date=testing_on_date)
  distancingfun = get_distancing_stepdown_fun(dayseq, distancing_on_date=distancing_on_date, distancing_stepdown_dates=distancing_stepdown_dates)

  interventions = list(lockdown=lockfun, schools=schoolsfun, testing=testingfun, distancing=distancingfun) 

  sir_results = lapply(1:nsim, function(i){
    res = run_sir_model(state0=state0, 
                        params=pars[[i]],  # note: effect_intvx is in params, and is not passed to run_sir_model separately 
                        region_adj=adj, 
                        populations=as.numeric(populations), 
                        tmax=tmax, 
                        interventions=interventions,
                        capacities=county_capacities)
    res$sim_id = i
    res
  })

  sir_results_all = ldply(sir_results, rbind)
  for(nm in c("Connecticut", region_names)){
    sir_results_all[, paste0("rHsum.", nm)] <-  sir_results_all[,paste0("rH.", nm)]+sir_results_all[,paste0("rHbar.", nm)]
  }
  sir_results_long <- melt(sir_results_all, id.vars = c("time", "sim_id"))
  sir_results_summary <- sir_results_long %>% group_by(variable, time) %>% 
			                     summarise(
                             mean = mean(value),
			                       lower = quantile(value, 0.05, na.rm=TRUE),
			                       upper = quantile(value, 0.95, na.rm=TRUE))
  return(list(raw_results=sir_results, summary=sir_results_summary))
}


plot_interventions = function(sir_results, daymax, stayhome_compares=FALSE, titles=NULL) {

  if(stayhome_compares){
    sir_results_full <- sir_results
    sir_results <- sir_results[[1]]
  }

  dayseq = seq(day0, daymax, by="day")
  tmax = as.numeric(difftime(daymax, day0, units="days"))

  monthseq = seq(day0, daymax, by="month")
  monthseq_lab = format(monthseq, "%b %Y")
  daymonthseq = difftime(monthseq, day0, units="days")

  #print(monthseq)
  #lockfun = get_state_lockdown_fun(dayseq, offdate=lockdown_end_date)
  #schoolsfun = get_school_in_session_fun(dayseq, schools_reopen_date=schools_reopen_date)
  #testingfun = get_testing_on_fun(dayseq, testing_on_date=testing_on_date)
  #interventions = list(lockdown=lockfun, schools=schoolsfun, testing=testingfun) 
  #stop("here")

  par(mar=c(3,3,3,0), bty="n")
  nn <- 4
  if(stayhome_compares) nn <- 3 + length(sir_results_full)
  #layout(matrix(c(1:nn),nrow=,nn), heights=rep(2, nn))

  plot(sir_results[[1]]$intervention_pattern, ylim=c(0,1), xlim=c(0,1.05*tmax), type="n", ylab="", xlab="", main="Overall contact intervention", axes=FALSE)
  axis(1, at=daymonthseq, lab=monthseq_lab)
  axis(2)
  polygon(c(1,1:(tmax+1), tmax+1), c(0,sir_results[[1]]$intervention_pattern, 0), col="orange", border=NA)
  abline(v=Sys.Date()-day0, col="gray", lty=2)

  #if(!stayhome_compares){
      #plot(sir_results[[1]]$intervention_lockdown, ylim=c(0,1), type="n", ylab="", xlab="", main="Stay-at-home order in place", axes=FALSE)
      #axis(1, at=daymonthseq, lab=monthseq_lab)
      #axis(2)
      #polygon(c(1,1:(tmax+1), tmax+1), c(0,sir_results[[1]]$intervention_lockdown, 0), col="orange", border=NA)
      #abline(v=Sys.Date()-day0, col="gray", lty=2)
  #}


  #plot(sir_results[[1]]$intervention_schools, ylim=c(0,1), type="n", ylab="", xlab="", main="Schools in session", axes=FALSE)
  #axis(1, at=daymonthseq, lab=monthseq_lab)
  #axis(2)
  #polygon(c(1,1:(tmax+1), tmax+1), c(0,sir_results[[1]]$intervention_schools, 0), col="orange", border=NA)
  #abline(v=Sys.Date()-day0, col="gray", lty=2)
  
  #if(!stayhome_compares){
      #plot(sir_results[[1]]$intervention_lockdown, ylim=c(0,1), type="n", ylab="", xlab="", main="Stay-at-home order in place", axes=FALSE)
      #axis(1, at=daymonthseq, lab=monthseq_lab)
      #axis(2)
      #polygon(c(1,1:(tmax+1), tmax+1), c(0,sir_results[[1]]$intervention_lockdown, 0), col="orange", border=NA)
      #abline(v=Sys.Date()-day0, col="gray", lty=2)
  #}

  #plot(sir_results[[1]]$intervention_pattern, ylim=c(0,1),  type="n", ylab="", xlab="", main="Relative reduction in transmission", axes=FALSE)
  #axis(1, at=daymonthseq, lab=monthseq_lab)
  #axis(2)
  #polygon(c(1,1:(tmax+1), tmax+1), c(0,sir_results[[1]]$intervention_pattern, 0), col="orange", border=NA)
  #abline(v=Sys.Date()-day0, col="gray", lty=2)

  #plot(sir_results[[1]]$intervention_testing, ylim=c(0,1), type="n", ylab="", xlab="", main="Expanded testing", axes=FALSE)
  #axis(1, at=daymonthseq, lab=monthseq_lab)
  #axis(2)
  #polygon(c(1,1:(tmax+1), tmax+1), c(0,sir_results[[1]]$intervention_testing, 0), col="orange", border=NA)
  #abline(v=Sys.Date()-day0, col="gray", lty=2)

  #if(stayhome_compares){
      #for(i in 1:length(sir_results_full)){
        #plot(sir_results_full[[i]][[1]]$intervention_lockdown, ylim=c(0,1), type="n", ylab="", xlab="", main=titles[[i]], axes=FALSE)
        #axis(1, at=daymonthseq, lab=monthseq_lab)
        #axis(2)
        #polygon(c(1,1:(tmax+1), tmax+1), c(0,sir_results_full[[i]][[1]]$intervention_lockdown, 0), col="orange", border=NA)
        #abline(v=Sys.Date()-day0, col="gray", lty=2)
      #}
  #}

}

####################################
# Print a vector of counts in the RMD file
# get_compartment(date="2020-07-01", toprint="rD.Connecticut")
get_compartment <- function(data=sir_results_summary, date, toprint, start_day = day0){
  sir_result_internal = data.frame(filter(data, variable%in%toprint))
  t = as.numeric(difftime(as.Date(date), start_day, unit='days'))
  sir_result_internal = subset(sir_result_internal, time%in%t)
  out  <- as.character(format(sir_result_internal$mean, digits=2, big.mark=","))
  if(length(date)>2){
    out <-  paste0(paste(out[-length(out)], collapse=", "), ", and ", out[length(out)])
  }else if(length(date)==2){
    out <-  paste0(out[1], " and ", out[2])
  }
  return(out)
}









# the ordering of counties is the standard throughout the code (in the global_var.R)
# make sure adj and initial conditions have the same ordering! 

########################
# draw random params

rparams = function(params) {
  params_tmp = params
  # sample new param values
  params_tmp$beta_pre = rtruncdist(1, mean=(params$beta_pre*0.9975), sd=params$sd_beta_pre, lower=params$lower_beta_pre, upper=params$upper_beta_pre)
  params_tmp$q_Is = rtruncdist(1, mean=(params$q_Is), sd=params$sd_q_Is, lower=params$lower_q_Is, upper=params$upper_q_Is)
  params_tmp$gamma_H = rtruncdist(1, mean=params$gamma_H, sd=params$sd_gamma_H, lower=params$lower_gamma_H, upper=params$upper_gamma_H)
  params_tmp$gamma_Hbar = params_tmp$gamma_H
  params_tmp$m_H = rtruncdist(1, mean=params$m_H, sd=params$sd_m_H, lower=params$lower_m_H, upper=params$upper_m_H)
  params_tmp$m_Hbar_mult = rtruncdist(1, mean=params$m_Hbar_mult, sd=params$sd_m_Hbar_mult, lower=params$lower_m_Hbar_mult, upper=params$upper_m_Hbar_mult)
  params_tmp$lockdown_effect = rtruncdist(1, mean=params$lockdown_effect, sd=params$sd_lockdown_effect, lower=params$lower_lockdown_effect, upper=params$upper_lockdown_effect)
  # params_tmp$q_Im = rtruncdist(1, mean=(params$q_Im), sd=params$sd_q_Im, lower=params$lower_q_Im, upper=params$upper_q_Im)
  # params_tmp$q_A = rtruncdist(1, mean=(params$q_A), sd=params$sd_q_A, lower=params$lower_q_A, upper=params$upper_q_A)
  # params_tmp$delta = rtruncdist(1, mean=params$delta, sd=params$sd_delta, lower=params$lower_delta, upper=params$upper_delta)
  return(params_tmp)
}



########################
# function: get state0 from initial confitions
get_state0 = function(init_file_csv){
      init <- read.csv(init_file_csv, stringsAsFactors=FALSE) 
      # double double check order
      init <- init[match(CT_NAMES,  init$county), ]

      E_init = init$E
      I_s_init = init$Is
      I_m_init = init$Im
      A_init = init$A
      H_init = init$H
      Hbar_init = rep(0,dim(init)[1])
      NH_init = init$NH
      NI_init = init$NI
      D_init = init$D
      R_init = init$R
      S_init = as.numeric(init$population) - (E_init + I_s_init + I_m_init + A_init + H_init + Hbar_init + NH_init + NI_init + D_init + R_init)

       # this is state0 to be passed to get_sir_results 
      state0 = c(S=S_init, E=E_init, I_s=I_s_init, I_m=I_m_init, A=A_init, H=H_init, Hbar=Hbar_init, NH=NH_init, NI=NI_init, D=D_init, R=R_init)
      return(state0)
}










########################
## draw params from joint posterior 

# some supporting functions

# initial numbers exposed in each county: used to generate state0
# E_init_state0 is a global variable

# get state0 for a given set of parameters
get_state0_params <- function(params, E_init_state0, populations, adj, county_capacities){

nregions = length(E_init_state0)
dmax = day0 + params$time_num + 10
dayseq = seq(day0, dmax, by="day")
tmax = params$time_num + 1

params_tmp <- params
params_tmp$school_closure_effect <- 0
params_tmp$lockdown_effect <- 0
params_tmp$distancing_effect <- 0
params_tmp$testing_effect_Im <- 0
params_tmp$testing_effect_A <- 0
params_tmp$H_lag <- 0
params_tmp$D_lag <- 0

lockdown_end_date <- ymd("2020-06-01")
schools_reopen_date <- ymd("3000-01-01") # never
testing_on_date <- ymd("2020-06-01")
distancing_on_date  = lockdown_end_date + 1 # distancing on at end of lockdown
distancing_stepdown_dates = seq(ymd(distancing_on_date+1), ymd(ymd("2020-06-01")+30), length.out=2)

lockfun = get_state_lockdown_fun(dayseq, offdate=lockdown_end_date)
schoolsfun = get_school_in_session_fun(dayseq, schools_reopen_date=schools_reopen_date)
testingfun = get_testing_on_fun(dayseq, testing_on_date=testing_on_date)
distancingfun = get_distancing_stepdown_fun(dayseq, distancing_on_date=distancing_on_date, distancing_stepdown_dates=distancing_stepdown_dates)

interventions = list(lockdown=lockfun, schools=schoolsfun, testing=testingfun, distancing=distancingfun) 

E_init = E_init_state0
I_s_init = rep(0,nregions)
I_m_init = rep(0,nregions)
A_init = rep(0,nregions)
H_init = rep(0,nregions)
Hbar_init = rep(0,nregions)
NH_init = rep(0,nregions)
NI_init = rep(0,nregions)
D_init = rep(0,nregions)
R_init = rep(0,nregions)
S_init = as.numeric(populations) - (E_init + I_s_init + I_m_init + A_init + H_init + Hbar_init + NH_init + NI_init + D_init + R_init)

state0 = c(S=S_init, E=E_init, I_s=I_s_init, I_m=I_m_init, A=A_init, H=H_init, Hbar=Hbar_init, NH=NH_init, NI=NI_init, D=D_init, R=R_init)

res = run_sir_model(state0=state0, 
                    params=params_tmp,  
                    region_adj=adj, 
                    populations=as.numeric(populations), 
                    tmax=tmax, 
                    interventions=interventions, 
                    capacities=county_capacities)

init = matrix(0, ncol=10, nrow=8)
compartments <- c("E", "I_s", "I_m", "A", "H", "NH", "NI", "D", "R" )
  
for (k in 1:length(compartments)){
  compartment = compartments[k]
  comp.init <- c(nregions)
  for (i in 1:nregions){  
    region = CT_NAMES[i]
    idx = paste(compartment, region, sep=".")
    comp.init[i] = res[params$time_num, idx]
  }
 init[,k] = comp.init 
}
  init = as.data.frame(init)
  colnames(init) <- c("E", "Is", "Im", "A", "H", "NH", "NI", "D", "R", "county")
  init$county <- CT_NAMES

      E_init = init$E
      I_s_init = init$Is
      I_m_init = init$Im
      A_init = init$A
      H_init = init$H
      Hbar_init = rep(0,dim(init)[1])
      NH_init = init$NH
      NI_init = init$NI
      D_init = init$D
      R_init = init$R
      S_init = as.numeric(populations) - (E_init + I_s_init + I_m_init + A_init + H_init + Hbar_init + NH_init + NI_init + D_init + R_init)

       # this is state0 to be passed to get_sir_results 
      state0 = c(S=S_init, E=E_init, I_s=I_s_init, I_m=I_m_init, A=A_init, H=H_init, Hbar=Hbar_init, NH=NH_init, NI=NI_init, D=D_init, R=R_init)
      return(state0)
}

 



# return a list of params and state0 for a random draw from joint posterior
## inputs:
# params: list of parameters representing a chosen scenario
# posterior: data frame with a sample from joint posterior 
rposterior = function(params, posterior){

par_smpl = posterior[sample(c(1:nrow(posterior)), size = 1), ]

params$q_Is = par_smpl$q_Is
params$q_Im = 1 - params$q_A - params$q_Is
params$q_H = par_smpl$q_H

params$beta_pre = par_smpl$beta_pre 
params$k_A = par_smpl$k_A
params$k_Is = par_smpl$k_Is

params$gamma_H = par_smpl$gamma_H
params$gamma_Hbar = par_smpl$gamma_H
params$gamma_NH = par_smpl$gamma_NH
   
params$m_H = par_smpl$m_H
params$m_NH_mult = par_smpl$m_NH_mult
params$m_Hbar_mult = par_smpl$m_Hbar_mult
   
params$H_lag = par_smpl$H_lag
params$D_lag = par_smpl$D_lag
  
params$school_closure_effect = par_smpl$school_closure_effect
params$lockdown_effect = par_smpl$lockdown_effect

params$time_num = par_smpl$time_num

state0 = get_state0_params(params, E_init_state0, CT_POPULATIONS, CT_ADJ, COUNTY_CAPACITIES)

return( list (params, state0) )
}








#######################
## Run the model 
# draw_rparams controls uncertainty simulation: 
# if draw_rparams = TRUE, parameters are drawn independently using rparams()
# if draw_rparams = TRUE, parameters are drawn from joint posterior supplied as 'posterior'

get_sir_results = function(daymax, 
                           lockdown_end_date, 
                           schools_reopen_date,
                           testing_on_date,
                           distancing_on_date,
                           distancing_stepdown_dates,
                           nsim=1,
                           params,
                           state0,
                           posterior,
                           draw_rparams = FALSE,
                           seed = NULL) {

 pars <- list()
 state0s <- list()
   
# parameters, set seed if given
  if(!is.null(seed)) set.seed(seed)
  
  if (draw_rparams == TRUE){
       if(nsim == 1){
     pars[[1]] = params 
     state0s[[1]] = state0
  } else { 
    for(i in 1:nsim){ 
       pars[[i]] = rparams(params)
       state0s[[i]] = state0
    }
  } 
} else {
  if(nsim == 1){
     pars[[1]] = params
     state0s[[1]] = state0
  } else { 
    for(i in 1:nsim){ 
       rpost_out = rposterior(params, posterior)
       pars[[i]] = rpost_out[[1]]
       state0s[[i]] = rpost_out[[2]]
    }
  } 
}   

  dayseq = seq(day0, daymax, by="day")
  tmax = as.numeric(difftime(daymax, day0, units="days"))

  lockfun = get_state_lockdown_fun(dayseq, offdate=lockdown_end_date)
  schoolsfun = get_school_in_session_fun(dayseq, schools_reopen_date=schools_reopen_date)
  testingfun = get_testing_on_fun(dayseq, testing_on_date=testing_on_date)
  distancingfun = get_distancing_stepdown_fun(dayseq, distancing_on_date=distancing_on_date, distancing_stepdown_dates=distancing_stepdown_dates)

  interventions = list(lockdown=lockfun, schools=schoolsfun, testing=testingfun, distancing=distancingfun) 

  sir_results = lapply(1:nsim, function(i){
    res = run_sir_model(state0=state0s[[i]], 
                        params=pars[[i]],  # note: effect_intvx is in params, and is not passed to run_sir_model separately 
                        region_adj=CT_ADJ, 
                        populations=CT_POPULATIONS, 
                        tmax=tmax, 
                        interventions=interventions,
                        capacities=COUNTY_CAPACITIES)
    res$sim_id = i
    res
  })

  sir_results_all = ldply(sir_results, rbind)
  for(nm in c("Connecticut", colnames(CT_ADJ))){
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














####################################
# Print a vector of counts in the RMD file
# get_compartment(data=res1$sir_results_summary, date="2020-07-01", toprint="rD.Connecticut",start_day = day0)
get_compartment <- function(data, date, toprint, start_day){
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


get_dashboard_output <- function(data, plot=TRUE, filename=NULL){
  toprint <- c("rD.Connecticut", "rHsum.Connecticut", "dailyI.Connecticut", "dailyH.Connecticut")  
  if(is.null(names(data))) names(data) <- paste("Scenario", 1:length(data))
  out <- NULL
  for(i in 1:length(data)){
      sir_result_internal <- filter(data[[i]]$summary, variable%in%toprint) %>% 
                             ungroup(variable) %>% 
                             mutate(variable=revalue(variable, c("rD.Connecticut"="Projected_Cumulative_Deaths",
                                                                 "rHsum.Connecticut"="Projected_Hospitalizations",
                                                                 "dailyI.Connecticut"="Projected_New_Infections",
                                                                 "dailyH.Connecticut"="Projected_New_Hospitalizations"))) %>%
                             mutate(Date = day0 + time, value = mean, Scenario=names(data)[i]) %>%
                             select(Date, variable, value, Scenario) 
       out <- rbind(out, sir_result_internal)
  }
  out$Scenario <- factor(out$Scenario,  levels = names(data))

  obs <- DAT_CT_STATE[, c("date", "deaths", "cur_hosp")]
  colnames(obs) <- c("Date", "Actual_Cumulative_Deaths", "Actual_Hospitalizations")
  # add back 3/1 to 3/8
  zeros <- data.frame(Date = day0 + c(0:(obs$Date[1] - day0 + 1)))
  zeros$Actual_Cumulative_Deaths <- 0
  zeros$Actual_Hospitalizations <- 0
  obs <- rbind(zeros, obs)

  out <- dcast(out, Date + Scenario ~ variable,  value.var = "value")
  out <- out %>% full_join(obs,  by = "Date") %>% arrange(Scenario, Date)
  if(plot){
    g <- ggplot(melt(out, id.vars=c("Date", "Scenario")), aes(x=Date, y=value, color=variable, group=variable)) + geom_line() + facet_wrap(~Scenario) + theme(legend.position = "bottom")
    print(g)
  }
  out$Produced_by <- "Crawford Lab"
  out$Produced_time <- Sys.time()
  if(!is.null(filename)){
    write.csv(out, file = filename, row.names=FALSE, quote=FALSE)
  }
  return(out)
}






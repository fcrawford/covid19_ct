# the ordering of counties is the standard throughout the code (in the global_var.R)
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
      Hbar_init = rep(0,dim(init)[1])
      D_init = init$D
      R_init = init$R
      S_init = as.numeric(init$population) - (E_init + I_s_init + I_m_init + A_init + H_init + Hbar_init + D_init + R_init)

       # this is state0 to be passed to get_sir_results 
      state0 = c(S=S_init, E=E_init, I_s=I_s_init, I_m=I_m_init, A=A_init, H=H_init, Hbar=Hbar_init, D=D_init, R=R_init)
      return(state0)
}

#######################
# Run the model 

get_sir_results = function(daymax, 
                           lockdown_end_date, 
                           schools_reopen_date,
                           testing_on_date,
                           distancing_on_date,
                           distancing_stepdown_dates,
                           nsim=1,
                           params,
                           state0,
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
                        populations=populations, 
                        tmax=tmax, 
                        interventions=interventions,
                        capacities=county_capacities)
    res$sim_id = i
    res
  })

  sir_results_all = ldply(sir_results, rbind)
  for(nm in c("Connecticut", colnames(adj))){
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









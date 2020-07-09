# the ordering of counties is the standard throughout the code (in the global_var.R)
# make sure adj and initial conditions have the same ordering! 

########################
# draw random params

rparams = function(params) {
  params_tmp = params
  # sample new param values
  params_tmp$beta_pre = rtruncdist(1, mean=(params$beta_pre*0.9975), sd=params$sd_beta_pre, lower=params$lower_beta_pre, upper=params$upper_beta_pre)
  #params_tmp$q_Is = rtruncdist(1, mean=(params$q_Is), sd=params$sd_q_Is, lower=params$lower_q_Is, upper=params$upper_q_Is)
  #params_tmp$gamma_H = rtruncdist(1, mean=params$gamma_H, sd=params$sd_gamma_H, lower=params$lower_gamma_H, upper=params$upper_gamma_H)
  #params_tmp$gamma_Hbar = params_tmp$gamma_H
  #params_tmp$m_H = rtruncdist(1, mean=params$m_H, sd=params$sd_m_H, lower=params$lower_m_H, upper=params$upper_m_H)
  #params_tmp$lockdown_effect = rtruncdist(1, mean=params$lockdown_effect, sd=params$sd_lockdown_effect, lower=params$lower_lockdown_effect, upper=params$upper_lockdown_effect)
  # post-lockdown
  #params_tmp$distancing_effect = rtruncdist(1, mean=params$distancing_effect, sd=params$sd_distancing_effect, lower=params$lower_distancing_effect, upper=params$upper_distancing_effect)
  #params_tmp$testing_effect_Im = rtruncdist(1, mean=params$testing_effect_Im, sd=params$sd_testing_effect_Im, lower=params$lower_testing_effect_Im, upper=params$upper_testing_effect_Im)
  #params_tmp$testing_effect_A = rtruncdist(1, mean=params$testing_effect_A, sd=params$sd_testing_effect_A, lower=params$lower_testing_effect_A, upper=params$upper_testing_effect_A)
  #params_tmp$m_Hbar_mult = rtruncdist(1, mean=params$m_Hbar_mult, sd=params$sd_m_Hbar_mult, lower=params$lower_m_Hbar_mult, upper=params$upper_m_Hbar_mult)
  # params_tmp$q_Im = rtruncdist(1, mean=(params$q_Im), sd=params$sd_q_Im, lower=params$lower_q_Im, upper=params$upper_q_Im)
  # params_tmp$q_A = rtruncdist(1, mean=(params$q_A), sd=params$sd_q_A, lower=params$lower_q_A, upper=params$upper_q_A)
  # params_tmp$delta = rtruncdist(1, mean=params$delta, sd=params$sd_delta, lower=params$lower_delta, upper=params$upper_delta)
  return(params_tmp)
}




# UPDATE THIS FUNCTION AS WE INCREASE THE NUMBER OF INTERVENTIONS
# create a list of contact intervention effects from params
get_intervention_effects = function(params) {
   
   return( as.list(c(params$school_closure_effect, 
                     params$lockdown_effect, 
                     params$ph1_release_effect, 
                     params$ph2_release_effect)) )
}


# create a list of intervention ramping times from params
get_ramping_times = function(params) {
   
   return( as.list(c(1,  # ramping time for school closure
                     params$ramp_time, # ramping time for lockdown
                     round(params$ramp_time/2, 0),  # ramping time for phase 1 release
                     round(params$ramp_time/2, 0)) ) ) # ramping time for phase 2 release
}





########################
# function: get state0 from initial confitions csv file
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










###################################
## draw params from joint posterior 
###################################
# some supporting functions

# get state0 on day0 for a given set of parameters
get_state0_params <- function(params, E_init_state0, interventions, populations, adj, county_capacities){

nregions = length(E_init_state0)
mytmax = params$time_num + 2

params_tmp <- params

params_tmp$school_closure_effect <- 0
params_tmp$lockdown_effect <- 0

params_tmp$ph1_release_effect <- 0
params_tmp$ph2_release_effect <- 0
params_tmp$ph3_release_effect <- 0

params_tmp$testing_effect <- 0

int_num = length(interventions$distancing_list)
int_effects_tmp = as.list(rep(0, int_num))

mobfun_tmp = approxfun(1:mytmax, rep(1,mytmax) ,method='linear', rule=2)
interventions$mobility = mobfun_tmp

params_tmp$H_lag <- 0
params_tmp$D_lag <- 0


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
                    tmax=mytmax, 
                    interventions=interventions, 
                    int_effects = int_effects_tmp,
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
rposterior = function(params, posterior, interventions){

par_smpl = as.list(posterior[sample(c(1:nrow(posterior)), size = 1), ])

var_names = c('q_Is', 
             'beta_pre',
             'k_A',
             'gamma_H',
             'm_H',
             'H_lag',
             'D_lag',
             'school_closure_effect',
             'lockdown_effect',
             'ph1_release_effect',
             'testing_effect',
             'E_init',
             'time_num', 
             'ramp_time',
             'm_Hbar_mult')

  for (k in var_names){
    params[[k]] = par_smpl[[k]]
  }

  params$q_Im = 1 - params$q_A - params$q_Is
  params$gamma_Hbar = params$gamma_H
  
  E_init_state0 = params$E_init * E_INIT_COUNTY

  state0 = get_state0_params(params, E_init_state0, interventions, CT_POPULATIONS, CT_ADJ, COUNTY_CAPACITIES)

  return( list (params, state0) )
}












#######################
## Run the model 
# draw_rparams controls uncertainty simulation: 
# if draw_rparams = TRUE, parameters are drawn independently using rparams()
# if draw_rparams = TRUE, parameters are drawn from joint posterior supplied as 'posterior'

get_sir_results = function(daymax, 
                           int_off_dates,
                           nsim=1,
                           params,
                           state0,
                           posterior,
                           draw_rparams = FALSE,
                           seed = NULL, 
                           CI = 0.95) {

  dayseq = seq(day0, daymax, by="day")
  tmax = as.numeric(difftime(daymax, day0, units="days"))

  # get a list of intervetion ramping times from default params
  int_ramping_times = get_ramping_times(params)
  
   # list of contact intervention functions 
   distancingfun_list = get_distancing_fun_list(dayseq, INT_START_DATES, int_off_dates, int_ramping_times)

   # mobility
   mobilityfun = get_mobility_fun(dayseq, MOB)

   # testing
   testingfun = get_testing_fun(dayseq, TESTING)

   # combine in the list of interventions 
   interventions = list(distancing_list=distancingfun_list, mobility=mobilityfun, testing=testingfun) 

 pars <- list()
 state0s <- list()
 interventions_list <- list()  
 
# parameters, set seed if given
  if(!is.null(seed)) set.seed(seed)
  
  if (draw_rparams == TRUE){
       if(nsim == 1){
     pars[[1]] = params 
     state0s[[1]] = state0
     interventions_list[[1]] = interventions
  } else { 
    for(i in 1:nsim){ 
       pars[[i]] = rparams(params)
       state0s[[i]] = state0
       interventions_list[[i]] = interventions
    }
  } 
} else {
  if(nsim == 1){
     pars[[1]] = params
     state0s[[1]] = state0
     interventions_list[[1]] = interventions
  } else { 
    for(i in 1:nsim){ 
       rpost_out = rposterior(params, posterior, interventions)
       pars[[i]] = rpost_out[[1]]
       state0s[[i]] = rpost_out[[2]]
       
      # update interventions with ramping times from posterior
       ramping_upd = get_ramping_times(pars[[i]])
       distancingfun_list_upd = get_distancing_fun_list(dayseq, INT_START_DATES, int_off_dates, ramping_upd)
       interventions_list[[i]] = list(distancing_list=distancingfun_list_upd, mobility=mobilityfun, testing=testingfun)
       }
  } 
}   

  sir_results = lapply(1:nsim, function(i){
    res = run_sir_model(state0=state0s[[i]], 
                        params=pars[[i]],   
                        region_adj=CT_ADJ, 
                        populations=CT_POPULATIONS, 
                        tmax=tmax, 
                        interventions=interventions_list[[i]],
                        int_effects = get_intervention_effects(pars[[i]]),
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
			                       lower = quantile(value, (1-CI)/2, na.rm=TRUE),
			                       upper = quantile(value, 1-(1-CI)/2, na.rm=TRUE))
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



get_snapshot <- function(data, date, where = "Connecticut"){
  toprint <- unique(data$variable)
  toprint <- toprint[grep(where, toprint)] 
  time.toprint <- as.numeric(as.Date(date, "%m/%d/%Y") - day0)
  out <- NULL

  sir_result_internal <- filter(data, variable%in%toprint) %>% 
                         filter(time == time.toprint) %>%
                         ungroup(variable) %>%
                         mutate(Date = day0 + time) %>%
                         mutate(variable = gsub(paste0(".",where), "", variable)) %>%
                         select(variable, time, mean, lower, upper) 
  return(data.frame(sir_result_internal))
  return(out)
}






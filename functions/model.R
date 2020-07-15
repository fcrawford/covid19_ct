
# state0 is the initial condition 
# params is a list of parameter values 
# region_adj is an adjacency matrix whose row/col labels are 


run_sir_model = function(state0, params, region_adj, populations, tmax, interventions, int_effects, capacities, deathfun) {

  if(is.null(interventions$distancing_list) || is.null(interventions$mobility) || is.null(interventions$testing) ||
     !is.list(interventions$distancing_list) || !is.function(interventions$mobility) || !is.function(interventions$testing)) { 
    stop("must specify intervention functions for distancing (list of functions), mobility and testing")
  }

  nregions = nrow(region_adj) # number of counties 

  if(length(state0) != nregions*11) stop("length of state0 must be nregions*11")

  region_names = rownames(region_adj)

  # indices 
  S_idx    = 1:nregions
  E_idx    = (  nregions+1):(2*nregions)
  I_s_idx  = (2*nregions+1):(3*nregions)
  I_m_idx  = (3*nregions+1):(4*nregions)
  A_idx    = (4*nregions+1):(5*nregions)
  H_idx    = (5*nregions+1):(6*nregions)
  Hbar_idx = (6*nregions+1):(7*nregions)
  NH_idx   = (7*nregions+1):(8*nregions)
  NI_idx   = (8*nregions+1):(9*nregions)
  D_idx    = (9*nregions+1):(10*nregions)
  R_idx    = (10*nregions+1):(11*nregions)

  # number of adjacent regions for each region: use to keep overall beta the same for each county
  region_adj_num <- rowSums(region_adj)

  beta_pre = params$beta_pre

  # here we define the intervention function: 
  # note that this function should be <= 1
  # and the effect of different interventions is additive, so coefficients must sum to <= 1.  
  # check that distancing functions list is the same length as intervention effects list
  if (length(interventions$distancing_list)!=length(int_effects)) {stop("list of effects has to be the same length as list of distancing functions")}
  
  # contact intervention function
  contact_intervention_fun = function(t_tmp) {
    
    distance = 0
    for (k in 1:length(interventions$distancing_list)){
      distance = distance + interventions$distancing_list[[k]](t_tmp) * int_effects[[k]]
    }
      return((1 - distance) * interventions$mobility(t_tmp))
    }

  # quick integrity check:
  if(any(contact_intervention_fun(1:(tmax+1))<0)) {
    plot(contact_intervention_fun(1:(tmax+1)))
    stop("contact_intervention_fun returns negative values. Check that intervention prameters sum to <= 1")
  }

  # region-wise beta transmission matrix 
  beta_matrix  = ( (1-params$k_n)*diag(1,nregions) + params$k_n*(1/region_adj_num)*region_adj ) * beta_pre 

  # CFR
  #params$m_Hbar = params$m_H * params$m_Hbar_mult
  #params$m_NH   = params$m_H * params$m_NH_mult

  # testing effects
  params$testing_effect_Im = params$testing_effect
  params$testing_effect_A = params$testing_effect * params$te_A_mult

  model <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      S = state[S_idx]
      E = state[E_idx]
      I_s = state[I_s_idx]
      I_m = state[I_m_idx]
      A = state[A_idx]
      H = state[H_idx]
      Hbar = state[Hbar_idx]
      NH = state[NH_idx]
      NI = state[NI_idx]
      D = state[D_idx]
      R = state[R_idx]

      # now the overall "effective" beta matrix is a product of the pre-intervention beta matrix 
      # and the intervention function evaluated at the current time
      beta = beta_matrix * contact_intervention_fun(time)

      # fatality ratios
      m_H_cur = m_H * deathfun(time)
      m_Hbar = m_H_cur * m_Hbar_mult
      m_NH = m_H_cur * m_NH_mult
  
      # effect of testing on early isolation / recovery for mild  and asymptomatic 
      a_t_Im = interventions$testing(time) * testing_effect_Im
      a_t_A  = interventions$testing(time) * testing_effect_A

      # update gamma_NI to account for earlier isolation
      red_im_dur = 1/alpha_Im - 1/( (1 + a_t_Im) * alpha_Im ) # reduction in duration of infectiousness before isolation
      gamma_NI = 1 / (1/gamma_NI + red_im_dur) # new gamma_NI accounts for increased duration of isolation, remains unchanged if a_t_Im = 0
      
      # the is no isolation compartment for asymptomatic: testing effect among asymptomatic moves them sooner directly to R
      
      # Hospital capacities
      actual_capacities = as.numeric(lapply(COUNTY_CAPACITIES, function(cap)cap(time)))

      smoothstep = function(x) 1/(1+exp(-0.5*x))
      Hospital_capacities_breached = smoothstep(H-actual_capacities)
      
      ##################
      
      # delta      1 / latency period
      
      # q_A        % infectious and asymptomatic
      # q_Is       % infectious and severe
      # q_Im       % infectious and mild
      # q_H        % among severe who will transition to either H or Hbar (severe outside or nursing homes, prisons, etc.)
      
      # alpha_A    recovery/removal rate of asymptomatic (1/ duration of infectiousness of asymptomatic)
      # alpha_Is   1 / time between end of latency and beginning of potential hospitalization (duration of infectiousness of severe cases) 
      # alpha_Im   rate at which mild transition to isolation (NI) (1/ duration of infectiousness of mild)
      
      # gamma_H    recovery/removal rate of hospitalized (1 / length of stay in hospital)
      # gamma_Hbar recovery/removal rate of hospital overflow group (1 / time from potential hospitalization to recovery or death among Hbar)
      # gamma_NH   1 / time between transitioning out of infectious state until recovery or death among severe in nursing homes and alike
      # gamma_NI   1 / duration of isolation of mild cases until recovery

      # m_H        hospital CFR
      # m_Hbar     hospital overflow CFR
      # m_NH       nursing homes CFR
      
      # k_A        relative transmission from asymptomatic infectives (compared to symptomatic)
      # k_Is       isolation coefficient among severe before they get admitted to a hospital: added to delay death without increasing FOI

      
      #################
  
      
      dS    <-            -S*( beta %*% ( k_Is * I_s + I_m + k_A * A )/populations )    # populations normalized contact rates
      
      dE    <-  -delta*E + S*( beta %*% ( k_Is * I_s + I_m + k_A * A )/populations ) 
      
      dI_s  <-  q_Is * delta*E - alpha_Is * I_s 
      
      dI_m  <-  (1 - q_A - q_Is) * delta*E - (1 + a_t_Im) * alpha_Im * I_m
      
      dA    <-  q_A * delta*E - (1 + a_t_A) * alpha_A * A
      
      dH    <-  q_H * alpha_Is * I_s * (1-Hospital_capacities_breached) - gamma_H * H
      
      dHbar <-  q_H * alpha_Is * I_s * Hospital_capacities_breached - gamma_Hbar * Hbar
      
      dNH   <-  (1 - q_H) * alpha_Is * I_s - gamma_NH * NH
      
      dNI   <-  (1 + a_t_Im) * alpha_Im * I_m - gamma_NI * NI 
      
      dD    <-  gamma_H * m_H_cur * H     + gamma_Hbar * m_Hbar * Hbar     + gamma_NH * m_NH * NH
      
      dR    <-  gamma_H * (1-m_H_cur) * H + gamma_Hbar * (1-m_Hbar) * Hbar + gamma_NH * (1-m_NH) * NH + gamma_NI * NI + (1 + a_t_A) * alpha_A * A

      return(list(c(dS,dE,dI_s,dI_m,dA,dH,dHbar,dNH,dNI,dD,dR)))
    })
  }


  ts = seq(0,tmax,by=1)

  out = as.data.frame(ode(y=state0, times=ts, func=model, parms=params, method="lsodes"))

  names(out) = c("time",
                 paste0("S.", region_names, sep=""),
                 paste0("E.", region_names, sep=""),
                 paste0("I_s.", region_names, sep=""),
                 paste0("I_m.", region_names, sep=""),
                 paste0("A.", region_names, sep=""),
                 paste0("H.", region_names, sep=""),
                 paste0("Hbar.", region_names, sep=""),
                 paste0("NH.", region_names, sep=""),
                 paste0("NI.", region_names, sep=""),
                 paste0("D.", region_names, sep=""),
                 paste0("R.", region_names, sep=""))

  out$S.Connecticut    = rowSums(out[,paste0("S.", region_names, sep="")])
  out$E.Connecticut    = rowSums(out[,paste0("E.", region_names, sep="")])
  out$I_s.Connecticut  = rowSums(out[,paste0("I_s.", region_names, sep="")])
  out$I_m.Connecticut  = rowSums(out[,paste0("I_m.", region_names, sep="")])
  out$A.Connecticut    = rowSums(out[,paste0("A.", region_names, sep="")])
  out$H.Connecticut    = rowSums(out[,paste0("H.", region_names, sep="")])
  out$Hbar.Connecticut = rowSums(out[,paste0("Hbar.", region_names, sep="")])
  out$NH.Connecticut    = rowSums(out[,paste0("NH.", region_names, sep="")])
  out$NI.Connecticut    = rowSums(out[,paste0("NI.", region_names, sep="")])
  out$D.Connecticut    = rowSums(out[,paste0("D.", region_names, sep="")])
  out$R.Connecticut    = rowSums(out[,paste0("R.", region_names, sep="")])

  
  ## track modified or combined compartments for CT
  
  # daily new infections 
  out$dailyI.Connecticut <- params$delta * out$E.Connecticut
  # daily new hospitalizations
  out$dailyH.Connecticut <- params$q_H * params$alpha_Is * out$I_s.Connecticut
  # daily new severe cases in nursing homes and ALFs
  out$dailyIs_NH.Connecticut <- (1 - params$q_H) * params$q_Is * params$delta * out$E.Connecticut
  
  # cumulative infections, including those who died
  out$cum_modI.Connecticut <- cumsum(out$dailyI.Connecticut)
  # cumulative hospitalizations
  out$cum_modH.Connecticut <- cumsum(out$dailyH.Connecticut)
  # cumulative severe cases in nursing homes and ALFs
  out$cum_Is_NH.Connecticut <- cumsum(out$dailyIs_NH.Connecticut)
  
  pop_ct = sum(populations)
  # number and proportion of alive cumulative incidence: sum of currently infected and recovered
  out$alive_cum_incid_num.Connecticut = pop_ct - out$S.Connecticut - out$E.Connecticut - out$D.Connecticut
  out$alive_cum_incid_prop.Connecticut = out$alive_cum_incid_num.Connecticut / (pop_ct - out$D.Connecticut)

  # current number and prevalence of infections outside of hospitals and nursing homes at the state level; includes asymptomatic and those isolated at home
  out$currentI.Connecticut = out$A.Connecticut + out$I_m.Connecticut + out$NI.Connecticut + params$q_H * out$I_s.Connecticut + out$Hbar.Connecticut
  out$currentI_prop.Connecticut = out$currentI.Connecticut / (pop_ct - out$D.Connecticut)
  
  
  # daily and cumulative hospitalizations for counties
  for(i in region_names){
    dailyH <- params$q_H * params$alpha_Is * out[[paste0("I_s.", i)]]
    out[[paste0("cum_modH.", i)]] <- cumsum(dailyH)
  }
  
  
  
  # add variables to plot lagged deaths and hospitalizations 
  lag_shift <- function(v,lg) { c(rep(0,lg), v[1:(length(v) - lg)]) }
  
  D_lag <- params$D_lag
  H_lag <- params$H_lag
  detection_lag <- params$detect_lag
  
  # time-shift
  out$time_D_lag <- out$time + D_lag
  out$time_H_lag <- out$time + H_lag

  # lagged compartments in Connecticut: deaths, current hospitalizations, cumulative hospitalizations
  out$rD.Connecticut    = lag_shift(out$D.Connecticu, D_lag)
  out$rH.Connecticut    = lag_shift(out$H.Connecticut, H_lag)
  out$rHbar.Connecticut = lag_shift(out$Hbar.Connecticut, H_lag)
  out$rcum_modH.Connecticut = lag_shift(out$cum_modH.Connecticut, H_lag)
  out$detectE.Connecticut = lag_shift(out$E.Connecticut, detection_lag)
  
  # lagged compartments for counties: deaths and current hospitalizations
  for(i in region_names){
   out[[paste0("rD.", i)]]    = lag_shift(out[[paste0("D.", i)]], D_lag) 
   out[[paste0("rH.", i)]]    = lag_shift(out[[paste0("H.", i)]], H_lag)
   out[[paste0("rHbar.", i)]] = lag_shift(out[[paste0("Hbar.", i)]], H_lag)
  }
  


  ## intervention patterns ##  
  out$intervention_pattern = contact_intervention_fun(1:(tmax+1))
  #out$intervention_schools = interventions$schools(1:(tmax+1))
  #out$intervention_lockdown = interventions$lockdown(1:(tmax+1))
  out$mobility = interventions$mobility(1:(tmax+1))
  out$intervention_testing = interventions$testing(1:(tmax+1))
  
  # dynamics of hospital CFR
  out$hCFR = params$m_H * deathfun(1:(tmax+1))
  
  # dynamics of rates of removal from infectious compartments I_m and A
  out$alpha_Im.t = params$alpha_Im * (1 + params$testing_effect_Im * out$intervention_testing)
  out$alpha_A.t = params$alpha_A * (1 + params$testing_effect_A * out$intervention_testing)
    
  # dynamic of aggregate FOI
  out$FOI = ( params$q_A * params$k_A / out$alpha_A.t 
              + (1 - params$q_Is - params$q_A) / out$alpha_Im.t 
              + params$q_Is * params$k_Is / params$alpha_Is  )
  
  # compute dynamics of R_eff at the state level
  # R_eff takes into account social distancing, testing effect in reducing FOI and depletion of susceptibles

  out$R_eff.Connecticut = params$beta_pre * out$FOI * out$intervention_pattern * out$S.Connecticut / (pop_ct - out$D.Connecticut)
  
  return(out)
}



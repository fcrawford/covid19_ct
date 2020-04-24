


# intervention functions
# these return functions, instead of values, because they need to be evaluated at non-integer points by the ODE solver


############################

get_state_lockdown_fun = function(dayseq, offdate) {
  if(offdate<state_lockdown_start) {stop("off date cannot be before lockdown start")}

  tmax = as.numeric(max(dayseq)-day0)

  lockdown = sapply(dayseq,function(dy) {
    if(dy<state_lockdown_order) {0}
    else if(dy<state_lockdown_start) {as.numeric(dy-state_lockdown_order)/as.numeric(state_lockdown_start-state_lockdown_order)}
    else if(dy<offdate) {1}
    else {0}
  })
  return(approxfun(1:(tmax+1), lockdown, method="linear", rule=2))
}


##########################


# note this is 1 when school is in session, and NOT closed

get_school_in_session_fun = function(dayseq, schools_reopen_date) {
  tmax = as.numeric(max(dayseq)-day0)
  school_in_session = sapply(dayseq,function(dy) {
    if(dy<state_schools_close) {1}
    else if(dy<schools_reopen_date) {0}
    else {1}
  })
  return(approxfun(1:(tmax+1), school_in_session, method="linear", rule=2))
}

###########################

# testing

get_testing_on_fun = function(dayseq, testing_on_date) {
  tmax = as.numeric(max(dayseq)-day0)
  testing_on = sapply(dayseq,function(dy) {
    if(dy<testing_on_date) {0}
    else {1}
  })
  return(approxfun(1:(tmax+1), testing_on, method="linear", rule=2))
}

################################
# distancing After lockdown

get_distancing_on_fun = function(dayseq, distancing_on_date) {
  tmax = as.numeric(max(dayseq)-day0)
  distancing_on = sapply(dayseq,function(dy) {
    if(dy<distancing_on_date) {0}
    else {1}
  })
  return(approxfun(1:(tmax+1), distancing_on, method="linear", rule=2))
}





###########################

get_R0 = function(params) {
  R0 = params$beta_pre * (  params$q_A * params$k_A / params$gamma_A 
                          + params$q_Im / params$gamma_Im 
                          + (1 - params$q_Im - params$q_A) * params$q_ins * params$k_Is_ins / params$alpha  
                          + (1 - params$q_Im - params$q_A) * (1 - params$q_ins) * params$k_Is_noins / params$gamma_Is ) 
  return(R0)
}



get_R_lockdown = function(params){
  
  R0 = get_R0(params)
  R_lockdown = R0 * (1 - params$school_closure_effect - params$lockdown_effect)
  
  return(R_lockdown)
}



######################### 


release_age_thresh <- function(threshold_age, lockdown_effect){

  dat <- read.csv('../data/severity.csv', stringsAsFactors=FALSE) 

  # current age proportions
  ageprops = dat$pr_age

  # avg_severity_prop = ageprops %*% dat$pr_severe # not used

  release_ageprops = ageprops * (dat$max_age<threshold_age) # zero out older proportions 
  prop_released = sum(release_ageprops) # the proportion being released

  # normalized age proportions to be released
  release_ageprops_normalized = release_ageprops / prop_released

  # new age proportion not locked down: a mixture of locked down ageprops and the released people in the rest, normalized
  new_ageprops = (lockdown_effect*ageprops + (1-lockdown_effect)*prop_released*release_ageprops_normalized) / (lockdown_effect + (1-lockdown_effect)*prop_released)

  # our lockdown effect gets bigger (worse)
  new_lockdown_effect = lockdown_effect + (1-lockdown_effect)*prop_released

  # new average severity for released
  new_avg_severity_prop = sum(new_ageprops * dat$pr_severe)

  return(list(new_lockdown_effect=new_lockdown_effect, new_avg_severity_prop=new_avg_severity_prop))

}


###########################
# testing

#par(mfrow=c(3,1), mar=c(3,4,3,0),bty="n")

#lockfun = get_state_lockdown_fun(dmy("01/06/2020"))
#plot(dayseq, lockfun(1:(tmax+1)), type="l", main="State shelter-in-place order", col="orange")

#schoolsfun = get_school_in_session_fun()
#plot(dayseq, schoolsfun(1:(tmax+1)), type="l", main="State schools in session", col="purple")

# for example: 
#relative_transmission = 1 - (0.3*lockfun(1:(tmax+1)) + 0.4*(1-schoolsfun(1:(tmax+1))))
#plot(dayseq, relative_transmission, type="l", main="Relative transmission", col="blue", ylim=c(0,1))






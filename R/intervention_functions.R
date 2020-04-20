


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




###########################


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


  #prop_adj <- severity_by_age$pr_age * 0.3

  #for (i in 3: nrow(severity_by_age)){
    #if(severity_by_age$max_age[i] < threshold_age) {
      #prop_adj[i] <- severity_by_age$pr_age[i] * 0.8
    #}  
  #}
  #
  #prop_adj <- prop_adj / sum(prop_adj)
  #
  #q_Is_new <- sum(prop_adj * severity_by_age$pr_severe)
  #return(q_Is_new)
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






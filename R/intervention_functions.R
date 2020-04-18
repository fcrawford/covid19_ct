


# intervention functions
# these return functions, instead of values, because they need to be evaluated at non-integer points by the ODE solver


############################

state_lockdown_order = dmy("20/03/2020") # order date
state_lockdown_start = dmy("23/03/2020") # actual start date

get_state_lockdown_fun = function(offdate=dmy("01/01/2021"), post_off_effect=0) {
  if(offdate<state_lockdown_start) {stop("off date cannot be before lockdown start")}

  lockdown = sapply(dayseq,function(dy) {
    if(dy<state_lockdown_order) {0}
    else if(dy<state_lockdown_start) {as.numeric(dy-state_lockdown_order)/as.numeric(state_lockdown_start-state_lockdown_order)}
    else if(dy<offdate) {1}
    else {post_off_effect} 
  })
  return(approxfun(1:(tmax+1), lockdown, method="linear", rule=2))
}


##########################

state_schools_close = dmy("13/03/2020")

# note this is 1 when school is in session, and NOT closed

get_school_in_session_fun = function(state_schools_reopen=dmy("01/09/2020")) {
  school_in_session = sapply(dayseq,function(dy) {
    if(dy<state_schools_close) {1}
    else if(dy<state_schools_reopen) {0}
    else {1}
  })
  return(approxfun(1:(tmax+1), school_in_session, method="linear", rule=2))
}

###########################


# not used:

get_relative_transmission_fun = function(lockdown_off_date, state_schools_reopen_date) {

  lockfun = get_state_lockdown_fun(dmy("01/06/2020"))
  schoolsfun = get_school_in_session_fun(dmy("01/09/2020"))
  approxfun(1:(tmax+1), 1 - (0.3*lockfun(1:(tmax+1)) + 0.4*(1-schoolsfun(1:(tmax+1)))))

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






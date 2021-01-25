


# intervention functions
# these return functions, instead of values, because they need to be evaluated at non-integer points by the ODE solver



##########################

# generic function for switching contact interventions on and off

get_contact_intervention_fun = function(dayseq, startdate, offdate, ramping) {
  if(offdate < startdate) {stop("off date cannot be before start date")}

  tmax = as.numeric(max(dayseq)-day0)
  ramping_end_date = as.Date(startdate + ramping)

  if(ramping_end_date>offdate) {stop("ramping value too large")}

  cont_int = sapply(dayseq,function(dy) {
    if(dy<startdate) {0}
    else if(dy<ramping_end_date) {as.numeric(dy-startdate)/ramping}
    else if(dy<offdate) {1}
    else {0}
  })
  #return(approxfun(1:(tmax+1), cont_int, method="linear", rule=2))
  return(approxfun(0:tmax, cont_int, method="linear", rule=2))
}





##########################

# list of functions to calculate aggregate distancing function

# startdates and offdates are lists of the same length

# keep the list function, may be useful in the future
# in the model with random effects, use a list of length 1
# schools reopening intervention added on 08/17 for projections

get_distancing_fun_list = function(dayseq, startdates, offdates, rampings) {
  if(length(startdates)!=length(offdates)) {stop("list of start dates has to be the same length as list of off dates")}
  if(length(startdates)!=length(rampings)) {stop("list of start dates has to be the same length as list of ramping times")}

  distancing_fun_list = list()
  for (k in 1:length(startdates)){
   distancing_fun_list[[k]] = get_contact_intervention_fun(dayseq, startdates[[k]], offdates[[k]], rampings[[k]])
  }
  return(distancing_fun_list)
}






#################################

# random effects


get_random_effect_fun = function(dayseq, random_effect_at, random_effect_values){
  
  re1_at = random_effect_at[1]
  random_effect_at = c((re1_at-14), random_effect_at)
  random_effect_values = c(0, random_effect_values)
  
  ff = approxfun(random_effect_at, random_effect_values, method="linear", rule=2)
  
  xx = as.data.frame(cbind(c(0:(length(dayseq)-1)), ff(c(0:(length(dayseq)-1) )) ) )
  colnames(xx) = c('time', 're')
  
  sp <- smooth.spline(xx$time[(re1_at-1):nrow(xx)], xx$re[(re1_at-1):nrow(xx)], nknots=(length(random_effect_at)) )
  xx$smooth.re = c(rep(0,(re1_at-2)),sp$y)
  xx$smooth.re[1:re1_at]=0
  xx$smooth.re[(re1_at+1):(re1_at+3)]=xx$re[(re1_at+1):(re1_at+3)]
  
  tmax = as.numeric(max(dayseq)-day0)

  re.smooth = sapply(c(0:tmax),function(dy) {
    if(dy<0) {0}
    else if(dy<tmax) {xx$smooth.re[which(xx$time==dy)] }
    else {xx$smooth.re[nrow(xx)]}
  })

   return(approxfun(0:tmax, re.smooth, method="linear", rule=2))
}








##########################
 
# mobility (close contact)

 # returns 1 for dates prior to 03/01 (first day of mobility data) 
 # smoothed mobility for dates when data is available
 # the latest available mobilily value is extrapolated into the future 

get_mobility_fun = function(dayseq, mob.data){
  
   tmax = as.numeric(max(dayseq)-day0)
   
   mob_day0 = ymd(mob.data$date[1])
   mob_daymax = ymd(mob.data$date[nrow(mob.data)])
   
   mob = sapply(dayseq,function(dy) {
     if (dy <= mob_day0) {1}
     else if(dy < mob_daymax) {mob.data$smooth[which( ymd(mob.data$date) == dy)]} 
     else {mob.data$smooth[nrow(mob.data)] }
   })
   
   #time0 = as.numeric(day0 - mob_day0)
   #return(approxfun(time0:(tmax+time0), mob, method="linear", rule=2))
   return(approxfun(0:tmax, mob, method="linear", rule=2))
}






##########################
 
# testing / contact tracing

 # returns 0 for dates prior to March 1 (first reported date) 
 # smoothed log(daily number of tests) - log(daily tests on the first available date) for dates when data is available
 # the latest available testing value is extrapolated into the future 

get_testing_fun = function(dayseq, test.data){
  
   tmax = as.numeric(max(dayseq)-day0)
   
   test_day0 = ymd(test.data$date[1])
   test_daymax = ymd(test.data$date[nrow(test.data)])
   
   log.smooth.test_day0 = log(test.data$smooth[1])
     
   testing = sapply(dayseq,function(dy) {
     if (dy < test_day0) {0}
     else if(dy < test_daymax) { max( (log(test.data$smooth[which( ymd(test.data$date) == dy)]) - log.smooth.test_day0) , 0) } 
     else { log(test.data$smooth[nrow(test.data)] ) - log.smooth.test_day0 }
   })
   
   #time0 = as.numeric(day0 - test_day0)
   #return(approxfun(time0:(tmax+time0), testing, method="linear", rule=2))
   return(approxfun(0:tmax, testing, method="linear", rule=2))
}









#################################
 
# relative change in hospital CFR

# returns first value in the data at the beginning
# smoothed daily death hazard for dates when data is available
# the latest available relative hazard into the future


get_death_fun = function(dayseq, dhaz.data){
  
  tmax = as.numeric(max(dayseq)-day0)
  
  dhaz_day0 = ymd(dhaz.data$date[1])
  dhaz_daymax = ymd(dhaz.data$date[nrow(dhaz.data)])
  
  dhaz = sapply(dayseq, function(dy) {
     if (dy <= dhaz_day0) {dhaz.data$smooth.rel_haz[1]}
     else if(dy < dhaz_daymax) {dhaz.data$smooth.rel_haz[which( ymd(dhaz.data$date) == dy)]} 
     else {dhaz.data$smooth.rel_haz[nrow(dhaz.data)] }
   })
   
   #time0 = as.numeric(d_day0 - day0)
 
   return(approxfun(0:tmax, dhaz, method="linear", rule=2))
}










#################################
 
# relative change in SEVERITY: proportion of cases > 60 y.o.

# returns first value in the data at the beginning (starting 05/01)
# smoothed prportion of cases > 60 y.o. dates when data is available
# the latest available relative value into the future


get_severity_fun = function(dayseq, dsev.data){
  
  tmax = as.numeric(max(dayseq)-day0)
  
  dsev_day0 = ymd(dsev.data$date[1])
  dsev_daymax = ymd(dsev.data$date[nrow(dsev.data)])
  
  dsev = sapply(dayseq, function(dy) {
     if (dy <= dsev_day0) {dsev.data$sev.measure[1]}
     else if(dy < dsev_daymax) {dsev.data$sev.measure[which( ymd(dsev.data$date) == dy)] } 
     else {dsev.data$sev.measure[nrow(dsev.data)] }
   })
  
   #time0 = as.numeric(d_day0 - day0)
 
   return(approxfun(0:tmax, dsev, method="linear", rule=2))
}






#################################

# rate of hospital departure: 1/hospital LOS (length of stay)

# smoothed monthly averages for hospital LOS, days 
# latest value is projected into the future


get_hosp_discharge_fun = function(dayseq, hlos.data){

  tmax = as.numeric(max(dayseq)-day0)
  
  dat_day0 = ymd(hlos.data$date[1])
  dat_daymax = ymd(hlos.data$date[nrow(hlos.data)])
  
  hlos.data$rate = 1/hlos.data$smooth.alos
    
  dhdisch = sapply(dayseq, function(dy) {
     if (dy <= dat_day0) {hlos.data$rate[1]}
     else if(dy < dat_daymax) {hlos.data$rate[which( ymd(hlos.data$date) == dy)] } 
     else {hlos.data$rate[nrow(hlos.data)] }
   })
  
   #time0 = as.numeric(d_day0 - day0)
 
   return(approxfun(0:tmax, dhdisch, method="linear", rule=2))
}










































############################################
### not used in the latest model version ###
############################################

######################################

# # current hospitalizations coming from congregate settings
# 
# get_current_congregate_hosp_fun = function(dayseq, hosp_cong.data){
#   
#   tmax = as.numeric(max(dayseq)-day0)
#   
#   dat_day0 = ymd(hosp_cong.data$date[1])
#   dat_daymax = ymd(hosp_cong.data$date[nrow(hosp_cong.data)])
#   
#   d = sapply(dayseq, function(dy){
#     if (dy < dat_day0) {0}
#     else if (dy < dat_daymax) {hosp_cong.data$cur_hosp_cong[ which (ymd(hosp_cong.data$date) == dy)] }
#     else {0}
#   })
#   
#   return( approxfun(0:tmax, d, method="linear", rule=2))
# }
# 
# 
# 
# 
# 
# 
# 
# 
# ######################################
# 
# # arrival hospitalizations coming from congregate settings
# 
# get_arrival_congregate_hosp_fun = function(dayseq, hosp_cong.data){
#   
#   tmax = as.numeric(max(dayseq)-day0)
#   
#   dat_day0 = ymd(hosp_cong.data$date[1])
#   dat_daymax = ymd(hosp_cong.data$date[nrow(hosp_cong.data)])
#   
#   d = sapply(dayseq, function(dy){
#     if (dy < dat_day0) {0}
#     else if (dy < dat_daymax) {hosp_cong.data$arrival_hosp_cong[ which (ymd(hosp_cong.data$date) == dy)] }
#     else {0}
#   })
#   
#   return( approxfun(0:tmax, d, method="linear", rule=2))
# }
# 
# 






###########################

# linear interpolation between the dates at which random effects are evaluated
# prior to the first random effect date, function should return 0

#get_random_effect_fun = function(dayseq, random_effect_at, random_effect_values){
  
#  re1_at = random_effect_at[1]
#  random_effect_at = c((re1_at-14), random_effect_at)
#  random_effect_values = c(0, random_effect_values)
  
#  return(approxfun(random_effect_at, random_effect_values, method="linear", rule=2))
# }


# R0 based on beta and FOI parameters

#get_R0 = function(params) {
#  R0 = params$beta_pre * ( params$q_A * params$k_A / params$alpha_A 
#                          + (1 - params$q_Is - params$q_A) / params$alpha_Im 
#                          + params$q_Is * params$k_Is / params$alpha_Is  ) 
#  return(R0)
#}




#get_state_lockdown_fun = function(dayseq, offdate) {
#  if(offdate<state_lockdown_start) {stop("off date cannot be before lockdown start")}

#  tmax = as.numeric(max(dayseq)-day0)

#  lockdown = sapply(dayseq,function(dy) {
#    if(dy<state_lockdown_order) {0}
#    else if(dy<state_lockdown_start) {as.numeric(dy-state_lockdown_order)/as.numeric(state_lockdown_start-state_lockdown_order)}
#    else if(dy<offdate) {1}
#    else {0}
#  })
#  return(approxfun(1:(tmax+1), lockdown, method="linear", rule=2))
#}

# ##########################
# 
# # note this is 1 when school is in session, and NOT closed
# 
# get_school_in_session_fun = function(dayseq, schools_reopen_date) {
#   tmax = as.numeric(max(dayseq)-day0)
#   school_in_session = sapply(dayseq,function(dy) {
#     if(dy<state_schools_close) {1}
#     else if(dy<schools_reopen_date) {0}
#     else {1}
#   })
#   return(approxfun(1:(tmax+1), school_in_session, method="linear", rule=2))
# }
# 
# ################################
# # distancing After lockdown
# # stepdown dates govern return to normal contact
# # e.g. get_distancing_on_fun(dayseq, ymd("2020-06-01"), ymd(c("2020-07-01", "2020-07-15", "2020-07-30")))
# 
# get_distancing_stepdown_fun = function(dayseq, distancing_on_date, distancing_stepdown_dates) {
#   tmax = as.numeric(max(dayseq)-day0)
#   if(any(distancing_stepdown_dates < distancing_on_date)) stop("distancing stepdown dates must be greater than on date")
#   distancing_stepdown_dates = sort(distancing_stepdown_dates) 
#   distancing_stepdown_days = as.numeric(distancing_stepdown_dates - day0)
#   distancing_on_day = as.numeric(distancing_on_date - day0)
#   distancing_stepdown_fun = approxfun(c(0,distancing_on_day, distancing_stepdown_days), c(0,seq(1,0,length.out=length(distancing_stepdown_days)+1)), method="constant", rule=2)
# 
# 
#   return(distancing_stepdown_fun)
# }
# 

# ##########################
# get_mobility_fun_travel = function(dayseq, mob.data)  {
#   
#   tmax = as.numeric(max(dayseq)-day0)
#   
#   mob_day0 = ymd(mob.data$date[1])
#   mob_daymax = ymd(mob.data$date[nrow(mob.data)])
#   
#   lo = lowess(mob.data$mobility ~ mob.data$time, f=0.1)
#   mob.data$sm_mobility = lo$y
#   
#   mob = sapply(dayseq,function(dy) {
#     if (dy <= mob_day0) {1}
#     else if(dy < mob_daymax) {1 + mob.data$sm_mobility[which( ymd(mob.data$date) == dy)]} 
#     else {1 + mob.data$sm_mobility[nrow(mob.data)] }
#   })
#   
#   time0 = as.numeric(day0 - mob_day0)
# 
#   return(approxfun(time0:(tmax+time0), mob, method="linear", rule=2))
# }
# 
# 
# ##########################
# get_mobility_fun_stayput = function(dayseq, mob.data)  {
#   
#   tmax = as.numeric(max(dayseq)-day0)
#   
#   mob_day0 = ymd(mob.data$date[1])
#   mob_daymax = ymd(mob.data$date[nrow(mob.data)])
# 
#   # baseline proportion staying put: average of the first 7 days of mobility data
#   # transform such that baseline mobility measure is 1 and 0 corresponds to no mobility
#   bl_mobile_prop = 1 - mean(mob.data$stay_put[1:7])
#   mob.data$mobile_prop = 1 - mob.data$stay_put
#   mob.data$relative_mobility = 1 - (bl_mobile_prop - mob.data$mobile_prop)/bl_mobile_prop
#   
#   # smooth
#   lo = lowess(mob.data$relative_mobility ~ mob.data$time, f=0.1)
#   mob.data$sm_relative_mobility = lo$y
#   
#   mob = sapply(dayseq,function(dy) {
#     if (dy <= mob_day0) {1}
#     else if(dy < mob_daymax) {mob.data$sm_relative_mobility[which( ymd(mob.data$date) == dy)]} 
#     else {mob.data$sm_relative_mobility[nrow(mob.data)] }
#   })
#   
#   time0 = as.numeric(day0 - mob_day0)
# 
#   return(approxfun(time0:(tmax+time0), mob, method="linear", rule=2))
# }
# 
###########################

# testing

# get_testing_on_fun = function(dayseq, testing_on_date) {
#   tmax = as.numeric(max(dayseq)-day0)
#   testing_on = sapply(dayseq,function(dy) {
#     if(dy<testing_on_date) {0}
#     else {1}
#   })
#   return(approxfun(1:(tmax+1), testing_on, method="linear", rule=2))
# }





###########################
# test how functions look

#par(mfrow=c(3,1), mar=c(3,4,3,0),bty="n")

#lockfun = get_state_lockdown_fun(dmy("01/06/2020"))
#plot(dayseq, lockfun(1:(tmax+1)), type="l", main="State shelter-in-place order", col="orange")

#schoolsfun = get_school_in_session_fun()
#plot(dayseq, schoolsfun(1:(tmax+1)), type="l", main="State schools in session", col="purple")

# for example: 
#relative_transmission = 1 - (0.3*lockfun(1:(tmax+1)) + 0.4*(1-schoolsfun(1:(tmax+1))))
#plot(dayseq, relative_transmission, type="l", main="Relative transmission", col="blue", ylim=c(0,1))






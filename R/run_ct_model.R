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

source("model.R")
source("get_data.R")
source("intervention_functions.R")

########################
# Dates

# starting day
day0 = ymd("2020-03-01")

# ending day
daymax = ymd("2020-09-01")

# to use in the model
tmax = as.numeric(difftime(daymax, day0, units="days"))

# dayseq is the calendar day sequence we want to model
dayseq = seq(day0, daymax, by="day")

# match dayseq to time in the model and data
date.time <- as.data.frame(cbind(format(as.Date(dayseq)),c(0:(length(dayseq)-1))))
colnames(date.time) <- c("date", 'time')

# the month sequences are for display only
monthseq = seq(day0, daymax, by="month")
monthseq_lab = format(monthseq, "%b %Y")
daymonthseq = difftime(monthseq, day0, units="days")


########################

##### get reported data from CT #####
data <- get_ct_data(day0=day0)
dat_ct_state <- data$dat_ct_state
dat_ct_county <- data$dat_ct_county


############ modeling #############

# load fixed parameters
params_init = yaml.load_file("params.yml") 

# load region adjacency matrix 
adj = read.csv("../map/CT_adj_matrix.csv", stringsAsFactors=FALSE)
rownames(adj) = adj$X
adj = adj[,-1]
adj = as.matrix(adj)

# populations and initial conditions
nregions = nrow(adj)

init <- read.csv('../data/ct_init.csv', stringsAsFactors=FALSE) 
region_names = init$county

#########################
# set up initial state0

E_init = init$E
I_s_init = init$Is
I_m_init = init$Im
A_init = init$A
H_init = init$H
Hbar_init = rep(0,nregions)
D_init = init$D
R_init = init$R
S_init = init$population - (E_init + I_s_init + I_m_init + A_init + H_init + Hbar_init + D_init + R_init)

state0 = c(S=S_init, E=E_init, I_s=I_s_init, I_m=I_m_init, A=A_init, H=H_init, Hbar=Hbar_init, D=D_init, R=R_init)

########################
# draw random params

rparams = function() {
  params_tmp = params_init
  # sample new param values, and initial conditions!
  params_tmp$beta_pre = rnorm(1,mean=params_init$beta_pre,sd=0.1)
  # add more here
  return(params_tmp)
}


#######################
# set intervention patterns

lockfun = get_state_lockdown_fun(offdate=dmy("01/06/2020"), post_off_effect=params_init$distancing_effect)

schoolsfun = get_school_in_session_fun(state_schools_reopen=dmy("01/09/2020"))


#######################
# run the sims

nsim = 1

sir_results = lapply(1:nsim, function(i){
  res = run_sir_model(state0=state0, 
                      params=rparams(),  # note: effect_intvx is in params, and is not passed to run_sir_model separately 
                      region_adj=adj, 
                      populations=init$population, 
                      tmax=tmax, 
                      interventions=list(lockdown=lockfun, schools=schoolsfun)) # intervention functions! 
  res$sim_id = i
  res
})

#######################
# aggregate and summarize results across sims

sir_results_all = ldply(sir_results, rbind)
sir_results_long <- melt(sir_results_all, id.vars = c("time", "sim_id"))
sir_results_summary <- sir_results_long %>% group_by(variable, time) %>% 
			                   summarise(
                           mean = mean(value),
			                     lower = quantile(value, 0.05),
			                     upper = quantile(value, 0.95))

####################
# plotting 

plot_ct_region = function(region_name) {
  #compartment_plot_labels = c("D")
  #compartment_plot_names = c("Deaths")
  #compartment_plot_colors = rainbow(length(compartment_plot_labels))

  par(mar=c(3,4,3,0), bty="n")

  sir_result_region = filter(sir_results_summary, variable==paste("D.",region_name,sep=""))

  plot(0, type="n", xlab="", ylab="People", main=region_name, col="black", 
       ylim=c(0,max(sir_result_region$mean)), xlim=c(0,1.1*tmax), axes=FALSE)
  axis(1,at=daymonthseq, lab=monthseq_lab)
  axis(2)


  abline(v=Sys.Date()-day0, col="gray", lty=2)

  polygon(c(sir_result_region$time, rev(sir_result_region$time)), c(sir_result_region$lower, rev(sir_result_region$upper)), col=rgb(1,0,0,alpha=0.5), border=NA)

  lines(sir_result_region$time, sir_result_region$mean, col="red")
  # label mean
  text(sir_result_region$time[tmax+1], sir_result_region$mean[tmax+1], format(sir_result_region$mean[tmax+1],digits=1), pos=4, col="red")
  text(sir_result_region$time[tmax+1], sir_result_region$lower[tmax+1], format(sir_result_region$lower[tmax+1],digits=1), pos=4, col=rgb(1,0,0,alpha=0.5))
  text(sir_result_region$time[tmax+1], sir_result_region$upper[tmax+1], format(sir_result_region$upper[tmax+1],digits=1), pos=4, col=rgb(1,0,0,alpha=0.5))

  if(region_name == "Connecticut") {
    points(dat_ct_state$time, dat_ct_state$deaths, pch=16, col=rgb(1,0,0,alpha=0.5)) 
  } else {
    obs.region <- subset(dat_ct_county, county == region_name)
    obs.region$date <- ymd(obs.region$date)
    first.region.time <- round(as.numeric(difftime(obs.region$date[1], day0, units="days")),0)
    obs.region$time <- c(first.region.time:(nrow(obs.region)+first.region.time-1)) # add time variable that indexes time 
    #obs.region <- merge(obs.region, date.time, by='date')
    points(obs.region$time, obs.region$deaths, pch=16, col=rgb(1,0,0,alpha=0.5))
  }

	#legend(0, max(dat_ct_state$deaths), compartment_plot_names, lty=1, col=compartment_plot_colors, bty="n")

  # return some useful info
  # capacity exceeded?
  # describe intvx
	region_summary = paste("On ", format(daymax, "%b %d, %Y"),
											 " projections show ", format(sir_result_region$mean[tmax+1], digits=2),
											 " deaths in ", region_name,
                       sep="")

  return(region_summary)
}

#####################################
#
# @param which.plot the prefix of the compartment to plot, e.g., S, E, I_s, I_m. If a vector of more than one specified, it takes the sum of the compartment (only the mean!)
mapplot_ct_region = function(which.plot = "D", label = "Cumulative Deaths", ...) {
  map = CTmap
  ncol=3 # customize based on date range? 
  region_names <- as.character(map$NAME10)
  toplot <- paste(rep(which.plot,each=length(region_names)),
                  rep(region_names, length(which.plot)),sep=".")
  sir_result_internal = data.frame(filter(sir_results_summary, variable%in%toplot))
  sir_result_internal$County <- sir_result_internal$variable
  for(i in which.plot) sir_result_internal$County <- gsub(paste0(i,"\\."), "", sir_result_internal$County)
  
  # take the sum if necessary, remove lower and upper though
  if(length(which.plot) >  1){
    sir_result_internal <- aggregate(mean ~ County+time, data=sir_result_internal, FUN=function(x){mean(x)})
  }
  sir_result_internal$Date <- format(sir_result_internal$time + day0, "%B")
  sir_result_internal$Date <- factor(sir_result_internal$Date, unique(sir_result_internal$Date))
  # Plot last day cumulative at each month? Or take diff?
  t.index <- unique(sir_result_internal$time[format(sir_result_internal$time + day0, "%d")=="01"]) 
  t.index <- (t.index - 1)[-1]

  g <- mapPlot(data=subset(sir_result_internal, time %in% t.index), geo=map,
    by.data="County", by.geo = "NAME10",
    variables = "Date", values = "mean", is.long=TRUE, 
    legend.label = label, ...) #, direction=-1,  ...) 
  suppressMessages(
    g <- g + scale_fill_distiller(label,  palette = "Reds", direction=1) + ggtitle(paste("Projected", tolower(gsub("\n"," ",label)), "in Connecticut counties by month"))
  )
  return(g)
}

#obs.region <- subset(dat_ct_county, county == "New Haven")
########################

# make correspondence with calendar dates

# plot state overall, and by region  


par(mfrow=c(3,3))

plot_ct_region("Connecticut")
sapply(region_names, plot_ct_region)

plot.new()
par(mfrow=c(1,1))
plot(sir_results[[1]]$intervention_pattern, ylim=c(0,1), main="intervention pattern", type="l", bty="n")

# Death plot
plot.new()
death.map = mapplot_ct_region() #map=CTmap, ncol=3)
print(death.map)

# Cumulative hospitalization plot
plot.new()
hos.map = mapplot_ct_region("cum_modH",  "Cumulative Hospitalizations")  
print(hos.map)


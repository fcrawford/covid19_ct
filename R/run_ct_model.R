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
daymax = ymd("2020-06-01")

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



#############################3
# hypothetical end date of lockdown

lockdown_end_date=dmy("01/06/2020") 


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
  # sample new param values, 
  params_tmp$beta_pre = max(0,rnorm(1,mean=params_init$beta_pre,  sd=params_init$sd_beta_pre))
  params_tmp$delta = 1/max(0.001,rnorm(1,mean=1/params_init$delta,  sd=params_init$sd_recip_delta))
  params_tmp$gamma_H = 1/max(0.001,rnorm(1,mean=1/params_init$gamma_H,  sd=params_init$sd_recip_gamma_H))
  params_tmp$m_H = max(0,rnorm(1,mean=params_init$m_H,  sd=params_init$sd_m_H))
  params_tmp$m_H = min(1,params_tmp$m_H)
  params_tmp$m_Hbar_mult <- max(0,rnorm(1,mean=params_init$m_Hbar_mult,  sd=params_init$sd_m_Hbar_mult))
  #params_tmp$lockdown_effect <- min(1, rnorm(1,mean=params_init$lockdown_effect,  sd=params_init$sd_lockdown_effect))
  params_tmp$lockdown_effect <- rnorm(1,mean=params_init$lockdown_effect,  sd=params_init$sd_lockdown_effect)
  #params_tmp$lockdown_effect <- max(0, params_tmp$lockdown_effect)
  # add sampling of initial conditions
  return(params_tmp)
}


#######################
# set intervention patterns


lockfun = get_state_lockdown_fun(offdate=lockdown_end_date, post_off_effect=params_init$distancing_effect)

state_schools_reopen_date = dmy("01/09/2020")
schoolsfun = get_school_in_session_fun(state_schools_reopen=state_schools_reopen_date)


#######################
# run the sims

nsim = 100

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
			                     lower = quantile(value, 0.05, na.rm=TRUE),
			                     upper = quantile(value, 0.95, na.rm=TRUE))

####################
# plotting 
# @param which.plot the prefix of the compartment to plot, e.g., S, E, I_s, I_m. If a vector of more than one specified, it plots multiple lines

plot_ct_region = function(region_name, which.plot = "D", add=NULL) {
  #compartment_plot_labels = c("D")
  #compartment_plot_names = c("Deaths")
  #compartment_plot_colors = rainbow(length(compartment_plot_labels))
  lab.table <- data.frame(compartment=c("D","rD",
                                        "H","rH",
                                        "Hbar","cum_modH","S","E","I_s","I_m","A"),
                          color=c('#e41a1c','#e41a1c','#377eb8','#377eb8',
                                  '#4daf4a','#984ea3','#ff7f00','#ffff33',
                                  '#a65628','#f781bf','#999999'),
                          labels=c("Cumulative Deaths","Cumulative Deaths",
                                    "Hospitalizations","Hospitalizations",
                                    "Hospital Overflow",
                                    "Cumulative Hospitalizations",
                                    "Susceptible Population","Exposed Population",
                                    "Severe Infections","Mild Infections",
                                    "Asymptomatic Infections"))

  if("rH" %in% which.plot) add <- dat_ct_capacity

  par(mar=c(3,4,3,0), bty="n")
  toplot <- paste(rep(which.plot,each=length(region_name)),
                  rep(region_name, length(which.plot)),sep=".")
  sir_result_region= filter(sir_results_summary, variable%in%toplot)

  title <- paste0(lab.table$labels[lab.table$compartment==which.plot[1]], " in ", region_name)
  ymax <- max(sir_result_region$mean[sir_result_region$time <= tmax], na.rm=TRUE)
  
  if(!is.null(add)){
    if(region_name != "Connecticut"){
        sub.add <- subset(add, County == region_name)
      }else{
        sub.add <- aggregate(value~Date, data=add, FUN=function(x){sum(x)})
      }
      ymax <- max(ymax, sub.add$value)
  }

  plot(0, type="n", xlab="", ylab="People", main=title, col="black", 
       ylim=c(0,1.05*ymax), xlim=c(0,1.1*tmax), axes=FALSE)
  axis(1,at=daymonthseq, lab=monthseq_lab)
  axis(2)

  abline(v=Sys.Date()-day0, col="gray", lty=2)


  lab.table$color <- as.character(lab.table$color)
  for(i in 1:length(which.plot)){
    col.line <- lab.table$color[which(lab.table$compartment==which.plot[i])]
    col.polygon <- adjustcolor(col.line, alpha.f = 0.5)
    sir_result_region_sub <- filter(sir_result_region, variable==paste0(which.plot[i],".",region_name))
    polygon(c(sir_result_region_sub$time, rev(sir_result_region_sub$time)), c(sir_result_region_sub$lower, rev(sir_result_region_sub$upper)), col=col.polygon, border=NA)

    lines(sir_result_region_sub$time, sir_result_region_sub$mean, col=col.line)
    text(sir_result_region_sub$time[tmax+1], sir_result_region_sub$mean[tmax+1], format(sir_result_region_sub$mean[tmax+1],digits=2, big.mark=","), pos=4, col=col.line)
    text(sir_result_region_sub$time[tmax+1], sir_result_region_sub$lower[tmax+1], format(sir_result_region_sub$lower[tmax+1],digits=2, big.mark=","), pos=4, col=col.polygon)
    text(sir_result_region_sub$time[tmax+1], sir_result_region_sub$upper[tmax+1], format(sir_result_region_sub$upper[tmax+1],digits=2, big.mark=","), pos=4, col=col.polygon)
    if(!is.null(add)){
      lines(as.numeric(as.Date(sub.add$Date) - day0), sub.add$value, col='gray30', lty  = 2,  lwd=1.2)
      tmp1 = as.numeric(today() - day0)
      tmp2 = as.numeric(daymax - day0)
      capvalue = sub.add$value[length(sub.add$value)]
      lines(c(tmp1,tmp2), rep(capvalue,2), col='gray30', lty  = 2,  lwd=1.2)
    }
  }

  # Add observed deaths
  col.line <- lab.table$color[which(lab.table$compartment=="D")]
  if(region_name == "Connecticut" && "rD" %in% which.plot) {
    points(dat_ct_state$time, dat_ct_state$deaths, pch=16, cex=0.6, col=col.line) 
  } else if("rD" %in% which.plot) {
    obs.region <- subset(dat_ct_county, county == region_name)
    obs.region$date <- ymd(obs.region$date)
    first.region.time <- round(as.numeric(difftime(obs.region$date[1], day0, units="days")),0)
    obs.region$time <- c(first.region.time:(nrow(obs.region)+first.region.time-1)) # add time variable that indexes time 
    #obs.region <- merge(obs.region, date.time, by='date')
    points(obs.region$time, obs.region$deaths, pch=16, cex=0.6, col=col.line)
  }

  # Add observed hospitalization
  col.line <- lab.table$color[which(lab.table$compartment=="H")]
  if(region_name == "Connecticut" && "rH" %in% which.plot) {
      points(dat_ct_state$time, dat_ct_state$cur_hosp, pch=16, cex=0.6, col=col.line) 
  } else if("rH" %in% which.plot) {
    obs.region <- subset(dat_ct_county, county == region_name)
    obs.region$date <- ymd(obs.region$date)
    first.region.time <- round(as.numeric(difftime(obs.region$date[1], day0, units="days")),0)
    obs.region$time <- c(first.region.time:(nrow(obs.region)+first.region.time-1)) # add time variable that indexes time 
    #obs.region <- merge(obs.region, date.time, by='date')
    points(obs.region$time, obs.region$cur_hosp, pch=16, cex=0.6, col=col.line)
  }

	#legend(0, max(dat_ct_state$deaths), compartment_plot_names, lty=1, col=compartment_plot_colors, bty="n")

  # return some useful info
  # capacity exceeded?
  # describe intvx
  region_summary <- NULL
  if("D" %in% which.plot || "rD" %in% which.plot){
      sir_result_region_sub <- filter(sir_result_region, variable==paste0(which.plot[1],".",region_name))
      count <- sir_result_region_sub$mean[sir_result_region_sub$time==tmax]
      count.min <- sir_result_region_sub$lower[sir_result_region_sub$time==tmax]
      count.max <- sir_result_region_sub$upper[sir_result_region_sub$time==tmax]
      region_summary = paste(region_summary, "On ", format(daymax, "%B %d"),
                       " projections show ", format(count, digits=2, big.mark=","),
                       " cumulative deaths reported in ", region_name,
                       " with 90% uncertainty interval between ", format(count.min, digits=2, big.mark=","),
                       " and ", format(count.max, digits=2, big.mark=","),
                       ".",
                       sep="")    
  }
  if("H" %in% which.plot || "rH" %in% which.plot){
      sir_result_region_sub <- filter(sir_result_region, variable==paste0(which.plot[1],".",region_name))
      count <- max(sir_result_region_sub$mean, na.rm=TRUE)
      peak <- which.max(sir_result_region_sub$mean)
      count.min <- sir_result_region_sub$lower[peak]
      count.max <- sir_result_region_sub$upper[peak]
      region_summary = paste(region_summary, "On ", format(daymax, "%B %d"),
                       " projections show a peak of ", format(count, digits=2, big.mark=","),
                       " hospitalizations reported in ", region_name,
                       " with 90% uncertainty interval between ", format(count.min, digits=2, big.mark=","), 
                       " and ", format(count.max, digits=2, big.mark=","),
                       ". The dashed line shows historical and projected hospital capacity. ",
                       sep="")    
  }

  return(region_summary)
}

#####################################
#
# @param which.plot the prefix of the compartment to plot, e.g., S, E, I_s, I_m. If a vector of more than one specified, it takes the sum of the compartment (only the mean!)
mapplot_ct_region = function(which.plot = "D", label = "Cumulative Deaths", palette="Reds", ...) {
  map = CTmap
  #ncol=4 # customize based on date range? 
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
    g <- g + scale_fill_distiller(label,  palette=palette, direction=1) + ggtitle(paste("Projected", tolower(gsub("\n"," ",label)), "in Connecticut counties by month"))
  )
  return(g)
}

####################################

plot_interventions = function() {

  par(mar=c(3,3,3,0), bty="n")
  layout(matrix(c(1,2,3),nrow=3), heights=c(2,2,3))

  plot(sir_results[[1]]$intervention_schools, ylim=c(0,1), type="n", ylab="", xlab="", main="Schools in session", axes=FALSE)
  axis(1, at=daymonthseq, lab=monthseq_lab)
  axis(2)
  polygon(c(1,1:(tmax+1), tmax+1), c(0,sir_results[[1]]$intervention_schools, 0), col="orange", border=NA)
  abline(v=Sys.Date()-day0, col="gray", lty=2)

  plot(sir_results[[1]]$intervention_lockdown, ylim=c(0,1), type="n", ylab="", xlab="", main="Stay-at-home order in place", axes=FALSE)
  axis(1, at=daymonthseq, lab=monthseq_lab)
  axis(2)
  polygon(c(1,1:(tmax+1), tmax+1), c(0,sir_results[[1]]$intervention_lockdown, 0), col="orange", border=NA)
  abline(v=Sys.Date()-day0, col="gray", lty=2)


  plot(sir_results[[1]]$intervention_pattern, ylim=c(0,1),  type="n", ylab="", xlab="", main="Relative reduction in transmission", axes=FALSE)
  axis(1, at=daymonthseq, lab=monthseq_lab)
  axis(2)
  polygon(c(1,1:(tmax+1), tmax+1), c(0,sir_results[[1]]$intervention_pattern, 0), col="orange", border=NA)
  abline(v=Sys.Date()-day0, col="gray", lty=2)

}

####################################
# get_Cumulative("2020-07-01", "rD.Connecticut")
get_Cumulative <- function(date, tosum){
  sir_result_internal = data.frame(filter(sir_results_summary, variable%in%tosum))
  t = as.numeric(difftime(as.Date(date), day0, unit='days'))
  sir_result_internal = subset(sir_result_internal, time%in%t)
  sir_result_internal$mean  <- as.character(round(sir_result_internal$mean, 0))
  if(length(date)==1) return(sir_result_internal$mean)
  if(length(date)>1) return(paste(sir_result_internal$mean, collapse=", "))
}


########################

# make correspondence with calendar dates

# plot state overall, and by region  


#par(mfrow=c(3,3))
#plot_ct_region("Connecticut")
#sapply(region_names, plot_ct_region)

# test plot multiple lines
#plot_ct_region("Connecticut", c("H"), add = dat_ct_capacity)

# for(i in region_names){
#   summary <- plot_ct_region(i, c("H"), add = dat_ct_capacity)
#   print(summary)
# }




#plot.new()
#plot_interventions()

# Death plot
#plot.new()
#death.map = mapplot_ct_region() #map=CTmap, ncol=3)
#print(death.map)

# Cumulative hospitalization plot
#plot.new()
#hos.map = mapplot_ct_region("cum_modH",  "Cumulative\nHospitalizations")  
#print(hos.map)


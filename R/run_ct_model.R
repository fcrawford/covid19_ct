library(deSolve)
library(yaml)
library(lubridate)
library(plyr)
library(reshape2)
library(dplyr)
library(rgdal)
library(SUMMER)
library(ggplot2)
#####################

source("model.R")
source("get_data.R")

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
# run the sims

nsim = 30

sir_results = lapply(1:nsim, function(i){
  res = run_sir_model(state0=state0, params=rparams(), region_adj=adj, populations=init$population, tmax=tmax,
                      effect_intvx=0.5, intvx_time=25)
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

  par(mar=c(4,4,3,4), bty="n")

  sir_result_region = filter(sir_results_summary, variable==paste("D.",region_name,sep=""))

  plot(0, type="n", xlab="Time", ylab="People", main=region_name, col="black", 
       ylim=c(0,max(sir_result_region$mean)), xlim=c(0,1.1*tmax), axes=FALSE)
  axis(1,at=daymonthseq, lab=monthseq_lab)
  axis(2)


  abline(v=Sys.Date()-day0, col="gray", lty=2)

  polygon(c(sir_result_region$time, rev(sir_result_region$time)), c(sir_result_region$lower, rev(sir_result_region$upper)), col=rgb(1,0,0,alpha=0.5), border=NA)
  lines(sir_result_region$time, sir_result_region$mean, col="red")
  text(sir_result_region$time[tmax+1], sir_result_region$mean[tmax+1], format(sir_result_region$mean[tmax+1],digits=1), pos=4, col="red")

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

mapplot_ct_region = function(map, ...) {
  region_names <- as.character(map$NAME10)
  sir_result_internal = data.frame(filter(sir_results_summary, variable%in%paste("D.",region_names,sep="")))
  sir_result_internal$County <- gsub("D.", "", sir_result_internal$variable)
  sir_result_internal$Date <- format(sir_result_internal$time + day0, "%B")
  sir_result_internal$Date <- factor(sir_result_internal$Date, unique(sir_result_internal$Date))
  # Plot last day cumulative at each month? Or take diff?
  t.index <- unique(sir_result_internal$time[format(sir_result_internal$time + day0, "%d")=="01"]) 
  t.index <- (t.index - 1)[-1]

  g <- mapPlot(data=subset(sir_result_internal, time %in% t.index), geo=map,
    by.data="County", by.geo = "NAME10",
    variables = "Date", values = "mean", is.long=TRUE, 
    legend.label = "Deaths", direction=-1,  ...) 
  g <- g + scale_fill_distiller("Deaths",  palette = "Blues", direction=1)
  return(g)
}

#obs.region <- subset(dat_ct_county, county == "New Haven")
########################

# make correspondence with calendar dates

# plot state overall, and by region  


par(mfrow=c(3,3))

plot_ct_region("Connecticut")
sapply(region_names, plot_ct_region)

CTmap <- readOGR("../map/wgs84/countyct_37800_0000_2010_s100_census_1_shp_wgs84.shp")
mapplot_ct_region(map=CTmap, ncol=3)


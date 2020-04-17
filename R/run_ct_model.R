library(deSolve)
library(yaml)
library(lubridate)

#####################

source("model.R")

#####################

######### observations #########

# Olya: could you please clean this up, give informative variable names, and move the data loading/processing out of this file? 
# you can create a get_ct_data() function or something
# I want to keep the global namespace as clean as possible 

# Get COVID-19 Data from NYT
nyt <- read.csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
nyt <- subset(nyt, state == "Connecticut")
date <- as.character(nyt[dim(nyt)[1],  1])

## caclulate state-wide observations ##
nyt2 <- subset(nyt, select=c(date, cases, deaths))
nyt_ct <- aggregate(. ~ date, data=nyt2, FUN=function(x){sum(x)})
nyt_ct$time <- c(0:(nrow(nyt_ct)-1)) # add time variable that indexes time starting 0
time.date <- subset(nyt_ct, select=c(date,time)) # for future matching use
# calculate daily new reported cases and deaths 
nyt_ct$new_cases <- nyt_ct$new_deaths <- 0
nyt_ct$new_cases[1] <- nyt_ct$cases[1] 
for (i in 2:nrow(nyt_ct)){
   nyt_ct$new_cases[i] <- nyt_ct$cases[i] - nyt_ct$cases[i-1]
   nyt_ct$new_deaths[i] <- nyt_ct$deaths[i] - nyt_ct$deaths[i-1]
}
## add state-level hospitalizations (current counts): this csv file needs to be updated manually ##
ct.hosp <- read.csv('../data/ct_hosp.csv')
nyt_ct <- merge(nyt_ct, ct.hosp, by='time')



########################
# new date scheme 
# Olya could you please make sure that the data you use 

# starting day
day0 = ymd("2020-03-01")

# ending day
daymax = ymd("2020-09-01")

# to use in the model
tmax = as.numeric(difftime(daymax, day0, units="days"))

# dayseq is the calendar day sequence we want to model
dayseq = seq(day0, daymax, by="day")

# the month sequences are for display only
monthseq = seq(day0, daymax, by="month")
monthseq_lab = format(monthseq, "%b %Y")
daymonthseq = difftime(monthseq, day0, units="days")




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


#pop.prop <- init$population/sum(init$population) # proportion of state populaiton in each county


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

#################
# run the model 

sir_result = run_sir_model(state0=state0, params=params_init, region_adj=adj, populations=init$population, tmax=tmax)



compartment_plot_labels = c("D")
compartment_plot_names = c("Deaths")
compartment_plot_colors = rainbow(length(compartment_plot_labels))


plot_ct_region = function(region_name) {

  par(mar=c(4,4,3,4), bty="n")

  plot(sir_result$time, sir_result[,paste("I_s",region_name,sep=".")], type="n", 
       xlab="Time", ylab="People", main=region_name, col="red", 
       ylim=c(0,max(nyt_ct$deaths)), xlim=c(0,1.1*tmax), 
       axes=FALSE)
  axis(1,at=daymonthseq, lab=monthseq_lab)
  axis(2)

	for(j in 1:length(compartment_plot_labels)) {
    cidx = paste(compartment_plot_labels[j], region_name, sep=".")
    lines(sir_result$time, sir_result[,cidx], col=compartment_plot_colors[j])
	  text(sir_result$time[tmax+1], sir_result[tmax+1,cidx], format(sir_result[tmax+1,cidx],digits=1), pos=4, col=compartment_plot_colors[j])
	}

  if(region_name == "Connecticut") {
    points(nyt_ct$time, nyt_ct$deaths, pch=16, col=rgb(1,0,0,alpha=0.5)) 
  } else {
    obs.region <- subset(nyt, county == region_name)
    obs.region <- merge(obs.region, time.date, by='date')
    points(obs.region$time, obs.region$deaths, pch=16, col=rgb(1,0,0,alpha=0.5))
  }

	legend(0, max(nyt_ct$deaths), compartment_plot_names, lty=1, col=compartment_plot_colors, bty="n")

  # return some useful info
  # capacity exceeded?
  # describe intvx
	region_summary = paste("On ", format(daymax, "%b %d, %Y"),
											 " projections show ", format(sir_result[tmax,paste("D",region_name,sep=".")], digits=2),
											 " deaths in ", region_name,
                       sep="")

  return(region_summary)
}




########################

# make correspondence with calendar dates

# plot state overall, and by region  


par(mfrow=c(3,3))

plot_ct_region("Connecticut")
sapply(region_names, plot_ct_region)



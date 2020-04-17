library(deSolve)
library(yaml)
library(lubridate)

#####################

source("model.R")
source("get_data.R")

#####################


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
adj = read.csv("../map/CT_adj_matrix.csv")
rownames(adj) = adj$X
adj = adj[,-1]
adj = as.matrix(adj)

# populations and initial conditions
nregions = nrow(adj)
region_names = rownames(adj)

init <- read.csv('../data/ct_init.csv') 
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
       xlab="Time", ylab="People", main=region_name, col="black", 
       ylim=c(0,max(dat_ct_state$deaths)), xlim=c(0,1.1*tmax), 
       axes=FALSE)
  axis(1,at=daymonthseq, lab=monthseq_lab)
  axis(2)

	for(j in 1:length(compartment_plot_labels)) {
    cidx = paste(compartment_plot_labels[j], region_name, sep=".")
    lines(sir_result$time, sir_result[,cidx], col=compartment_plot_colors[j])
	  text(sir_result$time[tmax+1], sir_result[tmax+1,cidx], format(sir_result[tmax+1,cidx],digits=1), pos=4, 
	       col=compartment_plot_colors[j])
	}

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


	legend(0, max(dat_ct_state$deaths), compartment_plot_names, lty=1, col=compartment_plot_colors, bty="n")

  # return some useful info
	region_summary = paste("In ", region_name,
                       " on date ", tmax,
											 " projections show ", format(sir_result[tmax,paste("D",region_name,sep=".")], digits=2),
											 " deaths.", sep="")
}



#obs.region <- subset(dat_ct_county, county == "New Haven")
########################

# make correspondence with calendar dates

# plot state overall, and by region  


par(mfrow=c(3,3))

plot_ct_region("Connecticut")
sapply(region_names, plot_ct_region)



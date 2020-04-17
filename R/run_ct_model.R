library(deSolve)
library(yaml)

#####################

source("model.R")

#####################

######### observations #########

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
pop.prop <- init$population/sum(init$population) # proportion of state populaiton in each county


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
sir_result = run_sir_model(state0=state0, params=params_init, region_adj=adj, populations=init$population, tmax=100)

# get aggregate state-level projections
sir_result$S <- rowSums(sir_result[,2:9])
sir_result$E <- rowSums(sir_result[,10:17])
sir_result$Is <- rowSums(sir_result[,18:25])
sir_result$Im <- rowSums(sir_result[,26:33])
sir_result$A <- rowSums(sir_result[,34:41])
sir_result$H <- rowSums(sir_result[,42:49])
sir_result$Hbar <- rowSums(sir_result[,50:57])
sir_result$D <- rowSums(sir_result[,58:65])
sir_result$R <- rowSums(sir_result[,66:73])

ct.proj <- subset(sir_result, select=c(time, S, E, Is, Im, A, H, Hbar, D, R))

# compute state-level daily and cumulative infections and hospitalizations 
ct.proj$dailyI <- params_init$delta * ct.proj$E
ct.proj$dailyH <- params_init$alpha * ct.proj$Is

ct.proj$cum_modI <- cumsum(ct.proj$dailyI)
ct.proj$cum_modH <- cumsum(ct.proj$dailyH)



#### plot state-wide and county-level cumulative deaths ####

# par for this plot
quartz(width=7.4, height=6.6)
layout(matrix(1:9,nrow=3,ncol=3,byrow=T), widths=c(1,1,1), heights=c(1,1,1))
par (oma=c(0,0,0,0))
par (mar=c(4.1,4.1,3.1,2.1))
#layout.show(9)

dur <- 45 # upper bound of time for plotting
ymax.ct <- 1.2*nyt_ct$deaths[nrow(nyt_ct)] # for state-level plotting

plot(x=ct.proj$time, y=ct.proj$D, type='l', ylim=c(0,ymax.ct), xlim=c(0,dur), main='deaths, state-wide CT', 
     xlab='days', ylab='cum. deaths')
par(new=T)
plot(x=nyt_ct$time, y=nyt_ct$deaths, type='p', col='red', ylim=c(0,ymax.ct), xlim=c(0,dur), xlab='', ylab='')

for(i in 1:nregions) {

  region = region_names[i]
  
   obs.region <- subset(nyt, county == region)
   obs.region <- merge(obs.region, time.date, by='date')
   ymax.region <- 1.2*obs.region$deaths[nrow(obs.region)]
   
  compartment = "D"
  idx = paste(compartment, region, sep=".")
  plot(sir_result$time, sir_result[,idx], type="l", xlab="days", ylab='cum. deaths', main=idx, ylim=c(0,ymax.region), xlim=c(0,dur))
  par(new=T)
  plot(x=obs.region$time, y=obs.region$deaths, type='p', col='red', ylim=c(0,ymax.region), xlim=c(0,dur), xlab='', ylab='')
}








# melt this? 
#compartment_labels = c("S", "I_s", "I_m", "A", "H", "Hbar", "D", "R")



par(mar=c(4,4,1,1), mfrow=c(2,4), bty="n")

for(i in 1:nregions) {

  region = region_names[i]

  compartment = "I_s"
  idx = paste(compartment, region, sep=".")
  plot(sir_result$time, sir_result[,idx], type="l", xlab="Time", ylab=idx)

}



########################

# make correspondence with calendar dates

# plot state overall, and by region  
















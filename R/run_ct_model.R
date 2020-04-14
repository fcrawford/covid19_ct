library(deSolve)
library(yaml)

#####################

source("model.R")

#####################


# load fixed parameters
params_init = yaml.load_file("params.yml") 

# load region adjacency matrix 
adj = read.csv("../map/CT_adj_matrix.csv")
rownames(adj) = adj$X
adj = adj[,-1]
adj = as.matrix(adj)

# populations? 

nregions = nrow(adj)

region_names = rownames(adj)


#########################
# set up initial state0

E_init = runif(nregions,min=0.01,max=0.15)
I_s_init = rep(0,nregions)
I_m_init = rep(0,nregions)
A_init = rep(0,nregions)
H_init = rep(0,nregions)
Hbar_init = rep(0,nregions)
D_init = rep(0,nregions)
R_init = rep(0,nregions)
S_init = 1 - (E_init + I_s_init + I_m_init + A_init + H_init + Hbar_init + D_init + R_init)

state0 = c(S=S_init, E=E_init, I_s=I_s_init, I_m=I_m_init, A=A_init, H=H_init, Hbar=Hbar_init, D=D_init, R=R_init)

#################
# run the model 

sir_result = run_sir_model(state0=state0, params=params_init, region_adj=adj, tmax=100)


# melt this? 

compartment_labels = c("S", "I_s", "I_m", "A", "H", "Hbar", "D", "R")



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
















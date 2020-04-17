
source("model.R")


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

tmax = 100

nsim = 10

region_names = c("Connecticut", region_names)
death_cols = paste0("D.", region_names, sep="")



sir_results = lapply(1:nsim, function(i){
  # sample new param values, or whatever
  params_init_i = params_init
  params_init_i$beta = rnorm(1,params_init$beta,sd=0.05)
  res = run_sir_model(state0=state0, params=params_init_i, region_adj=adj, populations=init$population, tmax=tmax)
  res$sim_id = i
  res
})

library(plyr)

# data frames stacked 
res = ldply(sir_results, rbind)

# get 5% and 95% quantiles for each variable across sim_ids 








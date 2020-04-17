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

tmax = 100

sir_result = run_sir_model(state0=state0, params=params_init, region_adj=adj, tmax=tmax)

# get whole-state result 

# make region-specific plotting function

# give it whole-state data to get whole-state results 

# melt this? 


compartment_plot_labels = c("S", "I_s", "I_m", "A", "H", "Hbar", "D", "R")
compartment_plot_names = c("Susceptible", "Severe infection", "Mild infection", "Asymptomatic", "Hospitalized", "Unhospitalized", "Deaths", "Recovered")
compartment_plot_colors = rainbow(length(compartment_plot_labels))


plot_ct_region = function(region_name) {

  par(mar=c(4,4,3,4), bty="n")

  plot(sir_result$time, sir_result[,paste("I_s",region_name,sep=".")], type="n", xlab="Time", ylab="People", main=region_name, col="red", ylim=c(0,max(sir_result[,-1])))

	for(j in 1:length(compartment_plot_labels)) {

    cidx = paste(compartment_plot_labels[j], region_name, sep=".")
    lines(sir_result$time, sir_result[,cidx], col=compartment_plot_colors[j])
	  text(sir_result$time[tmax], sir_result[tmax,cidx], format(sir_result[tmax,cidx],digits=1), pos=4, col=compartment_plot_colors[j])

	}

	legend(0.6*tmax, max(sir_result[,-1]), compartment_plot_names, lty=1, col=compartment_plot_colors, bty="n")

  # return some useful info
	region_summary = paste("In ", region_name,
                       " on date ", tmax,
											 " projections show ", format(sir_result[tmax,paste("D",region_name,sep=".")], digits=2),
											 " deaths.", sep="")
}




########################

# make correspondence with calendar dates

# plot state overall, and by region  
















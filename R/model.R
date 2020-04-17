
# state0 is the initial condition 
# params is a list of parameter values 
# region_adj is an adjacency matrix whose row/col labels are 

# need intervention coding 
# hosp capacity

run_sir_model = function(state0, params, region_adj, populations, tmax) {

  nregions = nrow(region_adj) # number of counties or towns, whatever. 

  if(length(state0) != nregions*9) stop("length of state0 must be nregions*9")

  region_names = rownames(region_adj)

  # indices 
  S_idx    = 1:nregions
  E_idx    = (  nregions+1):(2*nregions)
  I_s_idx  = (2*nregions+1):(3*nregions)
  I_m_idx  = (3*nregions+1):(4*nregions)
  A_idx    = (4*nregions+1):(5*nregions)
  H_idx    = (5*nregions+1):(6*nregions)
  Hbar_idx = (6*nregions+1):(7*nregions)
  D_idx    = (7*nregions+1):(8*nregions)
  R_idx    = (8*nregions+1):(9*nregions)


  # construct region-specific betas 
  # before/after

  beta_pre = params$beta_pre
  #beta_post = (1-params$effect_intvx)*beta_pre 
  
  # number of adjacent regions for each region: use to keep overall beta the same for each county
  region_adj_num <- rowSums(region_adj)
  
  beta_matrix_pre  = ( (1-params$k_n)*diag(1,nregions) + params$k_n*(1/region_adj_num)*region_adj ) * beta_pre

  params$m_Hbar = params$m_H * params$m_Hbar_mult

  #beta_matrix_post = ( (1-params$k_n)*diag(1,nregions) + params$k_n*region_adj ) * beta_post;


  model <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      S = state[S_idx]
      E = state[E_idx]
      I_s = state[I_s_idx]
      I_m = state[I_m_idx]
      A = state[A_idx]
      H = state[H_idx]
      Hbar = state[Hbar_idx]
      D = state[D_idx]
      R = state[R_idx]


      beta = beta_matrix_pre


      
      # slope: parameter that sets how close the sigmoids are to a true 0-1 switch
      
      # beta = (1- 1/(1+exp(slope*(time - intervention))))*(beta_matrix_pre- beta_matrix_post)+beta_matrix_post_pre; 
      
      # ^ smooth switch of the contact matrix at the intervention time using a function call
      
      ### Soheil 4/10: capacity overflow takes into account capacity increase
      
      Hospital_capacities_breached = 0 #(1- 1/(1+exp(slope*(H-(capacity * capacity_function(time)))))) 
      
      # ^ 67*1 vector telling us if hospitals in the counties are full (~ 0-1 output)
      
      ##################
      
      # k_A        relative transmission from asymptomatic infectives
      # delta      1/incubation period
      # q_A        % infected and asymptomatic
      # q_Im       % infected and mild
      # alpha      rate at which sever cases need hospitalization/care
      # gamma_A    recovery rate of asymptomatic
      # gamma_H    recovery rate of hosp
      # gamma_Hbar recovery rate of unhosp
      # m_H  r     % hospitalized that die
      # m_Hbar     % unhospitalized that die
      
      #################
  
      
      dS    <-            -S*( beta %*% (I_s + I_m + k_A*A)/populations )    # populations normalized contact rates
      
      dE    <-  -delta*E + S*( beta %*% (I_s + I_m + k_A*A)/populations )
      
      dI_s  <- (1 - q_Im - q_A)*delta*E - alpha*I_s
      
      dI_m  <- q_Im*delta*E - gamma_Im*I_m
      
      dA    <- q_A*delta*E - gamma_A*A
      
      dH    <- alpha*I_s*(1-Hospital_capacities_breached) - gamma_H*H
      
      dHbar <- alpha*I_s*Hospital_capacities_breached - gamma_Hbar*Hbar
      
      dD    <- gamma_H*m_H*H + gamma_Hbar*m_Hbar*Hbar
      
      dR    <- gamma_H*(1 - m_H)*H + gamma_Hbar*(1 - m_Hbar)*Hbar + gamma_Im*I_m + gamma_A*A

      return(list(c(dS,dE,dI_s,dI_m,dA,dH,dHbar,dD,dR)))
    })
  }


  ts = seq(0,tmax,by=1)

  out = as.data.frame(ode(y=state0, times=ts, func=model, parms=params, method="lsodes"))

  names(out) = c("time",
                 paste0("S.", region_names, sep=""),
                 paste0("E.", region_names, sep=""),
                 paste0("I_s.", region_names, sep=""),
                 paste0("I_m.", region_names, sep=""),
                 paste0("A.", region_names, sep=""),
                 paste0("H.", region_names, sep=""),
                 paste0("Hbar.", region_names, sep=""),
                 paste0("D.", region_names, sep=""),
                 paste0("R.", region_names, sep=""))

  return(out)
}



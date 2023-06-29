library(BDgraph)
library(pheatmap)
library(fda)
library(coda)

rm(list = ls()); 
setwd("C:/E/1_BFGM/Code -v7 replications/Helper_functions")
source('simulation_functions.R')
source('call_DBFGM.R')
source('MCMC_changepoint_DBFGM.R')
source('Helper_functions_mcmc.R')
source('performance_functions.R')


setwd("C:/E/1_BFGM/Code -v7 replications")
dir.create("Simulation_data")

### Generate data, 25 replications --------------------------------
## Set up parameters
n = 50   # sample size
p = 15;   # number of curves in each sample
T_data = 128 * 2;    # number of time points
changepoint_true = 128 + 1   # true changepoint used to generate the data

pii_local = 1
K_true = 5;   # number of basis used to generate the data
basis_type_true = 'polynomial'  # specify the basis used in data generation
continuous = TRUE    # continuous adjustment at the changepoint
sigma_epsilon_true = 0.05   # standard error of noise

## Generate synthetic data, one replication
set.seed(1)
data_changepoint = simulate_changepoint_data(p,
                                             K_true,
                                             n,
                                             T_data,
                                             pii_local,
                                             sigma_epsilon_true,
                                             basis_type_true,
                                             continuous = TRUE,
                                             changepoint_true)
folder_name = "Simulation_data"
file_name = paste("Simulation_data_changepoint.Rdata", sep = "")
save(data_changepoint, file=paste(folder_name, '/', file_name, sep = ""))
# Plot simulated data (Figure 3.1) -------------
pheatmap(data_changepoint$param_true$G_x_true[[1]] + 0, cluster_rows=FALSE, cluster_cols=FALSE, 
         color = c("white", "black"),breaks = c(0, 0.5, 1))
pheatmap(data_changepoint$param_true$G_x_true[[2]] + 0, cluster_rows=FALSE, cluster_cols=FALSE, 
         color = c("white", "black"),breaks = c(0, 0.5, 1))
matplot(data_changepoint$Y[1,,], type = 'l')

## Generate 25 replications 
nrep = 25
for (rep_ind in 1:nrep){
  random_seed = 12345 + rep_ind
  data_changepoint_rep = simulate_changepoint_data_replications(data_changepoint,
                                                                continuous = TRUE,
                                                                random_seed)
  folder_name = "Simulation_data"
  file_name = paste("Synthetic_data_rep", rep_ind, ".Rdata", sep = "")
  save(data_changepoint_rep, file=paste(folder_name, '/', file_name, sep = ""))
}


############### Run MCMC
### Model parameters ---------------------------
K = 5  # number of basis functions
basis_type = 'spline'   # type of basis functions
v0 = 0.02^2  # spike variance
h = 50^2;   
v1 = v0 * h   # slab variance
disp = TRUE    # display MCMC progress
a_pi = 10; b_pi = 40
changepoint_interval = c(changepoint_true - 20, changepoint_true + 20)  # constraint on the changepoint

### Run MCMC -------------------------------------------------------
## Set up MCMC
nburn = 3000; nsave = 2000; 
nrep = 25
folder_name = "Simulation_results_DBFGM"  # folder to save results
dir.create(folder_name)

rep_ind = 25
for (rep_ind in 1:2){ 
  print(rep_ind)
  
  ## Load data 
  folder_name = "Simulation_data"
  file_name = paste("Synthetic_data_rep", rep_ind, ".Rdata", sep = "")
  file=paste(folder_name, '/', file_name, sep = "")
  load(file)
  
  set.seed(12345 + rep_ind)
  ## Run
  mcmc_output = call_DBFGM(data_changepoint_rep,   
                                 K, basis_type,
                                 v0, v1, a_pi, b_pi,
                                 changepoint_interval,
                                 nburn, nsave,
                                 disp)
  ## Save results
  folder_name = "Simulation_results_DBFGM"
  file_name = paste("MCMC_output_DBFGM_rep", rep_ind, ".Rdata", sep = "")
  save(mcmc_output, file=paste(folder_name, '/', file_name, sep = ""))
  
  ## Delete more burn-in period 
  # mcmc_output = mcmc_output_DBFGM
  # nburn_perf = 3000
  # nsave_perf = 2000
  # mcmc_output$changepoint_save = mcmc_output_DBFGM$changepoint_save[(nburn_perf+1):(nburn_perf+nsave_perf)];
  # for (s_i in 1:2){
  #   mcmc_output[[s_i]]$C_save = mcmc_output_DBFGM[[s_i]]$C_save[,,(nburn_perf+1):(nburn_perf+nsave_perf)]
  #   mcmc_output[[s_i]]$adj_save = mcmc_output_DBFGM[[s_i]]$adj_save[,,(nburn_perf+1):(nburn_perf+nsave_perf)]
  #   mcmc_output[[s_i]]$Sig_save = mcmc_output_DBFGM[[s_i]]$Sig_save[,,(nburn_perf+1):(nburn_perf+nsave_perf)]
  #   mcmc_output[[s_i]]$pii_block_save = mcmc_output_DBFGM[[s_i]]$pii_block_save[,,(nburn_perf+1):(nburn_perf+nsave_perf)]
  #   mcmc_output[[s_i]]$B_save = mcmc_output_DBFGM[[s_i]]$B_save[,,(nburn_perf+1):(nburn_perf+nsave_perf)]
  #   mcmc_output[[s_i]]$sigma_epsilon_save = mcmc_output_DBFGM[[s_i]]$sigma_epsilon_save[(nburn_perf+1):(nburn_perf+nsave_perf)]
  # }
  ## Get model performances before and after the change point
  param_true = data_changepoint_rep$param_true
  for (s_i in 1:2){
    data_oneint = list()  # data and estimations inside one state/interval
    data_oneint$G_x_true = param_true$G_x_true[[s_i]]
    data_oneint$G_b_true = param_true$G_b_true[[s_i]]
    data_oneint$Omega_b_true = param_true$Omega_b_true[[s_i]]
    data_oneint$B_true = param_true$B_true[[s_i]]
    data_oneint$F_true = param_true$F_true
    data_oneint$K_true = param_true$K_true
    mcmc_output_oneint = mcmc_output[[s_i]]
    
    performance_graph = get_mcmc_perf_graph(mcmc_output_oneint, data_oneint, 
                                                       K, p, standardize = FALSE, block_thresh = 0.5, disp = FALSE)
    cat("\nGraph estimation performance of state ", s_i, "\n")
    print_mcmc_results(performance_graph, data_oneint)

    file_name = paste("MCMC_performance_DBFGM_s", s_i, "_rep", rep_ind, ".Rdata", sep = "")
    save(performance_graph, file=paste(folder_name, '/', file_name, sep = ""))
  }
  
}

## Run the partially separable functional graphical model (PS-FGM) ----------
# Given true change point
file.edit('Helper_functions/Simulation_PSFGM_rep.R')

## Run the functional Bayesian graphical lasso (BL-FGM-indp)-----------------
# given the true change point 
file.edit('Helper_functions/Simulations_BLFGM_rep.R')
# given full data, fit one model
file.edit('Helper_functions/Simulations_BLFGM_fulldata_rep.R')
file.edit('Helper_functions/Simulations_BLFGM_fulldata_rep_add_performance.R')

## Compare performance of all models --------------------
file.edit('Helper_functions/Simulations_analyze_results_rep.R')

library(BDgraph)
library(pheatmap)
library(fda)
library(coda)

rm(list = ls()); 
setwd("C:/E/1_BFGM/Code -v7 replications/Helper_functions")
source('simulation_functions.R')

setwd("C:/E/1_BFGM/Code -v7 replications/Replications")
dir.create("Simulation_data")

### Set up simulation parameters ---------------------------
## Parameters for data generation
n = 50   # sample size
p = 15;   # number of curves in each sample
T_data = 128 * 2;    # number of time points
changepoint_true = 128 + 1   # true changepoint used to generate the data

pii_local = 1
K_true = 5;   # number of basis used to generate the data
basis_type_true = 'polynomial'  # specify the basis used in data generation
continuous = TRUE    # continuous adjustment at the changepoint
sigma_epsilon_true = 0.05   # standard error of noise

### Generate synthetic data, one replication ---------------------
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

### Generate 25 replications ------------------------------------------------------
nrep = 25
for (rep_ind in 1:nrep){
  # set.seed(12345 + rep_ind)
  # data_changepoint_rn = simulate_jump_data_add_random_noise(p,
  #                                                           K_true,
  #                                                           n,
  #                                                           T_data,
  #                                                           pii_local,
  #                                                           sigma_epsilon_true,
  #                                                           basis_type_true,
  #                                                           continuous = TRUE,
  #                                                           changepoint_true,
  #                                                           data_changepoint)
  random_seed = 12345 + rep_ind
  data_changepoint_rep = simulate_changepoint_data_replications(data_changepoint,
                                                                continuous = TRUE,
                                                                random_seed)
  folder_name = "Simulation_data"
  file_name = paste("Synthetic_data_rep", rep_ind, ".Rdata", sep = "")
  save(data_changepoint_rep, file=paste(folder_name, '/', file_name, sep = ""))
}






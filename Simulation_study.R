library(BDgraph)
library(pheatmap)
library(fda)
library(coda)
getwd()
rm(list = ls()); 
setwd("Helper_functions")
source('Helper_functions_data_generation.R')
source('call_functions.R')
source('MCMC_algorithms.R')
source('Helper_functions_mcmc.R')
source('performance_functions.R')
setwd("..")
dir.create("Simulation_data")

### Generate data, 50 replications --------------------------------
## Set up parameters
n = 50   # sample size
p = 15;   # number of curves in each sample
T_data = 128 * 2;    # number of time points
changepoint_true = c(128 + 1)   # true changepoint used to generate the data

pii_local = 1
K_true = 5;   # number of basis used to generate the data
basis_type_true = 'polynomial'  # specify the basis used in data generation
continuous = TRUE    # continuous adjustment at the changepoint
sigma_epsilon_true = 0.05   # standard error of noise

## Generate synthetic data, one replication
set.seed(1)
data_changepoint = simulate_data(p,
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

## Plot simulated data (Figure 3.1) 
pheatmap(data_changepoint$param_true$G_x_true[[1]] + 0, cluster_rows=FALSE, cluster_cols=FALSE, 
         color = c("white", "black"),breaks = c(0, 0.5, 1))
pheatmap(data_changepoint$param_true$G_x_true[[2]] + 0, cluster_rows=FALSE, cluster_cols=FALSE, 
         color = c("white", "black"),breaks = c(0, 0.5, 1))
matplot(data_changepoint$Y[1,,], type = 'l')

## Generate 50 replications 
nrep = 50
for (rep_ind in 1:nrep){
  random_seed = 12345 + rep_ind
  data_changepoint_rep = simulate_data_replications(data_changepoint,
                                                    continuous = TRUE,
                                                    random_seed)
  folder_name = "Simulation_data"
  file_name = paste("Synthetic_data_rep", rep_ind, ".Rdata", sep = "")
  save(data_changepoint_rep, file=paste(folder_name, '/', file_name, sep = ""))
}


#### Run MCMC ---------------------------
### Model parameters
K = 5  # number of basis functions
basis_type = 'spline'   # type of basis functions
v0 = 0.02^2  # spike variance
h = 50^2;   
v1 = v0 * h   # slab variance
disp = TRUE    # display MCMC progress
a_pi = 2; b_pi = 7
changepoint_interval = matrix(NA, 1, 2)
changepoint_interval[1,] = c(109, 149)  

### Set up MCMC
nburn = 3000; nsave = 2000; 
nrep = 50
folder_name = "Simulation_results_DBFGM"  # folder to save results
dir.create(folder_name)

for (rep_ind in 1:nrep){ 
  print(rep_ind)
  
  ## Load data 
  folder_name = "Simulation_data"
  file_name = paste("Synthetic_data_rep", rep_ind, ".Rdata", sep = "")
  file=paste(folder_name, '/', file_name, sep = "")
  load(file)
  
  set.seed(12345 + rep_ind)
  ## Run MCMC for the proposed Dynamic Bayesian Function Graphical Model
  mcmc_output = call_DBFGM(data_changepoint_rep,
                            K,
                            nburn, nsave,
                            v0, v1, a_pi, b_pi,
                            basis_type,
                            changepoint_interval,
                            disp)
  ## Save results
  folder_name = "Simulation_results_DBFGM"
  file_name = paste("MCMC_output_DBFGM_rep", rep_ind, ".Rdata", sep = "")
  save(mcmc_output, file=paste(folder_name, '/', file_name, sep = ""))
  
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

## Run the functional Bayesian graphical lasso (BL-FGM)-----------------
# given the true change point 
file.edit('Helper_functions/Simulations_BLFGM_rep.R')
# given full data, fit one model
file.edit('Helper_functions/Simulations_BLFGM_fulldata_rep.R')
file.edit('Helper_functions/Simulations_BLFGM_fulldata_rep_add_performance.R')

## Run B-GGM ----------------
file.edit('Helper_functions/Simulations_BGGM_rep.R')

## Compare performance of all models --------------------
file.edit('Helper_functions/Simulations_analyze_results_rep.R')




########
## Sample from Prior -------------------
########
get_ind = function(p, K){
  p_all = p*K
  # get ind_upper_block and idx_upper
  ind_all = matrix(1:p_all^2, p_all, p_all)
  ind_upper_block = c()  # matrix index in and sorted by upper diagonal blocks
  for (i in 1:(p - 1)){
    for (j in (i+1): p){
      rows = ((i - 1) * K + 1) : (i * K)
      cols = ((j - 1) * K + 1) : (j * K)
      ind_upper_block = c(ind_upper_block, c(ind_all[rows, cols]))
    }
  } 
  return(ind_upper_block)
}

# Set value of parameters
p = 15;   # number of curves in each sample
K = 10
v0 = 0.02^2  # spike variance
h = 50^2;   
v1 = v0 * h   # slab variance
disp = TRUE    # display MCMC progress
a_pi = 2; b_pi = 7
print(paste0('prior mean of pi: ', a_pi/(a_pi + b_pi)))
print(paste0('prior variance of pi: ', a_pi*b_pi/(a_pi + b_pi)^2/(a_pi + b_pi + 1)))
x <- rbeta(10000, a_pi, b_pi)                              
hist(x, freq=F, col="grey", border="white",main="rbeta(10000, 2, 1)", xlab="x", ylab="f(x)")          
print(paste0('desired edge inclusion probability is: ', 2/(p-1)))
# compute desired b_pi
ind_temp = 2
b_pi = (a_pi - a_pi * e_p)/e_p
# MCMC initial values
pii_block = matrix(1/2, nrow = p, ncol = p)   
diag(pii_block) = 1
p_all = p*K
C = diag(p_all); adj = diag(TRUE, p_all)
Sig = C
## Sample from prior
nburn = 1000
nsave = 1000
disp = FALSE
set.seed(1)
output = MCMC_from_DBFGM_prior(p,K, 
                               v0, v1, 
                               a_pi, b_pi,
                               Sig, C, adj, pii_block,
                               nburn, nsave,
                               disp)
## Compute edge inclusion probability of the graph in the coefficient space
ppi_local = apply(output$adj_save, c(1,2), mean)
ind_upper = get_ind(p, K)
post_mean = mean(ppi_local[ind_upper])
p_output_vec = c(p_output_vec, post_mean)
print(paste0('mean of edge inclusion probability: ', mean(ppi_local[ind_upper])))
print(paste0('std of edge inclusion probability: ', sd(ppi_local[ind_upper])))



################
### Scalability study -----------------------------
################
### Generate data --------------------------------
## Set up parameters
n = 30   # sample size
p = 50;   # number of curves in each sample
T_data = 200 * 2;    # number of time points
changepoint_true = c(200 + 1)   # true changepoint used to generate the data

pii_local = 1
K_true = 5;   # number of basis used to generate the data
basis_type_true = 'polynomial'  # specify the basis used in data generation
continuous = TRUE    # continuous adjustment at the changepoint
sigma_epsilon_true = 0.05   # standard error of noise

## Generate synthetic data, one replication
set.seed(1)
data_changepoint = simulate_data(p,
                                 K_true,
                                 n,
                                 T_data,
                                 pii_local,
                                 sigma_epsilon_true,
                                 basis_type_true,
                                 continuous = TRUE,
                                 changepoint_true)
pheatmap(data_changepoint$param_true$G_x_true[[1]] + 0, cluster_rows=FALSE, cluster_cols=FALSE, 
         color = c("white", "black"),breaks = c(0, 0.5, 1))
pheatmap(data_changepoint$param_true$G_x_true[[2]] + 0, cluster_rows=FALSE, cluster_cols=FALSE, 
         color = c("white", "black"),breaks = c(0, 0.5, 1))
matplot(data_changepoint$Y[1,,5], type = 'l')

### Model parameters ---------------------------
K = 5  # number of basis functions
basis_type = 'spline'   # type of basis functions
v0 = 0.02^2  # spike variance
h = 50^2;   
v1 = v0 * h   # slab variance
disp = TRUE    # display MCMC progress
a_pi = 2; b_pi = 7
changepoint_interval = matrix(NA, 1, 2)
changepoint_interval[1,] = c(109, 149)  

### Run MCMC -------------------------------------
## Set up MCMC
nburn = 5000; nsave = 2000; 
set.seed(12345)
## Run
mcmc_output = call_DBFGM(data_changepoint,
                         K,
                         nburn, nsave,
                         v0, v1, a_pi, b_pi,
                         basis_type,
                         changepoint_interval,
                         disp)
## Save results
folder_name = "Simulation_results_DBFGM"
file_name = paste("Scalability_p50_n30.Rdata", sep = "")
save(mcmc_output, file=paste(folder_name, '/', file_name, sep = ""))

## Get model performances before and after the change point
param_true = data_changepoint$param_true
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
}



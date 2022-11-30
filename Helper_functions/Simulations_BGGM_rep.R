## SSSL on two time points

library(BDgraph)
library(pheatmap)
library(fda)
library(coda)

rm(list = ls()); 

setwd("C:/E/1_BFGM/Code -v7 replications/Helper_functions")
source('simulation_functions.R')
source('Helper_functions_mcmc.R')
source('performance_functions.R')
source('Call_functions_other.R')
source('MCMC_algorithms_other.R')

setwd("C:/E/1_BFGM/Code -v7 replications")


### Simulation parameters ---------------------------
nburn = 2000; nsave = 2000; 
## Hyper-parameters in the model 
p = 15
v0 = 0.02^2
h = 50^2;
v1 = v0 * h
disp = TRUE
pii = 3/(p-1)


### Set up MCMC ------------------------------------------------------
folder_name = "Simulation_results_BGGM"
dir.create(folder_name)


### Run MCMC -------------------------------------------------------

for (rep_ind in 1:25){  # 3:4
  print(rep_ind)
  
  mcmc_output_full = c(); mcmc_output = c()
  ### Change point Model ------------------------------------------------
  ## Load data 
  folder_name = "Simulation_data"
  file_name = paste("Synthetic_data_rep", rep_ind, ".Rdata", sep = "")
  file=paste(folder_name, '/', file_name, sep = "")
  load(file)
  changepoint_true = data_changepoint_rep$param_true$changepoint_true
  T_data = dim(data_changepoint_rep$Y)[2]
  t_list = c(256/4, 256/4*3)
  performance_graph = list()
  
  ## Run MCMC
  for (s_i in 1:2){
    Y = data_changepoint_rep$Y[,t_list[s_i],]
    set.seed(12345 + rep_ind)
    mcmc_output = call_SSSL(Y,
                            nburn, nsave, 
                            v0, v1, pii, lambda = 1,
                            disp = FALSE)
    
    ppi_block = apply(mcmc_output$adj_save, c(1,2), mean)
    adj_block = ppi_block > 0.5
    diag(adj_block) = TRUE
    
    temp = get_tp_fp_tn_fn(data_changepoint_rep$param_true$G_x_true[[s_i]]+0, adj_block+0)
    tpr = temp$tp / (temp$tp + temp$fn) # tpr
    fpr = temp$fp / (temp$fp + temp$tn) # fpr
    mcc = (temp$tp * temp$tn - temp$fp * temp$fn) /
      (sqrt(temp$tp + temp$fp) * sqrt(temp$tp + temp$fn) * sqrt(temp$tn + temp$fp) * sqrt(temp$tn + temp$fn))
    
    print(c(tpr,fpr,mcc))
    
  }
    
  folder_name = "Simulation_results_BGGM"
  file_name = paste("MCMC_performance_BGGM_rep", rep_ind, ".Rdata", sep = "")
  save(performance_graph, file=paste(folder_name, '/', file_name, sep = ""))
  
}

## Summary of all models --------------------------------
file.edit('Helper_functions/Simulations_analyze_results_rep.R')




library(BDgraph)
library(pheatmap)
library(fda)
library(coda)

rm(list = ls()); 
setwd("C:/E/1_BFGM/Code -v7 replications/Helper_functions")
source('simulation_functions.R')
source('performance_functions.R')
setwd("C:/E/1_BFGM/Code -v7 replications")

### Load results BDFGM--------------
nrep = 25
tpr = matrix(NA, nrep, 2); fpr = matrix(NA, nrep, 2); mcc = matrix(NA, nrep, 2)
for (rep_ind in 1:nrep){
  for (s_i in 1:2){
    folder_name = "Simulation_results_DBFGM"
    file_name = paste("MCMC_performance_DBFGM_s", s_i, "_rep", rep_ind, ".Rdata", sep = "")
    load(file=paste(folder_name, '/', file_name, sep = ""))
    tpr[rep_ind, s_i] = performance_graph$block_perf_v1$tpr_block
    fpr[rep_ind, s_i] = performance_graph$block_perf_v1$fpr_block
    mcc[rep_ind, s_i] = performance_graph$block_perf_v1$mcc_block
    
  }
}
apply(tpr, 2, mean)
apply(tpr, 2, sd)
apply(fpr, 2, mean)
apply(fpr, 2, sd)
apply(mcc, 2, mean)
apply(mcc, 2, sd)


### Load results PSFGM----------------------
# Simulation_PSFGM_rep.R
folder_name = "Simulation_results_PSFGM/"
nrep = 25
tpr = matrix(NA, nrep, 2); fpr = matrix(NA, nrep, 2); mcc = matrix(NA, nrep, 2)
for (rep_ind in 1:nrep){
  for (s_i in 1:2){
    file_name = paste("Output_PSFGM_rep", rep_ind, ".Rdata", sep = "")
    load(file=paste(folder_name, '/', file_name, sep = ""))
    tpr[rep_ind, s_i] = psfgm_output$performance[[s_i]][1]
    fpr[rep_ind, s_i] =  psfgm_output$performance[[s_i]][2]
    mcc[rep_ind, s_i] =  psfgm_output$performance[[s_i]][3]
    
  }
}
apply(tpr, 2, mean)
apply(tpr, 2, sd)
apply(fpr, 2, mean)
apply(fpr, 2, sd)
apply(mcc, 2, mean)
apply(mcc, 2, sd)


# Load results BLFGM ------------------
# Given change point
folder_name = "Simulation_results_BLFGM"
nrep = 25
tpr = matrix(NA, nrep, 2); fpr = matrix(NA, nrep, 2); mcc = matrix(NA, nrep, 2)
for (rep_ind in 1:nrep){
  file_name = paste("MCMC_output_BLFGM_alldata_rep", rep_ind, ".Rdata", sep = "")
  load(file=paste(folder_name, '/', file_name, sep = ""))
  for (s_i in 1:2){
    tpr[rep_ind, s_i] = blfgm_output$performance[[s_i]][1]
    fpr[rep_ind, s_i] =  blfgm_output$performance[[s_i]][2]
    mcc[rep_ind, s_i] =  blfgm_output$performance[[s_i]][3]
  }
}
# Given full data
# nrep = 25
# tpr = matrix(NA, nrep, 1); fpr = matrix(NA, nrep, 1); mcc = matrix(NA, nrep, 1)
# for (rep_ind in 1:nrep){
#     file_name = paste("MCMC_output_BLFGM_alldata_rep", rep_ind, ".Rdata", sep = "")
#     load(file=paste(folder_name, '/', file_name, sep = ""))
#     tpr[rep_ind] = blfgm_output$performance[1]
#     fpr[rep_ind] =  blfgm_output$performance[2]
#     mcc[rep_ind] =  blfgm_output$performance[3]
#   
# }
apply(tpr, 2, mean)
apply(tpr, 2, sd)
apply(fpr, 2, mean)
apply(fpr, 2, sd)
apply(mcc, 2, mean)
apply(mcc, 2, sd)


# BGGM -----------------
nrep = 25
tpr = matrix(NA, nrep, 2); fpr = matrix(NA, nrep, 2); mcc = matrix(NA, nrep, 2)
for (rep_ind in 1:nrep){
  folder_name = "Simulation_results_BGGM"
  file_name = paste("MCMC_performance_BGGM_rep", rep_ind, ".Rdata", sep = "")
  load(file=paste(folder_name, '/', file_name, sep = ""))

  for (s_i in 1:2){
    tpr[rep_ind, s_i] = performance_graph[[s_i]][1]
    fpr[rep_ind, s_i] =  performance_graph[[s_i]][2]
    mcc[rep_ind, s_i] =  performance_graph[[s_i]][3]
  }
  
}

apply(tpr, 2, mean)
apply(tpr, 2, sd)
apply(fpr, 2, mean)
apply(fpr, 2, sd)
apply(mcc, 2, mean)
apply(mcc, 2, sd)




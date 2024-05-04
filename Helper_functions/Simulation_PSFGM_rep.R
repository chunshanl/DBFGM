library(fgm)
library(mvtnorm)
library(fda)

rm(list = ls()); 
setwd("Helper_functions")
source('performance_functions.R')
setwd("..")
folder_name = "Simulation_results_PSFGM"  # folder to save results
dir.create(folder_name)

### Run MCMC -------------------------------------------------------
for (rep_ind in 1:50){
  
  print(rep_ind)
  psfgm_output = list()
  
  ## Load data 
  folder_name = "Simulation_data"
  file_name = paste("Synthetic_data_rep", rep_ind, ".Rdata", sep = "")
  file=paste(folder_name, '/', file_name, sep = "")
  load(file)
  
  for (s_i in 1:2){
    # Get data
    data = changepoint_to_static(data_changepoint_rep, s_i)
    p = dim(data$Y)[3]; T_data = dim(data$Y)[2]; n = dim(data$Y)[1]
    y <- list()
    for (j in 1:p){
      y[[j]] = data$Y[,,j]
    }
    # Run 
    set.seed(54321)
    #alpha_list = seq(0,1,length=10)  # tunning
    #gamma_list = seq(0,1,length=10)
    #for (ii in 1:10){
    #  for (jj in 1:10){
        G_est <- fgm(y, L = 5, alpha=0.3, gamma=0.7)$A
        temp = get_tp_fp_tn_fn(data$G_x_true+0, G_est+0)
        tpr = temp$tp / (temp$tp + temp$fn) # tpr
        fpr = temp$fp / (temp$fp + temp$tn) # fpr
        mcc = (temp$tp * temp$tn - temp$fp * temp$fn) /
          (sqrt(temp$tp + temp$fp) * sqrt(temp$tp + temp$fn) * sqrt(temp$tn + temp$fp) * sqrt(temp$tn + temp$fn))
   #     print(c(ii,jj))

    #  }
    #}
    # Organize results
    psfgm_output$G_est[[s_i]] = G_est
    psfgm_output$performance[[s_i]] = c(tpr, fpr,mcc)
    print(c(tpr, fpr,mcc))
  }
  
  ## Save results 
  folder_name = "Simulation_results_PSFGM"
  file_name = paste("Output_PSFGM_rep", rep_ind, ".Rdata", sep = "")
  save(psfgm_output, file=paste(folder_name, '/', file_name, sep = ""))
  
}

# add model performance to the model output blfgm_output

rm(list = ls()); 
library(wavelets)
setwd("Helper_functions/FunGraph-main")
#library(doParallel)
# library(foreach)
#registerDoParallel(24)
source('fgraph_cclasso.R')
source('Cmat_update.R')
source('Dmat_update.R')
source('Lam_update.R')
setwd("..")
source('performance_functions.R')
setwd("..")

### Set up MCMC
# Model parameters
MCMCspecs = list(B=1000,thin=5,burnin=1000,update=1000);
D_prior = list(a=0.1,b=0.1); lam_prior = list(a=0.1,b=1);

rep_ind = 1
### Run MCMC -------------------------------------------------------
for (rep_ind in 1:25){
  
  print(rep_ind)
  blfgm_output = list()
  
  ## Load data 
  folder_name = "Simulation_data"
  file_name = paste("Synthetic_data_rep", rep_ind, ".Rdata", sep = "")
  file=paste(folder_name, '/', file_name, sep = "")
  load(file)
  # Get data
  data = data_changepoint_rn
  p = dim(data$Y)[3]; T_data = dim(data$Y)[2]; n = dim(data$Y)[1]

  folder_name = "Simulation_results_BLFGM/"
  file_name = paste("MCMC_output_BLFGM_alldata_rep", rep_ind, ".Rdata", sep = "")
  load(file=paste(folder_name, '/', file_name, sep = ""))
  
  blfgm_output$performance = list()
  
  fgbay.est.all <- blfgm_output$fgbay.est
  # Estimation performance
  # before changepoint
  tmp <- matrix(FALSE,p,p); tmp[upper.tri(tmp)] <- TRUE
  uptri.ind <- array(tmp,c(p,p,data$param_true$changepoint_true-1))
  G1.true = array(NA, dim = c(p,p,data$param_true$changepoint_true-1)) 
  G1.neg = G1.true
  for (t in 1:(data$param_true$changepoint_true-1)){
    G1.true[,,t] = data$param_true$G_x_true[[1]] + 0
    G1.neg[,,t] = (!data$param_true$G_x_true[[1]])+ 0
  }
  fgbay.est = fgbay.est.all[,, 1:(data$param_true$changepoint_true-1) ]
  tp = sum((fgbay.est+G1.true==2)[uptri.ind])
  fp = sum((fgbay.est+G1.neg==2)[uptri.ind])
  tn = sum(((1-fgbay.est)+G1.neg==2)[uptri.ind])
  fn = sum(((1-fgbay.est)+G1.true==2)[uptri.ind])
  tpr = sum((fgbay.est+G1.true==2)[uptri.ind])/sum(G1.true[uptri.ind]) 
  fpr = sum((fgbay.est+G1.neg==2)[uptri.ind])/sum(G1.neg[uptri.ind])
  mcc = (tp * tn - fp * fn) / 
    (sqrt(tp + fp) * sqrt(tp + fn) * sqrt(tn + fp) * sqrt(tn + fn))
  print(c(tpr,fpr,mcc))
  blfgm_output$performance[[1]] = c(tpr, fpr,mcc)
  
  # after
  tmp <- matrix(FALSE,p,p); tmp[upper.tri(tmp)] <- TRUE
  uptri.ind <- array(tmp,c(p,p,T_data - data$param_true$changepoint_true+1))
  G1.true = array(NA, dim = c(p,p,T_data - data$param_true$changepoint_true+1)) 
  G1.neg = G1.true
  for (t in 1:(T_data - data$param_true$changepoint_true+1)){
    G1.true[,,t] = data$param_true$G_x_true[[2]] + 0
    G1.neg[,,t] = (!data$param_true$G_x_true[[2]])+ 0
  }
  fgbay.est = fgbay.est.all[,, (data$param_true$changepoint_true):T_data ]
  tp = sum((fgbay.est+G1.true==2)[uptri.ind])
  fp = sum((fgbay.est+G1.neg==2)[uptri.ind])
  tn = sum(((1-fgbay.est)+G1.neg==2)[uptri.ind])
  fn = sum(((1-fgbay.est)+G1.true==2)[uptri.ind])
  tpr = sum((fgbay.est+G1.true==2)[uptri.ind])/sum(G1.true[uptri.ind]) 
  fpr = sum((fgbay.est+G1.neg==2)[uptri.ind])/sum(G1.neg[uptri.ind])
  mcc = (tp * tn - fp * fn) / 
    (sqrt(tp + fp) * sqrt(tp + fn) * sqrt(tn + fp) * sqrt(tn + fn))
  print(c(tpr,fpr,mcc))
  blfgm_output$performance[[2]] = c(tpr, fpr,mcc)

  ## Save results 
  folder_name = "Simulation_results_BLFGM/"
  file_name = paste("MCMC_output_BLFGM_alldata_rep", rep_ind, ".Rdata", sep = "")
  save(blfgm_output, file=paste(folder_name, '/', file_name, sep = ""))
  
}

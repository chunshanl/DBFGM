# BL-FGM-indp in paper
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

folder_name = "Simulation_results_DBFGM"  # folder to save results
dir.create(folder_name)

### Set up MCMC
# Model parameters
MCMCspecs = list(B=1000,thin=5,burnin=1000,update=1000);
D_prior = list(a=0.1,b=0.1); lam_prior = list(a=0.1,b=1);
nrep = 25
# storage for the model performance measurements
performance_all = list(graph_est = vector("list", length = nrep))

### Run MCMC -------------------------------------------------------

for (rep_ind in 1:nrep){
  
  print(rep_ind)
  blfgm_output = list()
  
  ## Load data 
  folder_name = "Simulation_data"
  file_name = paste("Synthetic_data_rep", rep_ind, ".Rdata", sep = "")
  file=paste(folder_name, '/', file_name, sep = "")
  load(file)
  
  for (s_i in 1:2){
    # Get data
    data = changepoint_to_static(data_changepoint_rep, s_i)
    p = dim(data$Y)[3]; T_data = dim(data$Y)[2]; n = dim(data$Y)[1]
    y <- array(NA,c(n,p,T_data))
    for (i in 1:n){
      for (j in 1:p){
        y[i,j,] = data$Y[i,,j]
      }
    }
    # Generate basis
    tmp = diag(T_data)
    W = matrix(NA,T_data,T_data)
    for(t in 1:T_data) {
      a = dwt(tmp[t,],filter='haar',n.levels=5)
      W[t,] = c(a@V[[5]],unlist(a@W[5:1]))
    }
    Phi <- solve(W)
    # Run MCMC
    set.seed(54321)
    fgbay.samp <- fgraph_ccc(y,Phi,MCMCspecs,D_prior,lam_prior)
    # Graph estimation
    fgbay.est <- array(diag(p),c(p,p,T_data))
    for(t in 1:T_data) fgbay.est[,,t][upper.tri(diag(p))] <- apply(fgbay.samp$ccc[,t,],1,function(x) (quantile(x,0.025)>0 | quantile(x,0.975)<0)*1)
    # Estimation performance
    tmp <- matrix(FALSE,p,p); tmp[upper.tri(tmp)] <- TRUE
    uptri.ind <- array(tmp,c(p,p,T_data))
    G1.true = array(NA, dim = c(p,p,T_data)) 
    G1.neg = G1.true
    for (t in 1:T_data){
      G1.true[,,t] = data$G_x_true + 0
      G1.neg[,,t] = !data$G_x_true + 0
    }
    tp = sum((fgbay.est+G1.true==2)[uptri.ind])
    fp = sum((fgbay.est+G1.neg==2)[uptri.ind])
    tn = sum(((1-fgbay.est)+G1.neg==2)[uptri.ind])
    fn = sum(((1-fgbay.est)+G1.true==2)[uptri.ind])
    tpr = sum((fgbay.est+G1.true==2)[uptri.ind])/sum(G1.true[uptri.ind]) 
    fpr = sum((fgbay.est+G1.neg==2)[uptri.ind])/sum(G1.neg[uptri.ind])
    mcc = (tp * tn - fp * fn) / 
      (sqrt(tp + fp) * sqrt(tp + fn) * sqrt(tn + fp) * sqrt(tn + fn))
    # Organize results
    blfgm_output$fgbay.samp[[s_i]] = fgbay.samp
    blfgm_output$fgbay.est[[s_i]] = fgbay.est
    blfgm_output$performance[[s_i]] = c(tpr, fpr,mcc)
    print(c(tpr, fpr,mcc))
  }

  ## Save results 
  folder_name = "Simulation_results_BLFGM/"
  file_name = paste("MCMC_output_BLFGM_rep", rep_ind, ".Rdata", sep = "")
  save(blfgm_output, file=paste(folder_name, '/', file_name, sep = ""))
  
}

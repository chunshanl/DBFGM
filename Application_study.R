library(BDgraph)
library(pheatmap)
library(fda)
getwd()
rm(list = ls()); 
setwd("Helper_functions")
source('Helper_functions_data_generation.R')
source('call_functions.R')
source('MCMC_algorithms.R')
source('Helper_functions_mcmc.R')
source('performance_functions.R')
setwd("..")

### Get sst data
# script to process raw .nc data
#file.edit('Helper_functions/sst_data_processing.R')
# load data
load(file = 'sst_processed_update.Rdata')

#### Run the proposed Bayesian functional graphical model, static version -------------------
nburn = 3000; nsave = 2000; 
# Parameters
K = 10  # number of basis functions
basis_type = 'spline'   # type of basis functions
v0 = 0.02^2  # spike variance
h = 50^2;   
v1 = v0 * h   # slab variance
disp = TRUE    # display MCMC progress
a_pi = 10; b_pi = 40
disp = T

set.seed(12345)
mcmc_output_static = call_DBFGM_static(data,
                                       K,
                                       nburn, nsave,
                                       v0, v1, a_pi, b_pi,
                                       basis_type,
                                       disp)
### Save results
folder_name = "Application_results"
file_name = paste("MCMC_output_static.Rdata", sep = "")
save(mcmc_output_static, file=paste(folder_name, '/', file_name, sep = ""))

###### Run the proposed dynamic Bayesian functional graphical model --------------------
nburn = 3000; nsave = 2000; 
# Parameters
K = 10  # number of basis functions
basis_type = 'spline'   # type of basis functions
#basis_type = 'fourier'
v0 = 0.02^2  # spike variance
h = 50^2;   
v1 = v0 * h   # slab variance
disp = TRUE    # display MCMC progress
a_pi = 10; b_pi = 40
changepoint_interval = matrix(NA, 1, 2)
changepoint_interval[1,] = c(160, 200)  
disp = T
# Run MCMC
set.seed(123)
mcmc_output_DBFGM = call_DBFGM(data,
                               K,
                               nburn, nsave,
                               v0, v1, a_pi, b_pi,
                               basis_type,
                               changepoint_interval,
                               disp)
# Save results
folder_name = "Application_results"
file_name = paste("MCMC_output_DBFGM_10K.Rdata", sep = "")
save(mcmc_output_DBFGM, file=paste(folder_name, '/', file_name, sep = ""))


#### Run DBFGM, four seasons, three change points -----------------------------
nburn = 3000; nsave = 2000; 
# Parameters
K = 10  # number of basis functions
basis_type = 'spline'   # type of basis functions
v0 = 0.02^2  # spike variance
h = 50^2;   
v1 = v0 * h   # slab variance
disp = TRUE    # display MCMC progress
a_pi = 10; b_pi = 40
changepoint_interval = matrix(NA, 3, 2)
for (s_i in 1:3){
  changepoint_interval[s_i,] = c(364/4 * s_i-20, 364/4 * s_i+20)
}
changepoint_interval[2,] = c(140, 160)
set.seed(123)
mcmc_output_fourseasons = call_DBFGM(data,
                                     K,
                                     nburn, nsave,
                                     v0, v1, a_pi, b_pi,
                                     basis_type,
                                     changepoint_interval,
                                     disp)

### Save results
folder_name = "Application_results"
file_name = paste("MCMC_output_DBFGM_fourseasons.Rdata", sep = "")
save(mcmc_output_fourseasons, file=paste(folder_name, '/', file_name, sep = ""))

### Analyze MCMC outputs -------------------
file.edit('Helper_functions/Application_study_analyze_results.R')






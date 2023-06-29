library(BDgraph)
library(pheatmap)
library(fda)


setwd("./Helper_functions")
source('Helper_functions_mcmc.R')
#source('performance_functions.R')
#source('Helper_functions_sst.R')
source('call_DBFGM.R')
source('MCMC_changepoint_DBFGM.R')
source('Call_functions_other.R')
source('MCMC_algorithms_other.R')
setwd("..")

### Get sst data
# script to process raw .nc data
file.edit('Helper_functions/sst_data_processing.R')
load(file = 'sst_processed_update.Rdata')


#### Run the proposed Bayesian functional graphical model, static version -------------------
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
file_name = paste("MCMC_output_spline_static_update.Rdata", sep = "")
#file_name = paste("MCMC_output_fourier_static_5K.Rdata", sep = "")
save(mcmc_output_static, file=paste(folder_name, '/', file_name, sep = ""))

###### Run the proposed dynamic Bayesian functional graphical model --------------------
nburn = 3000; nsave = 2000; 
# Parameters
K = 6  # number of basis functions
basis_type = 'spline'   # type of basis functions
#basis_type = 'fourier'
v0 = 0.02^2  # spike variance
h = 50^2;   
v1 = v0 * h   # slab variance
disp = TRUE    # display MCMC progress
a_pi = 10; b_pi = 40
changepoint_interval = c(160, 200)  # constraint on the changepoint
disp = T
# Run MCMC
set.seed(123)
mcmc_output_DBFGM = call_DBFGM(data,   
                               K, basis_type,
                               v0, v1, a_pi, b_pi,
                               changepoint_interval,
                               nburn, nsave,
                               disp)
# Save results
folder_name = "Application_results"
file_name = paste("MCMC_output_spline_DBFGM_update.Rdata", sep = "")
# file_name = paste("MCMC_output_fourier_DBFGM.Rdata", sep = "")
save(mcmc_output_DBFGM, file=paste(folder_name, '/', file_name, sep = ""))


#### Run DBFGM, multiple change points -----------------------------
# Four seasons, three change points, constrained change point intervals
nburn = 3000; nsave = 2000; 
# Parameters
K = 5  # number of basis functions
basis_type = 'spline'   # type of basis functions
#basis_type = 'fourier'
v0 = 0.02^2  # spike variance
h = 50^2;   
v1 = v0 * h   # slab variance
disp = TRUE    # display MCMC progress
a_pi = 10; b_pi = 40
changepoint_interval = matrix(NA, 3, 2)
for (s_i in 1:3){
  changepoint_interval[s_i,] = c(364/4 * s_i-10, 364/4 * s_i+10)
}
changepoint_interval[2,] = c(140, 160)
set.seed(123)
mcmc_output_fourseasons = call_DBFGM_fourseasons(data,
                                                 K,
                                                 nburn, nsave,
                                                 v0, v1, a_pi, b_pi,
                                                 basis_type,
                                                 changepoint_interval,
                                                 disp)

### Save results
folder_name = "Application_results"
file_name = paste("MCMC_output_DBFGM_fourseasons_update.Rdata", sep = "")
save(mcmc_output_fourseasons, file=paste(folder_name, '/', file_name, sep = ""))

### Analyze MCMC outputs -------------------
file.edit('Helper_functions/Application_study_analyze_results.R')






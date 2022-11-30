library(anytime)
library(maps)
library(dplyr)
library(geosphere)
library(BDgraph)
library(pheatmap)
library(codetools)
library(fda)
library(coda)
library(fdapace)
library(matrixcalc)

### Load data --------------------------
# gctorture(FALSE)
rm(list = ls()); 
setwd("C:/E/1_BFGM/Code -v7 replications/Helper_functions")
source('simulation_functions.R')
source('Helper_functions_mcmc.R')
source('performance_functions.R')
source('Helper_functions_sst.R')
source('call_DBFGM.R')
source('MCMC_changepoint_DBFGM.R')
source('Call_functions_other.R')
source('MCMC_algorithms_other.R')

setwd("C:/E/1_BFGM/Code -v7 replications/")

## Read and process data -------
sst = read.csv('sst.csv')
## Select cols
p = 15
cols = round(seq(from = 2, to = ncol(sst), length.out = p))
cols[12] = cols[12]+1
## Get locations 
data_longitude = names(sst)[cols]  # data_longitude
temp = gsub('X.', '', data_longitude)
data_longitude = as.numeric(substr(temp,1,regexpr("\\.",temp)-1))
data_latitude = as.numeric(sst[1, cols] )  # data_latitude
locations = cbind(data_longitude, data_latitude)  # locations
save(locations, file = 'sst_locations.Rdata')
# # plot locations
# par(mar=c(0,0,0,0))
# map('world',
#     col="#f2f2f2", fill=TRUE, bg="white", lwd=0.05,
#     border=0,
# )
# points(x=locations[,1], y=locations[,2], col="slateblue", cex=0.5, pch=20)
## Get time
data_time = sst[-c(1:2),1]
data_time = anytime(data_time)
## Get temperature
sst_input = sst[-(1:2),cols]   # sst_input

## Moving average
smooth_data = TRUE
if (smooth_data){
  sst_input_smoothed = sst_input
  for (j in 1:p){
    sst_input_smoothed[[j]] = moving_average(sst_input[[j]], window_length = 13)
  }
  #plot(sst_input_smoothed[1:366,1],type = 'l')
  #lines(sst_input[1:366,1],type = 'l', col = 2)
  sst_input = sst_input_smoothed
}

## Reshape data into n*T*p
min(data_time)
max(data_time)
n = 2021-1979+1
T_data = 366
cur_year = format(data_time, '%Y')  # char
col_num = as.numeric(cur_year) - 1978
cur_day = (data_time - anytime(cur_year))/60/60/24
row_num = round(as.numeric(cur_day))

Y = array(0, dim = c(n, T_data, p))
for (i in 1:nrow(sst_input)){
  Y[col_num[i], row_num[i],] = as.numeric(sst_input[i,])
}
Y = Y[,-(365:366),]

# Start from March
Y_chop = array(NA, c(42,364,15))
ind_start = (31+29)
ind_end = 364*43 - (364 - 31 - 28)
for(j in 1:15){
  y_vec = c(t(Y[,,j]))
  y_vec = y_vec[ind_start:ind_end]
  y_matrix = matrix(y_vec, nrow = 364, ncol = 42)
  Y_chop[,,j]=t(y_matrix)
  
}
Y = Y_chop

# Standardize data
Y_sd = Y
n = dim(Y)[1]
for (i in 1:n){
  for(j in 1:p){
    Y_sd[i,,j] = (Y[i,,j] - mean(Y[i,,j]))/sd(Y[i,,j])
  }
}
data = list()
data$Y = Y_sd

save(data, file= 'sst_processed.Rdata')

## Plot data -----------------------------------------
i = 1
j = 1
matplot(t(Y[,,j]), type = 'l')
# xlab_temp = paste('Day of the year. Year of ', 1979 + j - 1)
# main_temp = paste('Longitude =', locations[j,1], ', Latitude=', locations[j,2])
# plot(sst_input[,j],xlab = '1979 - 2021', ylab = 'Temperature', main = main_temp,type = 'l')
# # plot(sst_input[1:364,j],xlab = '1979 - 2021', ylab = 'Temperature', main = main_temp)
#plot(data$Y[i,,j], xlab = '', ylab = 'Temperature', type = 'l')
#matplot(t(Y[,,j]), type = 'l',xlab = 'Day of the year, 1979 - 2021', ylab = 'Temperature')
# #plot(Y[5,,4],type = 'l')
# #plot(data_time[(365*4):(365*5)], y = sst[(365*4):(365*5),cols[4]], type = 'l')
# #plot(data_time[368:(368+364+10000)], y = sst[368:(368+364+10000),26], type = 'l')
# #matplot(t(Y[sample(1:n, n),,15]),type = 'l')

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
file_name = paste("MCMC_output_spline_static.Rdata", sep = "")
#file_name = paste("MCMC_output_fourier_static_5K.Rdata", sep = "")
save(mcmc_output_static, file=paste(folder_name, '/', file_name, sep = ""))

###### Run the proposed dynamic Bayesian functional graphical model --------------------
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
changepoint_interval = c(1, 364)  # constraint on the changepoint
disp = T
# Run MCMC
mcmc_output_DBFGM = call_DBFGM(data,   
                               K, basis_type,
                               v0, v1, a_pi, b_pi,
                               changepoint_interval,
                               nburn, nsave,
                               disp)
# Save results
folder_name = "Application_results"
file_name = paste("MCMC_output_spline_DBFGM.Rdata", sep = "")
#file_name = paste("MCMC_output_fourier_DBFGM.Rdata", sep = "")
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
set.seed(12345)
mcmc_output_fourseasons = call_DBFGM_fourseasons(data,
                                            K,
                                            nburn, nsave,
                                            v0, v1, a_pi, b_pi,
                                            basis_type,
                                            changepoint_interval,
                                            disp)

### Save results
folder_name = "Application_results"
file_name = paste("MCMC_output_DBFGM_fourseasons.Rdata", sep = "")
save(mcmc_output_full, file=paste(folder_name, '/', file_name, sep = ""))

### Analyze MCMC outputs -------------------
file.edit('Helper_functions/Application_study_analyze_results.R')




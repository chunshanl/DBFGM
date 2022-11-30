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

### Load data --------------------------
load('sst_processed.Rdata')
n = dim(data$Y)[1]; T_data = dim(data$Y)[2]; p = dim(data$Y)[3];  
load('sst_locations.Rdata')

## Get the static Bayesian functional graphical model ---------
folder_name = "Application_results"
file_name = paste("MCMC_output_spline_static.Rdata", sep = "")
# file_name = paste("MCMC_output_fourier_static.Rdata", sep = "")
load(file=paste(folder_name, '/', file_name, sep = ""))
str(mcmc_output_static)
ppi_local_static = apply(mcmc_output_static$adj_save, c(1,2), mean)

# plot 
ppi_local = ppi_local_static
adj_local = ppi_local > 0.7
#pheatmap(adj_local + 0, cluster_rows = F, cluster_cols = F)
#pheatmap(ppi_local, cluster_rows = F, cluster_cols = F)
adj_block = matrix(FALSE, p, p)
for (i in 1:(p-1)){
  for (j in (i+1):p){
    adj_block[i,j] = sum(adj_local[ ((i-1)*K+1):(i*K),  ((j-1)*K+1):(j*K) ]) >0
  }
}
diag(adj_block) = TRUE
adj_block = adj_block | t(adj_block)
# pheatmap(adj_block + 0, cluster_rows = F, cluster_cols = F)
# Plot blocked graph
par(mar=c(0,0,0,0))
map('world',
    col="#f2f2f2", fill=TRUE, bg="white", lwd=0.05,
    border=0,ylim=c(-20,20) , xlim = c(50, 200)
)
points(x=locations[,1], y=locations[,2], col="slateblue", cex=2, pch=20)
for (j_1 in 1:(p-1)){
  for (j_2 in (j_1+1):p){
    if (adj_block[j_1,j_2]){
      inter <- gcIntermediate(locations[j_1,], locations[j_2,] , addStartEnd=TRUE, breakAtDateLine=F)
      lines(inter, col="slateblue", lwd=1)
    }
  }
}

(sum(adj_block)-p)/2
# Fitted curves
B_est =  apply(mcmc_output_static$B_save, c(1,2), mean)
plot(data$Y[i,,j])
lines(mcmc_output_static$FLC%*% B_est[i,((j-1)*K + 1): (j*K)], type = 'l', col = 'red')
# Convergence
plot(mcmc_output_static$sigma_epsilon_save)
plot(mcmc_output_static$B_save[10,20,])


## Get the proposed DBFGM model ------------------------------------------------
folder_name = "Application_results"
file_name = paste("MCMC_output_spline_DBFGM.Rdata", sep = "")
# file_name = paste("MCMC_output_fourier_DBFGM.Rdata", sep = "")
load(paste(folder_name, '/', file_name, sep = ""))
str(mcmc_output_DBFGM)
# graphs
p = 15
K = 5
ppi_local_all = array(NA, dim = c(2,75,75))
for (s_i in 1:2){
  ppi_local_all[s_i,,] = apply(mcmc_output_DBFGM[[s_i]]$adj_save, c(1,2), mean)
}
for (s_i in 1:2){
  ppi_local = ppi_local_all[s_i,,]
  adj_local = ppi_local > 0.7
  #pheatmap(adj_local + 0, cluster_rows = F, cluster_cols = F)
  #pheatmap(ppi_local, cluster_rows = F, cluster_cols = F)
  adj_block = matrix(FALSE, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      adj_block[i,j] = sum(adj_local[ ((i-1)*K+1):(i*K),  ((j-1)*K+1):(j*K) ]) >0
    }
  }
  diag(adj_block) = TRUE
  adj_block = adj_block | t(adj_block)
  # pheatmap(adj_block + 0, cluster_rows = F, cluster_cols = F)
  # Plot blocked graph
  par(mar=c(0,0,0,0))
  map('world',
      col="#f2f2f2", fill=TRUE, bg="white", lwd=0.05,
      border=0,ylim=c(-20,20) , xlim = c(50, 200)
  )
  points(x=locations[,1], y=locations[,2], col="slateblue", cex=2, pch=20)
  for (j_1 in 1:(p-1)){
    for (j_2 in (j_1+1):p){
      if (adj_block[j_1,j_2]){
        inter <- gcIntermediate(locations[j_1,], locations[j_2,] , addStartEnd=TRUE, breakAtDateLine=F)
        lines(inter, col="slateblue", lwd=1)
      }
    }
  }
}

## Change point
plot(mcmc_output_DBFGM$changepoint_save)
mean(mcmc_output_DBFGM$changepoint_save)
sd(mcmc_output_DBFGM$changepoint_save)
# Fitted curves
B_est = list(); for (s_i in 1:2){B_est[[s_i]] = apply(mcmc_output_DBFGM[[s_i]]$B_save, c(1,2), mean)}
interval_ind = matrix(c(1, 181,182,364), nrow = 2, ncol = 2, byrow = TRUE)
X_est = compute_X(B_est, mcmc_output_DBFGM$FLC,interval_ind, p = 15)
i=1;j=1
plot(data$Y[i,,j])
lines(X_est[i,,j], col = 'red')


## Get the output of adding fourseasons ------------------------------------------------
folder_name = "Application_results"
file_name = paste("MCMC_output_DBFGM_fourseasons.Rdata", sep = "")
file_name = paste("MCMC_output_DBFGM_fourseasons_spline.Rdata", sep = "")
load(paste(folder_name, '/', file_name, sep = ""))
ppi_local_all = array(NA, dim = c(4,75,75))
for (s_i in 1:4){
  ppi_local_all[s_i,,] = apply(mcmc_output_fourseasons[[s_i]]$adj_save, c(1,2), mean)
}
# plot graphs
for (s_i in 1:4){
  ppi_local = ppi_local_all[s_i,,]
  adj_local = ppi_local > 0.7
  #pheatmap(adj_local + 0, cluster_rows = F, cluster_cols = F)
  #pheatmap(ppi_local, cluster_rows = F, cluster_cols = F)
  adj_block = matrix(FALSE, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      adj_block[i,j] = sum(adj_local[ ((i-1)*K+1):(i*K),  ((j-1)*K+1):(j*K) ]) >0
    }
  }
  diag(adj_block) = TRUE
  adj_block = adj_block | t(adj_block)
  # pheatmap(adj_block + 0, cluster_rows = F, cluster_cols = F)
  # Plot blocked graph
  par(mar=c(0,0,0,0))
  map('world',
      col="#f2f2f2", fill=TRUE, bg="white", lwd=0.05,
      border=0,ylim=c(-20,20) , xlim = c(50, 200)
  )
  points(x=locations[,1], y=locations[,2], col="slateblue", cex=2, pch=20)
  for (j_1 in 1:(p-1)){
    for (j_2 in (j_1+1):p){
      if (adj_block[j_1,j_2]){
        inter <- gcIntermediate(locations[j_1,], locations[j_2,] , addStartEnd=TRUE, breakAtDateLine=F)
        lines(inter, col="slateblue", lwd=1)
      }
    }
  }
 
}





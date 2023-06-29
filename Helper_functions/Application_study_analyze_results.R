library(ggplot2)
library(maps)
library(pheatmap)
rm(list = ls()); 

setwd("./Helper_functions")
#source('simulation_functions.R')
source('Helper_functions_mcmc.R')
source('performance_functions.R')
source('Helper_functions_sst.R')
setwd("..")

### Load data --------------------------
load('sst_processed_update.Rdata')
n = dim(data$Y)[1]; T_data = dim(data$Y)[2]; p = dim(data$Y)[3];  
selected_points = c(-130,30,
                    -160, 30,
                    150, 30,
                    170, 30,
                    -120,10,
                    -150, 10,
                    160, 10,
                    179, 10,
                    -90, -10,
                    -110, -10,
                    -140, -10,
                    -170, -10,
                    -90, -30,
                    -110, -30,
                    -140, -30,
                    -170, -30)
locations = matrix(selected_points, ncol = 2, byrow = T)

## Get the static Bayesian functional graphical model ---------
folder_name = "Application_results"
file_name = paste("MCMC_output_spline_static_update.Rdata", sep = "")
# file_name = paste("MCMC_output_fourier_static.Rdata", sep = "")
load(file=paste(folder_name, '/', file_name, sep = ""))
ppi_local_static = apply(mcmc_output_static$adj_save, c(1,2), mean)

# Get estimated graph
ppi_local = ppi_local_static
adj_local = ppi_local > 0.5
#pheatmap(adj_local + 0, cluster_rows = F, cluster_cols = F)
#pheatmap(ppi_local, cluster_rows = F, cluster_cols = F)
adj_block = matrix(FALSE, p, p)
K = dim(ppi_local_static)[1]/p
for (i in 1:(p-1)){
  for (j in (i+1):p){
    adj_block[i,j] = sum(adj_local[ ((i-1)*K+1):(i*K),  ((j-1)*K+1):(j*K) ]) >0
  }
}
diag(adj_block) = TRUE
adj_block = adj_block | t(adj_block)
#pheatmap(adj_block + 0, cluster_rows = F, cluster_cols = F)

# Plot estimated graph
the_plot <- function()
{
  plot_graphs(locations, adj_block)
}
png(
  "sst_s1_update.png",
  width     = 5,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
par(
  mgp=c(6,2,0),
  mar      = c(10, 8, 3, 3),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2.5
)
the_plot()
dev.off()
# number of edges
sum(adj_block)-p
(sum(adj_block)-p)/(p*(p-1))

# Fitted curves
B_est =  apply(mcmc_output_static$B_save, c(1,2), mean)
plot(data$Y[i,,j])
lines(mcmc_output_static$FLC%*% B_est[i,((j-1)*K + 1): (j*K)], type = 'l', col = 'red')
# Convergence
plot(mcmc_output_static$sigma_epsilon_save)
plot(mcmc_output_static$B_save[10,20,])


## Get results of the proposed DBFGM model -------------------------------------
# Load MCMC result
folder_name = "Application_results"
file_name = paste("MCMC_output_spline_DBFGM_update.Rdata", sep = "")
# file_name = paste("MCMC_output_fourier_DBFGM.Rdata", sep = "")
load(paste(folder_name, '/', file_name, sep = ""))
## Get the estimated graphs and plot the graphs
K = dim(mcmc_output_DBFGM$FLC)[2]
ppi_local_all = array(NA, dim = c(2,p*K,p*K))
for (s_i in 1:2){
  ppi_local_all[s_i,,] = apply(mcmc_output_DBFGM[[s_i]]$adj_save, c(1,2), mean)
}
for (s_i in 1:2){
  ppi_local = ppi_local_all[s_i,,]
  adj_local = ppi_local > 0.5
  adj_block = matrix(FALSE, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      adj_block[i,j] = sum(adj_local[ ((i-1)*K+1):(i*K),  
                                      ((j-1)*K+1):(j*K) ]) >0
    }
  }
  diag(adj_block) = TRUE
  adj_block = adj_block | t(adj_block)
  # Plot blocked graph
  plt = plot_graphs(locations, adj_block)
  print(plt)
  print(sum(adj_block)-p)
  print((sum(adj_block)-p)/(p*(p-1)))
}

## Change point
plot(mcmc_output_DBFGM$changepoint_save)
mean(mcmc_output_DBFGM$changepoint_save)
sd(mcmc_output_DBFGM$changepoint_save)
## Fitted curves
B_est = list(); for (s_i in 1:2){B_est[[s_i]] = apply(mcmc_output_DBFGM[[s_i]]$B_save, c(1,2), mean)}
interval_ind = matrix(c(1, 181,182,364), nrow = 2, ncol = 2, byrow = TRUE)
X_est = compute_X(B_est, mcmc_output_DBFGM$FLC,interval_ind, p = 15)
i=1;j=1
plot(data$Y[i,,j])
lines(X_est[i,,j], col = 'red')
## Convergence
s_i = 1
plot(mcmc_output_DBFGM[[s_i]]$sigma_epsilon_save)
plot(mcmc_output_DBFGM[[s_i]]$B_save[1,2,])
plot(mcmc_output_DBFGM[[s_i]]$C_save[2,5,])


## Get the output of adding fourseasons ------------------------------------------------
folder_name = "Application_results"
file_name = paste("MCMC_output_DBFGM_fourseasons.Rdata", sep = "")
#file_name = paste("MCMC_output_DBFGM_fourseasons_spline.Rdata", sep = "")
load(paste(folder_name, '/', file_name, sep = ""))
K = dim(mcmc_output_DBFGM$FLC)[2]
ppi_local_all = array(NA, dim = c(4,p*K,p*K))
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
  
  plt = plot_graphs(locations, adj_block)
  print(plt)
 
}





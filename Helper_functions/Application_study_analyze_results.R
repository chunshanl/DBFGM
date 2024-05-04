library(ggplot2)
library(maps)
library(pheatmap)
rm(list = ls()); 

setwd("./Helper_functions")
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

#############
## Get the static Bayesian functional graphical model result ---------
#############
folder_name = "Application_results"
file_name = paste("MCMC_output_spline_static.Rdata", sep = "")
load(file=paste(folder_name, '/', file_name, sep = ""))
K = 5. # 10
# Get estimated graph
ppi_local_static = apply(mcmc_output_static$adj_save, c(1,2), mean)
ppi_local = ppi_local_static
adj_local = ppi_local > 0.5
adj_block = matrix(FALSE, p, p)
K = dim(ppi_local_static)[1]/p
for (i in 1:(p-1)){
  for (j in (i+1):p){
    adj_block[i,j] = sum(adj_local[ ((i-1)*K+1):(i*K),  ((j-1)*K+1):(j*K) ]) >0
  }
}
diag(adj_block) = TRUE
adj_block = adj_block | t(adj_block)
# Plot estimated graph
the_plot <- function()
{
  plot_graphs(locations, adj_block)
}
png(
  "sst_static.png",
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

# Convergence plot
plot(mcmc_output_static$sigma_epsilon_save)
plot(mcmc_output_static$B_save[10,20,])

#############
## Get results of the proposed DBFGM model -------------------------------------
#############
# Load MCMC result
folder_name = "Application_results"
file_name = paste("MCMC_output_DBFGM_10K.RData", sep = "")
#file_name = paste("MCMC_output_DBFGM_10K_a2b7.RData", sep = "")
load(paste(folder_name, '/', file_name, sep = ""))
K = dim(mcmc_output_DBFGM$FLC)[2]
num_interval = length(mcmc_output_DBFGM$C)
n_changepoint = num_interval - 1

## Get the estimated graphs and plot the graphs
ppi_local_est = array(NA, dim = c(2,p*K,p*K))
for (s_i in 1:2){
  ppi_local_est[s_i,,] = apply(mcmc_output_DBFGM[[s_i]]$adj_save, c(1,2), mean)
}
adj_block_est = list()
for (s_i in 1:2){
  ppi_local = ppi_local_est[s_i,,]
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
  # Print number of edges
  print(sum(adj_block)-p)
  print((sum(adj_block)-p)/(p*(p-1)))
  # Save
  adj_block_est[[s_i]] = adj_block
}

## Change point
changepoint_est = round(mean(mcmc_output_DBFGM$changepoint_save))
print(paste0('Estimated changepoint: ', changepoint_est))
plot(mcmc_output_DBFGM$changepoint_save)

B_est = list(); for (s_i in 1:2){B_est[[s_i]] = apply(mcmc_output_DBFGM[[s_i]]$B_save, c(1,2), mean)}
sigma_epsilon_est = c() 
for (s_i in 1:num_interval){sigma_epsilon_est[s_i] = mean(mcmc_output_DBFGM[[s_i]]$sigma_epsilon_save)}
print(paste0('sigma epsilon est = ', sigma_epsilon_est))

## Convergence
s_i = 1
plot(mcmc_output_DBFGM[[s_i]]$sigma_epsilon_save)
plot(mcmc_output_DBFGM[[s_i]]$B_save[1,2,])
plot(mcmc_output_DBFGM[[s_i]]$C_save[2,5,])

## Save the plots
adj_block = adj_block_all[[1]]
plot_graphs(locations, adj_block)
the_plot <- function()
{
  plot_graphs(locations, adj_block)
}
png(
  "sst_s1_K5.png",  # CHANGE
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

adj_block = adj_block_all[[2]]
plot_graphs(locations, adj_block)
the_plot <- function()
{
  plot_graphs(locations, adj_block)
}
png(
  "sst_s2_k5.png",
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

## Compute DIC
DIC = compute_dic(num_interval, changeopint_est, data, B_est, mcmc_output_DBFGM$FLC, sigma_epsilon_est, mcmc_output_DBFGM)
print(DIC)

##############
## Get the output of adding fourseasons ------------------------------------------------
##############
folder_name = "Application_results"
file_name = paste("MCMC_output_DBFGM_fourseasons.Rdata", sep = "")
load(paste(folder_name, '/', file_name, sep = ""))
K = dim(mcmc_output_DBFGM$FLC)[2]
num_interval = length(mcmc_output_fourseasons$C)
n_changepoint = num_interval - 1
# Get estimated graphs
ppi_local_all = array(NA, dim = c(num_interval,p*K,p*K))
for (s_i in 1:num_interval){
  ppi_local_all[s_i,,] = apply(mcmc_output_fourseasons[[s_i]]$adj_save, c(1,2), mean)
}
for (s_i in 1:num_interval){
  ppi_local = ppi_local_all[s_i,,]
  adj_local = ppi_local > 0.5
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

## Change point
changepoint_est = round(apply(mcmc_output_fourseasons$changepoint_save,1,mean))
print(paste0('Estimated changepoint: ', changepoint_est))

## Fitted curves
B_est = list(); for (s_i in 1:num_interval){B_est[[s_i]] = apply(mcmc_output_fourseasons[[s_i]]$B_save, c(1,2), mean)}
sigma_epsilon_est = c() 
for (s_i in 1:num_interval){sigma_epsilon_est[s_i] = mean(mcmc_output_fourseasons[[s_i]]$sigma_epsilon_save)}
print(paste0('sigma epsilon est = ', sigma_epsilon_est))

## Compute AIC
DIC = compute_dic(num_interval, changeopint_est, data, B_est, mcmc_output_fourseasons$FLC, sigma_epsilon_est, mcmc_output_fourseasons)
print(DIC)








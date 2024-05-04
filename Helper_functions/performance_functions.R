get_tp_fp_tn_fn = function( truth, result ){
  # Given the true adjacency matrix and the estimated edges
  # calculate true positive, false positive, true negative, false negative 
  # and MCC
  
  true_up = truth[upper.tri(truth)]
  res_up = result[upper.tri(result)]
  
  output = list()
  
  output$tp = sum(true_up & res_up)
  output$fp = sum(!true_up & res_up)
  output$tn = sum(!true_up & !res_up)
  output$fn = sum(true_up & !res_up)
  
  return(output)
}


get_mcmc_perf_graph = function(mcmc_output_oneint, data_oneint, 
                               K, p, standardize, block_thresh, disp = TRUE){
  # Compute performance metrics of a precision matrix (graph) estimation
  
  blocked_SSSL_performance = list()
  # Posterior mean
  ppi_local = apply(mcmc_output_oneint$adj_save, c(1,2), mean)
  adj_local = ppi_local > block_thresh
  Omega_b = apply(mcmc_output_oneint$C_save, c(1,2), mean)
  if (standardize){
    sdY = apply(data_oneint$Y,2,sd)
    Omega_b = diag(1/sdY) %*% Omega_b %*% diag(1/sdY)
  }
  Omega_b[!adj_local]=0
  
  # Blocked graph performance
  adj_block = matrix(FALSE, p, p)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      adj_block[i,j] = sum(adj_local[ ((i-1)*K+1):(i*K),  ((j-1)*K+1):(j*K) ]) >0
    }
  }
  diag(adj_block) = TRUE
  adj_block = adj_block | t(adj_block)
  adj_block_v1 = adj_block
  output = get_tp_fp_tn_fn(data_oneint$G_x_true, adj_block_v1)
  
  perf = list()
  perf$tpr_block = output$tp / (output$tp + output$fn) # tpr
  perf$fpr_block = output$fp / (output$fp + output$tn) # fpr
  perf$mcc_block = (output$tp * output$tn - output$fp * output$fn) /
    (sqrt(output$tp + output$fp) * sqrt(output$tp + output$fn) * sqrt(output$tn + output$fp) * sqrt(output$tn + output$fn))
  blocked_SSSL_performance$block_perf = perf
  
  # Save MCMC summaries
  blocked_SSSL_performance$ppi_local = ppi_local
  blocked_SSSL_performance$Omega_b = Omega_b
  blocked_SSSL_performance$adj_local = adj_local
  blocked_SSSL_performance$adj_block = adj_block
  
  return(blocked_SSSL_performance)
  
}


print_mcmc_results = function(performance_graph, data){
  cat('Graph Estimation: TPR =', round(performance_graph$block_perf$tpr_block,2),
      ', FPR =', round(performance_graph$block_perf$fpr_block,2), 
      ', MCC =', round(performance_graph$block_perf$mcc_block,2),'\n')
  # if (K == data$K_true){
  #   cat('Coefficient space Graph Estimation: TPR =', round(performance_graph$local_perf$tpr,2),
  #       ', FPR =', round(performance_graph$local_perf$fpr,2), 
  #       ', MCC =', round(performance_graph$local_perf$mcc,2))
  #}
  
}

changepoint_to_static = function(data_changepoint, s_i){
  
  data = list()
  data$G_x_true = data_changepoint$param_true$G_x_true[[s_i]]
  data$G_b_true = data_changepoint$param_true$G_b_true[[s_i]]
  data$Omega_b_true = data_changepoint$param_true$Omega_b_true[[s_i]]
  data$B_true = data_changepoint$param_true$B_true[[s_i]]
  cp = data_changepoint$param_true$changepoint_true
  T_data = dim(data_changepoint$Y)[2]
  if (s_i == 1){
    interval_ind = 1:(cp - 1)
  }else{
    interval_ind = cp:T_data
  }
  data$F_true = data_changepoint$param_true$F_true[interval_ind, ]
  data$Y = data_changepoint$Y[, interval_ind, ]
  data$sigma_epsilon_true = data_changepoint$param_true$sigma_epsilon_true
  data$pii_local = data_changepoint$param_true$pii_local
  
  return(data)
}


compute_log_p = function(Y, B, FLC, interval_ind, p, sigma_epsilon_vec){
  n = dim(B[[1]])[1]
  T_data = dim(FLC)[1]
  K = dim(FLC)[2]
  X = array(NA, c(n, T_data, p))
  log_p = array(NA, c(n, T_data, p))
  for (i in 1:p){
    for (s_i in 1:dim(interval_ind)[1]){
      b = B[[s_i]][, ((i-1)*K+1):(i*K)]  # B is n*pk, b is of size n * K
      temp = interval_ind[s_i,1]:interval_ind[s_i,2]
      X[, temp, i] = b %*% t(FLC[temp, ])
      log_p[, temp, i] = -(Y[, temp, i] - X[, temp, i])^2/sigma_epsilon_vec[s_i]^2/2 -
        log(sigma_epsilon_vec[s_i]) - 1/2*log(2*pi)
    }
  }
  return(log_p) # n by T by p
}

compute_dic = function(num_interval, changeopint_est, data, B_est, FLC, sigma_epsilon_est, mcmc_output_DBFGM){
  n = dim(data$Y)[1]; T_data = dim(data$Y)[2]; p = dim(data$Y)[3];  
  interval_ind = matrix(NA, num_interval, 2)
  interval_ind[1:length(changepoint_est), 2] = changepoint_est -1
  interval_ind[2:num_interval, 1] = changepoint_est
  interval_ind[1,1] = 1
  interval_ind[num_interval,2] = T_data
  log_p_theta_hat = compute_log_p(data$Y, B_est, FLC, interval_ind, p, sigma_epsilon_est)
  
  log_p_iter = c()
  nsave = mcmc_output_DBFGM$nsave
  B_est_iter = list()
  sigma_epsilon_est_iter = c()
  for (iter in 1:nsave){
    for (j in 1:num_interval){
      B_est_iter[[j]] = mcmc_output_DBFGM[[j]]$B_save[,,iter]
      sigma_epsilon_est_iter[j] =  mcmc_output_DBFGM[[j]]$sigma_epsilon_save[iter]
    }
    log_p_iter[iter] = sum(compute_log_p(data$Y, B_est_iter, FLC, interval_ind, p, sigma_epsilon_est_iter))
  }
  
  
  DIC = -2*sum(log_p_theta_hat) + 2 * var(log_p_iter)
    #2*num_interval*(p*K)^2 + num_interval * prod( dim(B_est[[1]])) + num_interval
  print(paste0('DIC = ', DIC))
}


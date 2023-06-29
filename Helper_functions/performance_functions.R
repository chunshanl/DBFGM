fb_loss = function( A_true, A_est){
  temp1 = sqrt(sum((A_est - A_true)^2))
  temp2 = sqrt(sum(A_true^2))
  loss = temp1 / temp2
  return(loss);
}

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


get_perf_graph = function(true_adj, ppi, true_C, est_C ){
  # Compute performance metrics of a precision matrix (graph) estimation
  # no block structure
  #'@input
  # todo: Now only p by p, no third dimension
  # true_adj: p x p x S array of true adjacency matrices 
  # ppi: p x p x S array of posterior edge inclusion probabilities
  # true_C: p x p x S array of true precision matrices 
  # est_C: p x p x S array of estimated precision matrices 
  
  #' @return a list containing
  # tpr: true positive rate of edge estimation
  # fpr: false positive rate of edge estimation
  # mcc: MCC, Matthews correlation coefficient of edge estimation
  # auc: AUC, area under the curve of edge estimation
  # roc_curve: first column is false positive rates and second column is true
  # positive rates of edge estimation
  # fl: Frobenius loss from comparing true and estimated precision matrix
  
  perf = list()
  #----------------------------------------------------------------------------
  # Compute tp, fp, tn, fn, mcc
  median_model = ppi > .5 
  output = get_tp_fp_tn_fn( true_adj, median_model )
  perf$tpr = output$tp / (output$tp + output$fn);
  perf$fpr = output$fp / (output$fp + output$tn);
  perf$mcc = (output$tp * output$tn - output$fp * output$fn) / 
    (sqrt(output$tp + output$fp) * sqrt(output$tp + output$fn) * 
       sqrt(output$tn + output$fp) * sqrt(output$tn + output$fn))
  
  #----------------------------------------------------------------------------
  # Compute ROC
  opts = 0:1000 / 1000
  tpr_roc = matrix(0, 1001, 1); fpr_roc = matrix(0, 1001, 1)
  
  for (i in 1:1001){
    cur_threshold = opts[i]
    cur_adj = ppi > cur_threshold
    output = get_tp_fp_tn_fn( true_adj, cur_adj)
    tpr_roc[i] = output$tp / (output$tp + output$fn)
    fpr_roc[i] = output$fp / (output$fp + output$tn)
  }
  # plot(c(fpr_roc), c(tpr_roc), ylim = c(0,1), xlim = c(0,1))
  
  # if (tpr_roc[1]!=1 | fpr_roc[1]!=1){
  #   fpr_roc = c(1, fpr_roc)
  #   tpr_roc = c(1, tpr_roc)
  #   warning('Warning: cannot reach full graph when threshold = 0')
  # }
  
  perf$auc = sum( (fpr_roc[1:1000] - fpr_roc[2:1001]) * 
                    (tpr_roc[2:1001] + tpr_roc[1:1000]) / 2)
  
  perf$roc_curve = cbind(fpr_roc, tpr_roc)
  
  #----------------------------------------------------------------------------
  # Compute Frobenius loss
  est_C[!median_model] = 0
  perf$fl = fb_loss(true_C, est_C)
  
  return(perf)
}


get_mcmc_perf_blocked_graph_v2 = function(mcmc_output_oneint, data_oneint, 
                                          K, p, standardize, block_thresh, disp = TRUE){
  blocked_SSSL_performance = list()
  
  ppi_local = apply(mcmc_output_oneint$adj_save, c(1,2), mean)
  adj_local = ppi_local > block_thresh
  # pheatmap(adj_local + 0, cluster_rows = F, cluster_cols = F)
  Omega_b = apply(mcmc_output_oneint$C_save, c(1,2), mean)
  if (standardize){
    sdY = apply(data_oneint$Y,2,sd)
    Omega_b = diag(1/sdY) %*% Omega_b %*% diag(1/sdY)
  }
  Omega_b[!adj_local]=0

  # Blocked graph performance - version 1
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
  blocked_SSSL_performance$block_perf_v1 = perf
  if (disp){
    pheatmap(adj_block_v1+0, cluster_rows = FALSE, cluster_cols = FALSE, main = 'Estimated block graph')
    pheatmap(data_oneint$G_x_true+0, cluster_rows = FALSE, cluster_cols = FALSE, main = 'True block graph' )
  }
  # Block graph performance - version 2
  adj_block_save = array(NA, dim = c(p, p, dim(mcmc_output_oneint$adj_save)[3]))
  for (iter in 1:(dim(mcmc_output_oneint$adj_save)[3])){
    adj_block = matrix(FALSE, nrow = p, ncol = p)
    adj_local = mcmc_output_oneint$adj_save[,,iter]
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        adj_block[i,j] = sum(adj_local[ ((i-1)*K+1):(i*K),  ((j-1)*K+1):(j*K) ]) >0
      }
    }
    diag(adj_block) = TRUE
    adj_block = adj_block | t(adj_block)
    adj_block_save[,,iter] = adj_block
  }
  ppi_block_v2 = apply(adj_block_save, c(1,2), mean)
  adj_block_v2 = ppi_block_v2 > 0.5
  # output = get_tp_fp_tn_fn(data_oneint$G_x_true, adj_block_v2)
  perf = list()
  perf = get_perf_graph(data_oneint$G_x_true, ppi_block_v2, NA, NA)
  blocked_SSSL_performance$block_perf_v2 = perf
  # if (disp){
  #   pheatmap(pii_block, cluster_rows = FALSE, cluster_cols = FALSE, main = 'Local ppi')
  #   pheatmap(adj_block + 0, cluster_rows = FALSE, cluster_cols = FALSE, main = 'Local ppi')
  # }
  # Block-wise edge inclusion probabilities
  pii_block = apply(mcmc_output_oneint$pii_block_save, c(1,2), mean)
  adj_block_v3 = pii_block > (2/(K^2))
  output = get_tp_fp_tn_fn(data_oneint$G_x_true, adj_block_v3)
  perf = list()
  perf$tpr_block = output$tp / (output$tp + output$fn) # tpr
  perf$fpr_block = output$fp / (output$fp + output$tn) # fpr
  perf$mcc_block = (output$tp * output$tn - output$fp * output$fn) /
    (sqrt(output$tp + output$fp) * sqrt(output$tp + output$fn) * sqrt(output$tn + output$fp) * sqrt(output$tn + output$fn))
  blocked_SSSL_performance$block_perf_v3 = perf
  if (disp){
    pheatmap(pii_block, cluster_rows=FALSE, cluster_cols=FALSE, main = 'Posterior mean of pi')
    pheatmap(adj_block_v3 + 0, cluster_rows=FALSE, cluster_cols=FALSE, main = 'Graph form posterior mean of pi')
  }
  # Local graph performance
  if (K == data_oneint$K_true){
    perf = list()
    perf = get_perf_graph(data_oneint$G_b_true, ppi_local, data_oneint$Omega_b_true, Omega_b)
    blocked_SSSL_performance$local_perf = perf
    if (disp){
      plot(perf$roc_curve[,1], perf$roc_curve[,2], main = 'ROC of local graph', xlab = '', ylab = '')
      pheatmap(ppi_local, cluster_rows = FALSE, cluster_cols = FALSE, main = 'Local ppi')
      pheatmap(adj_local+0, cluster_rows = FALSE, cluster_cols = FALSE, main = 'Estimated local graph')
      pheatmap(data_oneint$G_b_true+0, cluster_rows = FALSE, cluster_cols = FALSE, main = 'True local graph' )
      pheatmap(Omega_b, cluster_rows=FALSE, cluster_cols=FALSE, main = 'Estimated Omega')
      pheatmap(data_oneint$Omega_b_true, cluster_rows=FALSE, cluster_cols=FALSE, main = 'True Omega')
    }
  }
  
  # Save MCMC summaries
  blocked_SSSL_performance$ppi_local = ppi_local
  blocked_SSSL_performance$Omega_b = Omega_b
  blocked_SSSL_performance$adj_local = adj_local
  blocked_SSSL_performance$adj_block = adj_block
  blocked_SSSL_performance$pii_block = pii_block
  
  return(blocked_SSSL_performance)
  
}


get_mcmc_perf_blocked_graph = function(mcmc_output, data, 
                                       K, p, disp = TRUE){
  
  output = list()
  
  adj_block_save = mcmc_output$adj_block_save; adj_local_save = mcmc_output$adj_local_save
  Omega_b_save = mcmc_output$C_save; 
  adj_block_true = data$G_x_true; adj_local_true = data$G_b_true; Omega_b_true =data$Omega_b_true
  
  # Performance - Block graph
  ppi_block = apply(adj_block_save, c(1,2), mean)
  # ppi_block[lower.tri(ppi_block, diag = TRUE)] = 0
  adj_block = ppi_block > 0.5
  if (disp){
    pheatmap(adj_block_true + 1, cluster_rows=FALSE, cluster_cols=FALSE, main = 'True block')
    pheatmap(adj_block+1, cluster_rows=FALSE, cluster_cols=FALSE, main = 'Estimated block')
    pheatmap(ppi_block, cluster_rows=FALSE, cluster_cols=FALSE, main = 'Block ppi')
  }
  
  temp = get_tp_fp_tn_fn(adj_block_true, adj_block )
  output$adj_block$tpr = temp$tp / (temp$tp + temp$fn) # tpr
  output$adj_block$fpr = temp$fp / (temp$fp + temp$tn) # fpr
  output$adj_block$mcc = (temp$tp * temp$tn - temp$fp * temp$fn) / 
    (sqrt(temp$tp + temp$fp) * sqrt(temp$tp + temp$fn) * sqrt(temp$tn + temp$fp) * sqrt(temp$tn + temp$fn))
  
  # Performance - local graph
  ppi_local =  apply(adj_local_save, c(1,2), mean)
  # ppi_local[lower.tri(ppi_local, diag = TRUE)] = 0
  adj_local = ppi_local > 0.5
  if(disp){
    pheatmap(adj_local_true+1-1, cluster_rows=FALSE, cluster_cols=FALSE, main = 'True local')
    pheatmap(adj_local+1-1, cluster_rows=FALSE, cluster_cols=FALSE, main = 'Estimated local')
    pheatmap(ppi_local, cluster_rows=FALSE, cluster_cols=FALSE, main = 'ppi local')
  }

  # diag(adj_block) = TRUE
  adj_block_span =  matrix(as.logical(kronecker(adj_block, matrix(1, K, K))), 
                           p*K, p*K)
  adj_prod = adj_local & adj_block_span
  ppi_prod = ppi_local * adj_prod
  
  if (disp){
    pheatmap(adj_prod+1-1, cluster_rows=FALSE, cluster_cols=FALSE, main = 'Estimated prod')
    pheatmap(ppi_prod, cluster_rows=FALSE, cluster_cols=FALSE, main = 'ppi prod')
  }
  
  if (dim(adj_local_save)[1] == dim(adj_local_true)[1]){
    Omega_b_est = apply(Omega_b_save, c(1,2), mean)
    # sdY = apply(data$Y,2,sd)
    # Omega_b_est = diag(1/sdY) %*% Omega_b_est %*% diag(1/sdY)
    Omega_b_est[!adj_prod] = 0
    # sdB = apply(data$B_true,2,sd)
    Omega_b_est = diag(1/mcmc_output$SDsd) %*% Omega_b_est %*% diag(1/mcmc_output$SDsd)
    if (disp){
      pheatmap(Omega_b_est, cluster_rows=FALSE, cluster_cols=FALSE, main = 'Estimated Omega')
    }
    perf = get_perf_graph(adj_local_true, ppi_prod, Omega_b_true, Omega_b_est)
    output$local = perf
  }
  
  return(output)
}


compute_performance_coefficients = function(B_save, B_true, FLC, time_intex, disp){
  
  # Not finished!
  
  B_est = apply(B_save, c(1,2), mean)
  fb_loss(B_true, B_est)
  
  n = dim(B_save)[1]
  K = dim(FLC)[2]
  p = dim(B_save)[2]/K
  i = sample(1:n, 1)
  j = sample(1:p, 1)
  
  
  # fitted curve
  fitted = c(data$F_true[time_index,] %*% B_est[i,((j-1)*K + 1): (j*K)],
             data$F_true[data$jump:T_all,] %*% B_est_after[i,((j-1)*K + 1): (j*K)])
  plot(time_intex, fitted, type = 'l', col = 'red')
  # true curve
  truef = c(data$F_true[1:(data$jump-1),] %*% data$B_true_before[i,((j-1)*K + 1): (j*K)],
            data$F_true[data$jump:T_all,] %*% data$B_true_after[i,((j-1)*K + 1): (j*K)])
  lines(1:T_all, truef)
  # data
  points(data$Y[i, ,j])
  
  # values
  B_est_before[i,((j-1)*K + 1): (j*K)]
  data$B_true_before[i,((j-1)*K + 1): (j*K)]
  
  # trace plot
  plot(mcmc_output$B_save_before[10,1,])
  
  
}


continuous_adjust = function(FLC, jump, Sig_before, Sig_after,
                             B_before_noadjust, B_after_noadjust){
  
  K = dim(FLC)[2]
  p = dim(B_before_noadjust)[2]/K
  
  temp = FLC[jump,]
  A = kronecker(diag(p), t(temp))
  A = cbind(A, -A)
  
  Sig_all = matrix(0, p*K*2, p*K*2)
  Sig_all[1:(p*K), 1:(p*K)] = Sig_before
  Sig_all[(p*K+1):(p*K*2), (p*K+1):(p*K*2)] = Sig_after
  
  transform_matrix = Sig_all %*% t(A) %*% 
    solve( A %*% Sig_all %*% t(A)) %*% A
  
  temp = cbind(B_before_noadjust, B_after_noadjust) - cbind(B_before_noadjust, B_after_noadjust) %*% 
    t(transform_matrix)
  B_before = temp[,1:(p*K)]
  B_after = temp[,(p*K+1):(2*p*K)]
  # F_true[jump_true,]%*% B_true_before_adjusted
  
  adjusted_B = list()
  adjusted_B$B_before = B_before
  adjusted_B$B_after = B_after
  
  return(adjusted_B)
}


test_convergence = function(data, mcmc_output, disp){
  
  p = dim(data$G_x_true)[1]
  K_true = dim(data$G_b_true)[1]/p
  K = dim(mcmc_output$adj)[1]/p
  n_save = dim(mcmc_output$Sig_save)[3]
  output = list()

  # Total number of block edges
  adj_block_save = array(NA, dim = c(p, p, n_save))
  for (iter in 1:(dim(mcmc_output$adj_save)[3])){
    adj_block = matrix(FALSE, nrow = p, ncol = p)
    adj_local = mcmc_output$adj_save[,,iter]
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        adj_block[i,j] = sum(adj_local[ ((i-1)*K+1):(i*K),  ((j-1)*K+1):(j*K) ]) >0
      }
    }
    diag(adj_block) = TRUE
    adj_block = adj_block | t(adj_block)
    adj_block_save[,,iter] = adj_block
  }
  output$adj_block_save = adj_block_save
  num_edges = (colSums(adj_block_save,dims = 2) - p)/2
  z_score = geweke.diag(num_edges)$z
  p_value = 2*pnorm(abs(z_score), lower.tail = FALSE)
  output$num_edges_block_p_value = p_value
  if (disp){
    plot(num_edges, main = '# edges in the functional space across iterations'); abline(h = mean(num_edges), col = 'red')
  }
   
  # Block edges
  p_value_block_adj = matrix(NA, ncol = p, nrow = p)
  for (i in 1:(p)){
    for (j in i:(p)){
      # z_score = geweke.diag(mcmc_output_SSSL$adj_save[i,j,]+0)$z
      z_score = geweke.diag(output$adj_block_save[i,j,]+0)$z
      p_value_block_adj[i,j] = 2*pnorm(abs(z_score), lower.tail = FALSE)
    }
  }
  output$block_edge_p_value = p_value_block_adj
  if (disp){
    temp = p_value_block_adj
    pheatmap(temp, cluster_rows=FALSE, cluster_cols=FALSE, main = 'p-value of edges in the functional space')
  }
  temp = sum(p_value_block_adj > 0.05, na.rm = TRUE)/ (p*(p-1)/2)
  output$block_edge_conv_pct = temp
  # cat('percentage of edges in the functional space that converges:', temp)
  # plot(mcmc_output$adj_save[1,2,])
  
  # Total number of local edges
  num_edges = (colSums(mcmc_output$adj_save,dims = 2) - p*K)/2
  z_score = geweke.diag(num_edges)$z
  p_value = 2*pnorm(abs(z_score), lower.tail = FALSE)
  output$num_edges_local_p_value = p_value
  if (disp){
    plot(num_edges, main = '# edges in the coeff space across iterations')
    abline(h = mean(num_edges), col = 'red')
  }
  
  # Local edges
  p_value_adj = matrix(NA, ncol = p*K, nrow = p*K)
  for (i in 1:(p*K)){
    for (j in i:(p*K)){
      # z_score = geweke.diag(mcmc_output_SSSL$adj_save[i,j,]+0)$z
      z_score = geweke.diag(mcmc_output$adj_save[i,j,]+0)$z
      p_value_adj[i,j] = 2*pnorm(abs(z_score), lower.tail = FALSE)
    }
  }
  # output$local_edge_p_value = p_value_adj
  # if (disp){
  #   temp = p_value_adj
  #   temp[!data$G_b] = NA
  #   pheatmap(temp, cluster_rows=FALSE, cluster_cols=FALSE, main = 'p-value (true edges)')
  #   temp = p_value_adj
  #   temp[data$G_b] = NA
  #   pheatmap(temp, cluster_rows=FALSE, cluster_cols=FALSE, main = 'p-value (no edge)')
  # }
  temp = sum(p_value_adj > 0.05, na.rm = TRUE)/ (p*K*(p*K-1)/2)
  output$local_edge_conv_pct = temp
  # cat('percentage of edges in the coefficient space that converges:', temp)

  # block-wise edge inclusion prob 
  p_value_pi = matrix(NA, ncol = p, nrow = p)
  for (i in 1:p){
    for (j in i:p){
      z_score = geweke.diag(mcmc_output$pii_block_save[i,j,])$z
      p_value_pi[i,j] = 2*pnorm(abs(z_score), lower.tail = FALSE)
    }
  }
  output$pi_p_value = p_value_pi
  if (disp){
    temp = p_value_pi
    temp[!data$G_x] = NA
    pheatmap(temp, cluster_rows=FALSE, cluster_cols=FALSE, main = 'p-value of pi (true edges)')
    temp = p_value_pi
    temp[data$G_x] = NA
    pheatmap(temp, cluster_rows=FALSE, cluster_cols=FALSE, main = 'p-value of pi (no edge)')
  }
  temp = sum(p_value_pi > 0.1, na.rm = TRUE)/(p*(p-1)/2)
  output$pi_conv_pct = temp
  # cat('percentage of pi that converges:', temp)
  # plot(mcmc_output$pii_block_save[1,12,])
  #plot(mcmc_output$pii_block_save[1,3,])
  #plot(mcmc_output$adj_save[5,12,])
  
  # C
  p_value_C = matrix(NA, ncol = p*K, nrow = p*K)
  for (i in 1:(p*K)){
    for (j in i:(p*K)){
      # z_score = geweke.diag(mcmc_output_SSSL$C_save[i,j,])$z
      z_score = geweke.diag(mcmc_output$C_save[i,j,])$z
      p_value_C[i,j] = 2*pnorm(abs(z_score), lower.tail = FALSE)
    }
  }
  output$C_p_value = p_value_C
  # if (disp){
  #   temp = p_value_C
  #   temp[!data$G_x] = NA
  #   pheatmap(temp, cluster_rows=FALSE, cluster_cols=FALSE, main = 'p-value (true edges)')
  #   temp = p_value_C
  #   temp[data$G_x] = NA
  #   pheatmap(temp, cluster_rows=FALSE, cluster_cols=FALSE, main = 'p-value (no edge)')
  # }
  temp = sum(p_value_C > 0.05, na.rm = TRUE)/ (p*K*(p*K-1)/2)
  output$C_conv_pct = temp
  # cat('percentage of elements in Omega that converges:', temp)
  
  return(output)
  
}

print_results = function(performance_graph, data){
  cat('Graph Estimation Version 1: TPR =', round(performance_graph$block_perf_v1$tpr_block,2),
        ', FPR =', round(performance_graph$block_perf_v1$fpr_block,2), 
        ', MCC =', round(performance_graph$block_perf_v1$mcc_block,2),'\n')
  cat('Graph Estimation Version 2: TPR =', round(performance_graph$block_perf_v2$tpr,2),
        ', FPR =', round(performance_graph$block_perf_v2$fpr,2), 
        ', MCC =', round(performance_graph$block_perf_v2$mcc,2), '\n')
  cat('Graph Estimation Version 3: TPR =', round(performance_graph$block_perf_v3$tpr_block,2),
        ', FPR =', round(performance_graph$block_perf_v3$fpr_block,2), 
        ', MCC =', round(performance_graph$block_perf_v3$mcc_block,2), '\n')
  if (K == data$K_true){
    cat('Coefficient space Graph Estimation: TPR =', round(performance_graph$local_perf$tpr,2),
        ', FPR =', round(performance_graph$local_perf$fpr,2), 
        ', MCC =', round(performance_graph$local_perf$mcc,2))
  }

}




get_mcmc_perf_graph = function(mcmc_output_oneint, data_oneint, 
                                          K, p, standardize, block_thresh, disp = TRUE){
  blocked_SSSL_performance = list()
  
  ppi_local = apply(mcmc_output_oneint$adj_save, c(1,2), mean)
  adj_local = ppi_local > block_thresh
  # pheatmap(adj_local + 0, cluster_rows = F, cluster_cols = F)
  Omega_b = apply(mcmc_output_oneint$C_save, c(1,2), mean)
  if (standardize){
    sdY = apply(data_oneint$Y,2,sd)
    Omega_b = diag(1/sdY) %*% Omega_b %*% diag(1/sdY)
  }
  Omega_b[!adj_local]=0
  
  # Blocked graph performance - version 1
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

  # Local graph performance
  if (K == data_oneint$K_true){
    perf = list()
    perf = get_perf_graph(data_oneint$G_b_true, ppi_local, data_oneint$Omega_b_true, Omega_b)
    blocked_SSSL_performance$local_perf = perf
  }
  
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
  if (K == data$K_true){
    cat('Coefficient space Graph Estimation: TPR =', round(performance_graph$local_perf$tpr,2),
        ', FPR =', round(performance_graph$local_perf$fpr,2), 
        ', MCC =', round(performance_graph$local_perf$mcc,2))
  }
  
}


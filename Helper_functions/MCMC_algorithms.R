MCMC_DBFGM = function(nburn, nsave, 
                      Y, 
                      K, 
                      v0, v1, a_pi, b_pi,
                      FLC, 
                      changepoint_interval,
                      changepoint_vec,
                      B, Sig, C, adj, pii_block,
                      sigma_epsilon,
                      disp){
  
  # Compute the dimensions:
  n = dim(Y)[1];  T_data = dim(Y)[2]; p = dim(Y)[3];  
  p_all = p * K   # dimension for Omega
  # Construct intervals
  num_interval = length(changepoint_vec) + 1
  interval_ind = matrix(NA, num_interval, 2)
  interval_ind[1:length(changepoint_vec), 2] = changepoint_vec -1
  interval_ind[2:num_interval, 1] = changepoint_vec
  interval_ind[1,1] = 1
  interval_ind[num_interval,2] = T_data
  
  # MCMC values ----------------------------------------------------------------------------
  temp = compute_mcmc_values(FLC, Y, K)
  tFF = temp$tFF; tFy = temp$tFy
  temp = c()
  
  # Parameters and stuff used in the SSSL Hao Wang prior
  lambda = 1
  V0 = v0 * matrix(1, p_all, p_all);
  V1 = v1 * matrix(1, p_all, p_all);
  tau = list()
  for (s_i in 1:num_interval){ tau[[s_i]] = V0; tau[[s_i]][adj[[s_i]]] = v1}
  
  # get ind_upper_block and idx_upper
  ind_all = matrix(1:p_all^2, p_all, p_all)
  ind_upper_block = c()  # matrix index in and sorted by upper diagonal blocks
  idx_all = matrix(1:p^2, p, p)
  idx_upper = c()  # vector index of the upper diagonal in the p-by-p matrix
  for (i in 1:(p - 1)){
    for (j in (i+1): p){
      idx_upper = c(idx_upper, idx_all[i,j])
      rows = ((i - 1) * K + 1) : (i * K)
      cols = ((j - 1) * K + 1) : (j * K)
      ind_upper_block = c(ind_upper_block, c(ind_all[rows, cols]))
    }
  } 
  # get ind_noi_all
  ind_noi_all = compute_ind_noi(p_all)
  
  # Expand pii_block
  pii_block_expand = list()
  for (s_i in 1:num_interval){ pii_block_expand[[s_i]] =  kronecker(pii_block[[s_i]], matrix(1, K, K))}
  
  # Store the MCMC output----------------------------------------------------------------------------
  mcmc_output = list()
  for (s_i in 1:num_interval){
    mcmc_output[[s_i]] = list()
    mcmc_output[[s_i]]$C_save = array(NA, c(p_all, p_all, nsave)) 
    mcmc_output[[s_i]]$adj_save = array(NA, c(p_all, p_all, nsave)) 
    mcmc_output[[s_i]]$pii_block_save = array(NA, c(p, p, nsave))
    mcmc_output[[s_i]]$B_save = array(NA, c(n, p_all, nsave))
    mcmc_output[[s_i]]$sigma_epsilon_save = rep(NA, nsave)
  }
  mcmc_output$changepoint_save = matrix(NA, length(changepoint_vec), nsave)
  
  #----------------------------------------------------------------------------
  # Run the MCMC ----------------------------------------------------------------------------
  
  # Total number of MCMC simulations:
  nmc = nburn + nsave
  timer0 = proc.time()[3] # For timing the sampler
  for(iter in 1:nmc){
    
    # Display result
    if( iter %% 500 == 0 ){
      cat('iter =', iter, '\n')
      print(changepoint_vec)
    }  
    
    # Sample concentration matrix and adj ----------------------------------------------------------------------------
    for (s_i in 1:num_interval){
      Scov = t(B[[s_i]]) %*% B[[s_i]] 
      output = sample_C(Scov, C[[s_i]], Sig[[s_i]], adj[[s_i]], tau[[s_i]],
                        pii_block_expand[[s_i]], V0, V1, lambda, ind_noi_all, n)
      C[[s_i]] = output$C; Sig[[s_i]] = output$Sig
      adj[[s_i]] = output$adj; tau[[s_i]] = output$tau
      
      # Sample block pii -----------------------------------------------------
      output = sample_pii_block(pii_block[[s_i]], adj[[s_i]], ind_upper_block, idx_upper, K, p, a_pi, b_pi)
      pii_block[[s_i]] = output
      pii_block_expand[[s_i]] = kronecker(pii_block[[s_i]], matrix(1, K, K))
      #pheatmap(pii_block[[s_i]], cluster_rows=FALSE, cluster_cols=FALSE)
      #pheatmap(data$param_true$G_x_true[[s_i]]+0, cluster_rows=FALSE, cluster_cols=FALSE)
    }
    # pheatmap(adj[[4]]+0, cluster_rows=FALSE, cluster_cols=FALSE)
    # pheatmap(pii_block[[4]]+0, cluster_rows=FALSE, cluster_cols=FALSE)
    
    # Sample factors ----------------------------------------------------------------------------
    B = sample_B(tFF, tFy, FLC, interval_ind, sigma_epsilon, C, p_all, n, T_data)
    
    ## Sample change point -------------------------------------
    # compute (Y-X)^ 2
    kernel_sum_list = list()
    for (s_i in 1:num_interval){
      kernel = array(NA, c(n, T_data, p))
      for (i in 1:p){
        b = B[[s_i]][, ((i-1)*K+1):(i*K)]  # B is n*pK, b is of size n * K
        kernel[,, i] = (Y[,, i] - b %*% t(FLC))^2
      }
      kernel_sum_list[[s_i]] = rowSums(colSums(kernel))/sigma_epsilon[s_i]^2
    }
    # Compute likelihood of data
    for (point_i in 1:length(changepoint_vec)){
      # for each point in the changepoint range, save the kernel sum in kernal_sum_all_changepoints
      kernal_sum_all_changepoints = rep(NA, T_data)
      changepoint_range = changepoint_interval[point_i,1]:changepoint_interval[point_i,2]
      
      interval_ind_temp = interval_ind
      for (changepoint_temp in changepoint_range){
        interval_ind_temp[point_i, 2] = changepoint_temp-1
        interval_ind_temp[point_i+1, 1] = changepoint_temp
        kernal_sum_temp = 0
        for (s_i in 1:num_interval){
          kernal_sum_temp = kernal_sum_temp + sum(kernel_sum_list[[s_i]][interval_ind_temp[s_i,1]:interval_ind_temp[s_i,2]])
        }
        kernal_sum_all_changepoints[changepoint_temp] = kernal_sum_temp
      }
      
      w1 = - 0.5 * kernal_sum_all_changepoints
      w1_max = max(w1, na.rm = TRUE)
      w = exp(w1 - w1_max)
      w = w[changepoint_range]
      changepoint_vec[point_i] = sample(changepoint_range, 1, prob=w)
      interval_ind[point_i, 2] = changepoint_vec[point_i]-1
      interval_ind[point_i+1, 1] = changepoint_vec[point_i]
    }
    
    # Sample sigma epsilon -------------------------------
    X = compute_X(B, FLC, interval_ind, p)
    for (s_i in 1:num_interval){
      time_index = interval_ind[s_i,1]:interval_ind[s_i,2]
      sigma_epsilon[s_i] = sample_sigma_epsilon(Y[,time_index,], X[,time_index,])
    }
    
    # Store output----------------------------------------------------------------------------
    if (iter > nburn){
      for (s_i in 1:num_interval){
        mcmc_output[[s_i]]$C_save[,,iter - nburn] = C[[s_i]]
        mcmc_output[[s_i]]$adj_save[,,iter - nburn] = adj[[s_i]]
        mcmc_output[[s_i]]$pii_block_save[,,iter - nburn] = pii_block[[s_i]]
        mcmc_output[[s_i]]$B_save[,,iter-nburn] = B[[s_i]]
        mcmc_output[[s_i]]$sigma_epsilon_save[iter-nburn] = sigma_epsilon[s_i]
      }
      mcmc_output$changepoint_save[,iter - nburn] = changepoint_vec
    }
    
    
  }  # end iteration
  
  running_time = proc.time()[3] - timer0  # 
  print(paste('Total time: ', round(running_time/60) , 'minutes'))
  mcmc_output$running_time = running_time
  
  return (mcmc_output);
}


SSSL = function(nburn, nsave, Y, 
                v0, v1, pii,   # parameters 
                Sig, C, adj, # initialization
                disp_result){
  
  #' @param Y 
  #' @param v0
  #' @param v1
  #' @param pii
  #' @param lambda
  #' @param Sig  initial guess of Sigma
  
  
  # Compute the dimensions:
  p = ncol(Y);  n = nrow(Y)
  Scov = t(Y) %*% Y
  
  #----------------------------------------------------------------------------
  # Parameters and stuff used in the SSSL Hao Wang prior
  #----------------------------------------------------------------------------
  V0 = v0 * matrix(1, p, p);
  V1 = v1 * matrix(1, p, p);
  lambda = 1
  
  tau = V0;
  tau[adj] = v1;
  
  ind_noi_all = matrix(0, p, p-1);
  for(i in 1:p){
    if(i==1){
      ind_noi = 2:p
    }else if(i==p){
      ind_noi = 1:(p-1) 
    }else{
      ind_noi = c(1:(i-1), (i+1):p)
    }
    ind_noi_all[i,] = ind_noi
  }
  ind_noi_all = t(ind_noi_all)
  
  
  #----------------------------------------------------------------------------
  # Store the MCMC output
  #----------------------------------------------------------------------------
  C_save = array(NA, c(p, p, nsave)) 
  Sig_save = C_save; adj_save = C_save
  mcmc_output = list()
  
  #----------------------------------------------------------------------------
  # Run the MCMC
  #----------------------------------------------------------------------------
  
  # Total number of MCMC simulations:
  nmc = nburn + nsave
  
  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(iter in 1:nmc){
    
    # Display result
    if( iter %% 500 == 0 ){
      cat('iter =', iter, '\n')
    }  
    
    for (i in 1:p){
      
      #----------------------------------------------------------------------------
      # Sample the concentration matrix      
      ind_noi = ind_noi_all[,i]
      
      tau_temp = tau[ind_noi,i]
      
      Sig11 = Sig[ind_noi,ind_noi]; Sig12 = Sig[ind_noi,i]
      
      invC11 = Sig11 - Sig12 %*% t(Sig12)/Sig[i,i]
      
      Ci = (Scov[i,i] + lambda) * invC11 + diag(1/tau_temp)
      Ci = (Ci + t(Ci))/2
      
      Ci_chol = chol(Ci)
      beta = backsolve(Ci_chol, forwardsolve(t(Ci_chol), - Scov[ind_noi,i])
                       + rnorm(p-1))
      
      C[ind_noi,i] = beta
      C[i,ind_noi] = beta
      
      a_gam = 0.5 * n + 1
      b_gam = (Scov[i,i] + lambda) * 0.5
      gam = rgamma( n = 1, shape = a_gam, rate = b_gam)
      
      c = beta %*% invC11 %*% beta
      C[i,i] = gam + c
      
      # Below updating Covariance matrix according to one-column change of precision matrix
      invC11beta = invC11 %*% beta
      Sig[ind_noi,ind_noi] = invC11 + invC11beta %*% t(invC11beta)/gam
      Sig12 = -invC11beta/gam
      Sig[ind_noi,i] = Sig12
      Sig[i,ind_noi] = t(Sig12)
      Sig[i,i] = 1/gam
      
      #----------------------------------------------------------------------------
      # Bernoulli update of the adjacency matrix
      v0 = V0[ind_noi,i]
      v1 = V1[ind_noi,i]
      
      w1 = - 0.5 * log(v0) - 0.5 * beta^2 / v0 + log(1-pii)
      w2 = - 0.5 * log(v1) - 0.5 * beta^2/v1 + log(pii)
      w_max = apply( cbind(w1,w2), 1, max)
      w = exp( w2 - w_max) /
        rowSums( exp( cbind(w1,w2)- cbind(w_max, w_max) ) )
      
      z = runif(p-1) < w
      
      v = v0
      v[z] = v1[z]
      tau[ind_noi,i] = v
      tau[i,ind_noi] = v
      
      adj[ind_noi,i] = z
      adj[i,ind_noi] = z
      
    }
    
    #----------------------------------------------------------------------------
    # Store the MCMC output:
    if (iter > nburn){
      Sig_save[,,iter - nburn] = Sig
      C_save[,,iter - nburn] = C
      adj_save[,,iter - nburn] = adj
    }   
  }
  
  # Store the parameters and results
  mcmc_output$Sig_save = Sig_save
  mcmc_output$C_save = C_save
  mcmc_output$adj_save = adj_save
  
  mcmc_output$nburn = nburn
  mcmc_output$v0 = v0
  mcmc_output$v1 = v1
  mcmc_output$pii = pii
  
  running_time = round((proc.time()[3] - timer0)/60) 
  print(paste('Total time: ', running_time, 'minutes'))
  mcmc_output$running_time = running_time
  
  return (mcmc_output);
}



MCMC_from_DBFGM_prior = function(p,K, 
                                 v0, v1, 
                                 a_pi, b_pi,
                                 Sig, C, adj, pii_block,
                                 nburn, nsave,
                                 disp){
  
  #' @param K
  #' @param v0
  #' @param v1
  #' @param pii_local
  #' @param pii_block
  #' @param lambda
  #' @param Sig  
  
  # Compute the dimensions:
  p_all = p * K   # dimension for Omega in the coefficient space
  
  # MCMC values ----------------------------------------------------------------------------
  # Parameters and values used in the precision matrix prior
  lambda = 1
  V0 = v0 * matrix(1, p_all, p_all);
  V1 = v1 * matrix(1, p_all, p_all);
  tau = V0; tau[adj] = v1
  # get ind_upper_block and idx_upper
  ind_all = matrix(1:p_all^2, p_all, p_all)
  ind_upper_block = c()  # matrix index in and sorted by upper diagonal blocks
  idx_all = matrix(1:p^2, p, p)
  idx_upper = c()  # vector index of the upper diagonal in the p-by-p matrix
  for (i in 1:(p - 1)){
    for (j in (i+1): p){
      idx_upper = c(idx_upper, idx_all[i,j])
      rows = ((i - 1) * K + 1) : (i * K)
      cols = ((j - 1) * K + 1) : (j * K)
      ind_upper_block = c(ind_upper_block, c(ind_all[rows, cols]))
    }
  } 
  # get ind_noi_all
  ind_noi_all = compute_ind_noi(p_all)
  # expand pii_block
  pii_block_expand = kronecker(pii_block, matrix(1, K, K))
  
  # Store the MCMC output----------------------------------------------------------------------------
  C_save = array(NA, c(p_all, p_all, nsave)) 
  pii_block_save = array(NA, c(p, p, nsave))
  #Sig_save = C_save; 
  adj_save = C_save
  
  # Run the MCMC ----------------------------------------------------------------------------
  nmc = nburn + nsave  # total number of MCMC simulations:
  timer0 = proc.time()[3] # For timing the sampler
  for(iter in 1:nmc){
    
    # Display result
    if( disp & (iter %% 500 == 0) ){
      cat('iter =', iter, '\n')
    }  
    
    # Sample precision matrix and graph ----------------------------------------------------------------------------
    output = sample_C_from_prior(C, Sig, adj, 
                                   tau, pii_block_expand, 
                                   V0, V1, lambda, ind_noi_all)
    C = output$C; Sig = output$Sig
    adj = output$adj; tau = output$tau
    
    # Sample block-wise edge inclusion probabilities -----------------------------------------------------
    output = sample_pii_block(pii_block, adj, ind_upper_block, idx_upper, K, p, a_pi, b_pi)
    pii_block = output
    pii_block_expand = kronecker(pii_block, matrix(1, K, K))
    
    # Store output----------------------------------------------------------------------------
    if (iter > nburn){
      C_save[,,iter - nburn] = C
      adj_save[,,iter - nburn] = adj
      pii_block_save[,,iter - nburn] = pii_block
    }
  }  # end iteration
  
  running_time = proc.time()[3] - timer0  # 
  print(paste('Total time: ', round(running_time/60) , 'minutes'))
  
  mcmc_output = list()
  mcmc_output$C_save = C_save
  mcmc_output$adj_save = adj_save
  mcmc_output$pii_block_save = pii_block_save
  return (mcmc_output)
}


MCMC_DBFGM_static = function(nburn, nsave, 
                             Y, 
                             K, 
                             v0, v1, a_pi, b_pi,
                             FLC, 
                             B, Sig, C, adj, pii_block,
                             sigma_epsilon,
                             disp){
  
  # Compute the dimensions:
  n = dim(Y)[1];  T_data = dim(Y)[2]; p = dim(Y)[3];  
  p_all = p * K   # dimension for Omega
  
  # MCMC values ----------------------------------------------------------------------------
  X = array(NA, c(n, T_data, p))
  
  tFF.sum = matrix(0, p_all, p_all)
  for (t in 1:T_data){
    temp = kronecker(diag(p), t(FLC[t,]))
    tFF.sum = tFF.sum + t(temp) %*% temp
  }
  tFF.sum[abs(tFF.sum)< 1e-7] = 0
  
  tFy.sum = matrix(0, p_all, n)
  for (t in 1:T_data){
    temp1 = t(FLC[t,])
    temp2 = kronecker(diag(p), temp1)
    Yt = Y[,t,]
    tFy.sum = tFy.sum + t(Yt %*% temp2)
  }
  
  # Parameters and stuff used in the blocked-SSSL prior ----------------------
  lambda = 1
  
  V0 = v0 * matrix(1, p_all, p_all);
  V1 = v1 * matrix(1, p_all, p_all);
  tau = V0;
  tau[adj] = v1;
  
  ind_noi_all = compute_ind_noi(p_all)
  
  ind_all = matrix(1:p_all^2, p_all,p_all)
  ind_upper_block = c()  # matrix index in and sorted by upper diagonal blocks
  idx_all = matrix(1:p^2, p, p)
  idx_upper = c()  # matrix index the upper diagonal
  for (i in 1:(p - 1)){
    for (j in (i+1): p){
      idx_upper = c(idx_upper, idx_all[i,j])
      
      rows = ((i - 1) * K + 1) : (i * K)
      cols = ((j - 1) * K + 1) : (j * K)
      ind_upper_block = c(ind_upper_block, c(ind_all[rows, cols]))
    }
  } 
  
  ## Expand pii_block
  pii_block_expand = kronecker(pii_block, matrix(1, K, K))
  
  # Store the MCMC output -------------------------------------
  C_save = array(NA, c(p_all, p_all, nsave)) 
  Sig_save = C_save; adj_save = C_save
  pii_block_save = array(NA, c(p, p, nsave))
  B_save = array(NA, c(n, p_all, nsave))
  sigma_epsilon_save = array(NA, nsave)
  mcmc_output = list()
  
  #----------------------------------------------------------------------------
  # Run the MCMC
  #----------------------------------------------------------------------------
  
  # Total number of MCMC simulations:
  nmc = nburn + nsave
  
  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(iter in 1:nmc){
    
    # Display result
    if( (iter %% 500 == 0) & disp ){
      cat('iter =', iter, '\n')
    }  
    
    # ----------------------------------------------------------------------------
    # Sample concentration matrix and the adjacency matrix ---------------------------------
    Scov = t(B) %*% B
    output = sample_C(Scov, C, Sig,adj, tau, 
                      pii_block_expand, V0, V1, lambda, ind_noi_all, n)
    C = output$C; Sig = output$Sig 
    adj = output$adj; tau = output$tau
    
    # Sample block pii
    output = sample_pii_block(pii_block, adj, ind_upper_block, idx_upper, K, p, a_pi, b_pi)
    pii_block = output
    pii_block_expand = kronecker(pii_block, matrix(1, K, K))
    
    # -----------------------------------------------------------------------------
    # Sample factors --------------------------------------------------------------
    Q = 1 / sigma_epsilon^2 * tFF.sum + C
    l = 1 / sigma_epsilon^2 * tFy.sum
    temp = array(rnorm(p_all*n), c(p_all,n))
    Q_chol = chol(Q)
    B = backsolve(Q_chol, forwardsolve(t(Q_chol),l) +
                    temp)
    B = t(B)
    
    # --------------------------------------------------
    # Sample variance of errors
    for (i in 1:p){
      b = B[, ((i-1)*K+1):(i*K)]  # b is of size n * K
      X[,,i] = b %*% t(FLC)
    }
    gamma_shape = n * T_data * p / 2 + 1
    gamma_rate = sum((Y-X)^2)/2
    temp = rgamma(n = 1, shape = gamma_shape, rate = gamma_rate)
    sigma_epsilon = 1 / sqrt(temp)
    
    #----------------------------------------------------------------------------
    # Store the MCMC output:
    if (iter > nburn){
      Sig_save[,,iter - nburn] = Sig
      C_save[,,iter - nburn] = C
      adj_save[,,iter - nburn] = adj
      pii_block_save[,,iter - nburn] = pii_block
      B_save[,,iter-nburn] = B
      sigma_epsilon_save[iter-nburn] = sigma_epsilon
    }   
  }
  
  # Store the parameters and results
  mcmc_output$Sig_save = Sig_save
  mcmc_output$C_save = C_save
  mcmc_output$adj_save = adj_save
  mcmc_output$pii_block_save = pii_block_save
  mcmc_output$B_save = B_save
  mcmc_output$sigma_epsilon_save = sigma_epsilon_save
  
  running_time = round((proc.time()[3] - timer0)/60) 
  print(paste('Total time: ', running_time, 'minutes'))
  mcmc_output$running_time = running_time
  
  return (mcmc_output);
}


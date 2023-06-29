MCMC_changepoint_DBFGM = function(Y, 
                                  K, FLC, 
                                  v0, v1, a_pi, b_pi,
                                  changepoint_interval,
                                  changepoint,
                                  B, 
                                  Sig, C, adj, pii_block,
                                  sigma_epsilon,
                                  nburn, nsave, 
                                  disp){
                                 
  
  
  #' @param B 
  #' @param K
  #' @param v0
  #' @param v1
  #' @param pii_local
  #' @param pii_block
  #' @param lambda
  #' @param Sig  
  
  # Compute the dimensions:
  n = dim(Y)[1];  T_data = dim(Y)[2]; p = dim(Y)[3];  
  p_all = p * K   # dimension for Omega in the coefficient space
  interval_ind = matrix(NA, 2, 2)
  interval_ind[1,] = c(1, changepoint - 1)
  interval_ind[2,] = c(changepoint, T_data)
  
  # MCMC values ----------------------------------------------------------------------------
  temp = compute_mcmc_values(FLC, Y, K)
  tFF = temp$tFF; tFy = temp$tFy
  temp = c()
  
  # Parameters and values used in the precision matrix prior
  lambda = 1
  V0 = v0 * matrix(1, p_all, p_all);
  V1 = v1 * matrix(1, p_all, p_all);
  tau = list()
  for (s_i in 1:2){ tau[[s_i]] = V0; tau[[s_i]][adj[[s_i]]] = v1}
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
  pii_block_expand = list()
  for (s_i in 1:2){ pii_block_expand[[s_i]] = kronecker(pii_block[[s_i]], matrix(1, K, K))}
  
  # Store the MCMC output----------------------------------------------------------------------------
  C_save = list(); pii_block_save = list(); B_save = list(); sigma_epsilon_save = list()
  for (s_i in 1:2){
    C_save[[s_i]] = array(NA, c(p_all, p_all, nsave)) 
    pii_block_save[[s_i]] = array(NA, c(p, p, nsave))
    B_save[[s_i]] = array(NA, c(n, p_all, nsave))
    sigma_epsilon_save[[s_i]] = rep(NA, nsave)
  }
  #Sig_save = C_save; 
  adj_save = C_save
  changepoint_save = rep(NA, nsave)
  
  # Run the MCMC ----------------------------------------------------------------------------
  nmc = nburn + nsave  # total number of MCMC simulations:
  timer0 = proc.time()[3] # For timing the sampler
  for(iter in 1:nmc){
    
    # Display result
    if( disp & (iter %% 500 == 0) ){
      cat('iter =', iter, '\n')
      print(changepoint)
    }  
    
    # Sample precision matrices and graphs ----------------------------------------------------------------------------
    for (s_i in 1:2){
      Scov = t(B[[s_i]]) %*% B[[s_i]] 
      output = sample_C(Scov, C[[s_i]], Sig[[s_i]], adj[[s_i]], tau[[s_i]],
                        pii_block_expand[[s_i]], V0, V1, lambda, ind_noi_all, n)
      C[[s_i]] = output$C; Sig[[s_i]] = output$Sig
      adj[[s_i]] = output$adj; tau[[s_i]] = output$tau
      
      # Sample block-wise edge inclusion probabilities -----------------------------------------------------
      output = sample_pii_block(pii_block[[s_i]], adj[[s_i]], ind_upper_block, idx_upper, K, p, a_pi, b_pi)
      pii_block[[s_i]] = output
      pii_block_expand[[s_i]] = kronecker(pii_block[[s_i]], matrix(1, K, K))

    }

    
    # Sample factors ----------------------------------------------------------------------------
    B = sample_B(tFF, tFy, FLC, interval_ind, sigma_epsilon, C, p_all, n, T_data)
    
    ## Sample change point -------------------------------------
    # compute (Y-X)^ 2
    kernel_list = list()
    kernel = array(NA, c(n, T_data, p))
    for (s_i in 1:2){
      for (i in 1:p){
        b = B[[s_i]][, ((i-1)*K+1):(i*K)]  # B is n*pK, b is of size n * K
        kernel[,, i] = (Y[,, i] - b %*% t(FLC))^2
      }
      kernel_list[[s_i]] = kernel
    }

    kernel_1 = rowSums(colSums(kernel_list[[1]]))/sigma_epsilon[1]^2
    kernel_2 = rowSums(colSums(kernel_list[[2]]))/sigma_epsilon[2]^2
    kernal_sum_all_changepoints = rep(NA, T_data)
    for (changepoint_temp in changepoint_interval[1]:changepoint_interval[2]){
      kernal_sum_all_changepoints[changepoint_temp] =
        sum(kernel_1[1:(changepoint_temp-1)]) +
        sum(kernel_2[changepoint_temp:T_data])
    }
    kernal_sum_all_changepoints = kernal_sum_all_changepoints[changepoint_interval[1]:changepoint_interval[2]]
    
    w1 = - 0.5 * kernal_sum_all_changepoints
    w1_max = max(w1, na.rm = TRUE)
    w = exp(w1 - w1_max)
    changepoint = sample(changepoint_interval[1]:changepoint_interval[2], 1, prob=w)
    interval_ind = matrix(NA, 2, 2)
    interval_ind[1,] = c(1, changepoint - 1)
    interval_ind[2,] = c(changepoint, T_data)
    
    # Sample sigma epsilon -------------------------------
    X = compute_X(B, FLC, interval_ind, p)
    sigma_epsilon[1] = sample_sigma_epsilon(Y[,1:(changepoint - 1),], X[,1:(changepoint - 1),])
    sigma_epsilon[2] = sample_sigma_epsilon(Y[,changepoint:T_data,], X[,changepoint:T_data,])

    
    # Store output----------------------------------------------------------------------------
    if (iter > nburn){
      for (s_i in 1:2){
        #Sig_save[[s_i]][,,iter - nburn] = Sig[[s_i]]
        C_save[[s_i]][,,iter - nburn] = C[[s_i]]
        adj_save[[s_i]][,,iter - nburn] = adj[[s_i]]
        pii_block_save[[s_i]][,,iter - nburn] = pii_block[[s_i]]
        B_save[[s_i]][,,iter-nburn] = B[[s_i]]
        sigma_epsilon_save[[s_i]][iter-nburn] = sigma_epsilon[s_i]
      }

      changepoint_save[iter - nburn] = changepoint
      
    }
  }  # end iteration
  
  # Store the parameters and results
  mcmc_output = list()
  for (jj in 1:2){
    mcmc_output[[jj]] = list()
    mcmc_output[[jj]]$C_save = C_save[[jj]]
    #mcmc_output[[jj]]$Sig_save = Sig_save[[jj]]
    mcmc_output[[jj]]$adj_save = adj_save[[jj]]
    mcmc_output[[jj]]$pii_block_save = pii_block_save[[jj]]
    mcmc_output[[jj]]$B_save = B_save[[jj]]
    mcmc_output[[jj]]$sigma_epsilon_save = sigma_epsilon_save[[jj]]
  }
  mcmc_output$changepoint_save = changepoint_save
  
  running_time = proc.time()[3] - timer0  # 
  print(paste('Total time: ', round(running_time/60) , 'minutes'))
  mcmc_output$running_time = running_time
  
  return (mcmc_output);
}



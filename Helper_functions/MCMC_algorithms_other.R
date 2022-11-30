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
  # nmc = nburn+(nskip+1)*(nsave)
  nmc = nburn + nsave
  # skipcount = 0; isave = 0 # For counting
  
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


blocked_SSSL_v2 = function(burnin, nmc, Y, 
                           K,
                           v0, v1,         # parameters 
                           pii_block,
                           Sig, C, adj,   # initializations
                           disp_result){
  
  #' @param Y 
  #' @param v0
  #' @param v1
  #' @param pii
  #' @param lambda
  #' @param Sig  initial guess of Sigma
  
  
  # Compute the dimensions:
  p = ncol(Y);  n = nrow(Y);  nblock = p/K
  Scov = t(Y) %*% Y
  
  #----------------------------------------------------------------------------
  # Parameters and stuff used in the SSSL Hao Wang prior
  #----------------------------------------------------------------------------
  V0 = v0 * matrix(1, p, p);
  V1 = v1 * matrix(1, p, p);
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
  
  ind_all = matrix(1:p^2, p,p)
  ind_upper_block = c()  # matrix index in and sorted by upper diagonal blocks
  idx_all = matrix(1:nblock^2, nblock, nblock)
  idx_upper = c()  # matrix index the upper diagonal
  for (i in 1:(nblock - 1)){
    for (j in (i+1): nblock){
      idx_upper = c(idx_upper, idx_all[i,j])
      
      rows = ((i - 1) * K + 1) : (i * K)
      cols = ((j - 1) * K + 1) : (j * K)
      ind_upper_block = c(ind_upper_block, c(ind_all[rows, cols]))
    }
  } 
  
  # ind_diag_block = c()  # matrix index in and sorted by upper on-diagonal blocks
  # temp = ind_all
  # temp[lower.tri(temp, diag = T)] = NA
  # for(i in 1:nblock){
  #   rows = ((i - 1) * K + 1) : (i * K)
  #   cols = rows
  #   ind_diag_block = c(ind_diag_block, c(temp[rows, cols]))
  # }
  # ind_diag_block = ind_diag_block[!is.na(ind_diag_block)]
  
  ## Expand pii_block
  pii_block_expand = kronecker(pii_block, matrix(1, K, K))
  
  
  #----------------------------------------------------------------------------
  # Store the MCMC output
  #----------------------------------------------------------------------------
  C_save = array(NA, c(p, p, nsave)) 
  Sig_save = C_save; adj_save = C_save
  pii_block_save = array(NA, c(nblock, nblock, nsave))
  mcmc_output = list()
  
  #----------------------------------------------------------------------------
  # Run the MCMC
  #----------------------------------------------------------------------------
  
  # Total number of MCMC simulations:
  # nmc = nburn+(nskip+1)*(nsave)
  nmc = nburn + nsave
  # skipcount = 0; isave = 0 # For counting
  
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
      v0 = V0[ind_noi, i]
      v1 = V1[ind_noi, i]
      pii = pii_block_expand[ind_noi, i]
      
      w1 = - 0.5 * log(v0) - 0.5 * beta^2 / v0 + log(1-pii)
      w2 = - 0.5 * log(v1) - 0.5 * beta^2 / v1 + log(pii)
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
    # Sample edge inclusion probabilities
    ## Sample off-diagonal block pii
    adj_upper_block = array(adj[ind_upper_block], c(K^2, (nblock-1)*nblock/2))
    temp = colSums(adj_upper_block)
    temp1 = temp + 1
    temp2 = K^2 - temp + 1
    pii_sample = rbeta(n = (nblock - 1) * nblock/2, shape1 =  temp1, shape2 = temp2)
    pii_block_temp = matrix(0, nrow = nblock, ncol = nblock)
    pii_block_temp[idx_upper] = pii_sample
    pii_block_temp = pii_block_temp + t(pii_block_temp)
    ## Sample on-diagonal block pii
    # Fixed value
    diag(pii_block_temp) = diag(pii_block)
    ## 
    pii_block = pii_block_temp
    pii_block_expand = kronecker(pii_block, matrix(1, K, K))
    
    
    #----------------------------------------------------------------------------
    # Store the MCMC output:
    if (iter > nburn){
      Sig_save[,,iter - nburn] = Sig
      C_save[,,iter - nburn] = C
      adj_save[,,iter - nburn] = adj
      pii_block_save[,,iter - nburn] = pii_block
    }   
  }
  
  # Store the parameters and results
  mcmc_output$Sig_save = Sig_save
  mcmc_output$C_save = C_save
  mcmc_output$adj_save = adj_save
  mcmc_output$pii_block_save = pii_block_save
  
  running_time = round((proc.time()[3] - timer0)/60) 
  print(paste('Total time: ', running_time, 'minutes'))
  mcmc_output$running_time = running_time
  
  return (mcmc_output);
}



MCMC_known_FLC_v2_random_sigma = function(nburn, nsave, 
                                          Y, 
                                          K, 
                                          v0, v1, 
                                          FLC, 
                                          B, Sig, C, adj, pii_block,
                                          sigma_epsilon,
                                          disp){
  
  #' @param Y
  #' @param v0
  #' @param v1
  #' @param pii
  #' @param lambda
  #' @param Sig  initial guess of Sigma
  
  
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
  
  ind_noi_all = matrix(0, p_all, p_all-1);
  for(i in 1:p_all){
    if(i==1){
      ind_noi = 2:p_all
    }else if(i==p_all){
      ind_noi = 1:(p_all-1) 
    }else{
      ind_noi = c(1:(i-1), (i+1):p_all)
    }
    ind_noi_all[i,] = ind_noi
  }
  ind_noi_all = t(ind_noi_all)
  
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
    
    # -----------------------------------------------------------------------------
    # Sample factors --------------------------------------------------------------
    # Q = 1 / sigma_epsilon^2 * tFF.sum /n^2 + C/n^2
    # l = 1 / sigma_epsilon^2 * tFy.sum /n^2
    Q = 1 / sigma_epsilon^2 * tFF.sum + C
    l = 1 / sigma_epsilon^2 * tFy.sum
    temp = array(rnorm(p_all*n), c(p_all,n))
    Q_chol = chol(Q)
    B = backsolve(Q_chol, forwardsolve(t(Q_chol),l) +
                    temp)
    B = t(B)
    # B[10,((j*K)+1): (j*(K+1))]
    # data$B_true[10,((j*K)+1): (j*(K+1))]
    
    
    # ----------------------------------------------------------------------------
    # Sample concentration matrix and the adjacency matrix ---------------------------------
    Scov = t(B) %*% B
    
    for (i in 1:p_all){
      
      # Sample the concentration matrix -------------------------------------------------------
      ind_noi = ind_noi_all[,i]
      
      tau_temp = tau[ind_noi,i]
      
      Sig11 = Sig[ind_noi,ind_noi]; Sig12 = Sig[ind_noi,i]
      
      invC11 = Sig11 - Sig12 %*% t(Sig12)/Sig[i,i]
      
      Ci = (Scov[i,i] + lambda) * invC11 + diag(1/tau_temp)
      Ci = (Ci + t(Ci))/2
      
      Ci_chol = chol(Ci)
      beta = backsolve(Ci_chol, forwardsolve(t(Ci_chol), - Scov[ind_noi,i])
                       + rnorm(p_all-1))
      
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
      v0 = V0[ind_noi, i]
      v1 = V1[ind_noi, i]
      pii = pii_block_expand[ind_noi, i]
      
      w1 = - 0.5 * log(v0) - 0.5 * beta^2 / v0 + log(1-pii)
      w2 = - 0.5 * log(v1) - 0.5 * beta^2 / v1 + log(pii)
      w_max = apply( cbind(w1,w2), 1, max)
      w = exp( w2 - w_max) /
        rowSums( exp( cbind(w1,w2)- cbind(w_max, w_max) ) )
      
      z = runif(p_all-1) < w
      
      v = v0
      v[z] = v1[z]
      tau[ind_noi,i] = v
      tau[i,ind_noi] = v
      
      adj[ind_noi,i] = z
      adj[i,ind_noi] = z
      
    }
    
    #----------------------------------------------------------------------------
    # Sample edge inclusion probabilities
    ## Sample off-diagonal block pii
    adj_upper_block = array(adj[ind_upper_block], c(K^2, (p-1)*p/2))
    temp = colSums(adj_upper_block)
    temp1 = temp + 1
    temp2 = K^2 - temp + 1
    pii_sample = rbeta(n = (p - 1) * p/2, shape1 =  temp1, shape2 = temp2)
    pii_block_temp = matrix(0, nrow = p, ncol = p)
    pii_block_temp[idx_upper] = pii_sample
    pii_block_temp = pii_block_temp + t(pii_block_temp)
    ## Sample on-diagonal block pii
    # Fixed value
    diag(pii_block_temp) = diag(pii_block)
    ## 
    pii_block = pii_block_temp
    pii_block_expand = kronecker(pii_block, matrix(1, K, K))
    
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


MCMC_known_FLC_v2_random_sigma_packed = function(nburn, nsave, 
                                          Y, 
                                          K, 
                                          v0, v1, a_pi, b_pi,
                                          FLC, 
                                          B, Sig, C, adj, pii_block,
                                          sigma_epsilon,
                                          disp){
  
  #' @param Y
  #' @param v0
  #' @param v1
  #' @param pii
  #' @param lambda
  #' @param Sig  initial guess of Sigma
  
  
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
    # sample precision matrix and graph
    #output = sample_C(B, C, Sig,adj, tau, 
    #                  pii_block_expand, V0, V1, lambda, ind_noi_all)
    Scov = t(B) %*% B
    # pheatmap(Scov, cluster_rows  = F, cluster_cols = F)
    # set.seed(iter) ####################################!!!!!!!!!!!!!!!
    output = sample_C(Scov, C, Sig,adj, tau, 
                      pii_block_expand, V0, V1, lambda, ind_noi_all, n)
    C = output$C; Sig = output$Sig 
    adj = output$adj; tau = output$tau
    # pheatmap(adj+0, cluster_rows = F, cluster_cols = F)
    
    # Sample block pii
    output = sample_pii_block(pii_block, adj, ind_upper_block, idx_upper, K, p, a_pi, b_pi)
    pii_block = output
    pii_block_expand = kronecker(pii_block, matrix(1, K, K))
    # pheatmap(pii_block, cluster_rows = F, cluster_cols = F)
    
    # -----------------------------------------------------------------------------
    # Sample factors --------------------------------------------------------------
    # set.seed(iter) ####################################!!!!!!!!!!!!!!!
    # Q = 1 / sigma_epsilon^2 * tFF.sum /n^2 + C/n^2
    # l = 1 / sigma_epsilon^2 * tFy.sum /n^2
    Q = 1 / sigma_epsilon^2 * tFF.sum + C
    l = 1 / sigma_epsilon^2 * tFy.sum
    temp = array(rnorm(p_all*n), c(p_all,n))
    Q_chol = chol(Q)
    B = backsolve(Q_chol, forwardsolve(t(Q_chol),l) +
                    temp)
    B = t(B)
    # B[10,((j*K)+1): (j*(K+1))]
    # data$B_true[10,((j*K)+1): (j*(K+1))]
    
    # --------------------------------------------------
    # Sample variance of errors
    for (i in 1:p){
      b = B[, ((i-1)*K+1):(i*K)]  # b is of size n * K
      X[,,i] = b %*% t(FLC)
    }
    gamma_shape = n * T_data * p / 2 + 1
    gamma_rate = sum((Y-X)^2)/2
    # set.seed(iter) ############3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

MCMC_DBFGM_fourseasons = function(nburn, nsave, 
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
  interval_ind = matrix(NA, 4, 2)
  interval_ind[1,] = c(1, changepoint_vec[1] - 1)
  interval_ind[2,] = c(changepoint_vec[1], changepoint_vec[2] - 1)
  interval_ind[3,] = c(changepoint_vec[2], changepoint_vec[3] - 1)
  interval_ind[4,] = c(changepoint_vec[3], T_data)
  
  # MCMC values ----------------------------------------------------------------------------
  temp = compute_mcmc_values(FLC, Y, K)
  tFF = temp$tFF; tFy = temp$tFy
  temp = c()
  
  # Parameters and stuff used in the SSSL Hao Wang prior
  lambda = 1
  V0 = v0 * matrix(1, p_all, p_all);
  V1 = v1 * matrix(1, p_all, p_all);
  tau = list()
  for (s_i in 1:4){ tau[[s_i]] = V0; tau[[s_i]][adj[[s_i]]] = v1}
  
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
  for (s_i in 1:4){ pii_block_expand[[s_i]] =  kronecker(pii_block[[s_i]], matrix(1, K, K))}
  
  # Store the MCMC output----------------------------------------------------------------------------
  C_save = list(); pii_block_save = list(); B_save = list(); sigma_epsilon_save = list()
  for (s_i in 1:4){
    C_save[[s_i]] = array(NA, c(p_all, p_all, nsave)) 
    pii_block_save[[s_i]] = array(NA, c(p, p, nsave))
    B_save[[s_i]] = array(NA, c(n, p_all, nsave))
    sigma_epsilon_save[[s_i]] = rep(NA, nsave)
  }
  # Sig_save = C_save; 
  adj_save = C_save
  
  changepoint_save = matrix(NA, 3, nsave)
  
  #----------------------------------------------------------------------------
  # Run the MCMC ----------------------------------------------------------------------------
  
  # Total number of MCMC simulations:
  nmc = nburn + nsave
  timer0 = proc.time()[3] # For timing the sampler
  for(iter in 1:nmc){
    
    # Display result
    if( iter %% 100 == 0 ){
      cat('iter =', iter, '\n')
      print(changepoint_vec)
    }  
    
    # Sample concentration matrix and adj ----------------------------------------------------------------------------
    for (s_i in 1:4){
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
    kernel_all = matrix(NA, 4, T_data)
    kernel = array(NA, c(n, T_data, p))
    for (s_i in 1:4){
      for (i in 1:p){
        b = B[[s_i]][, ((i-1)*K+1):(i*K)]  # B is n*pK, b is of size n * K
        kernel[,, i] = (Y[,, i] - b %*% t(FLC))^2
      }
      kernel_all[s_i, ] = rowSums(colSums(kernel))
    }
    # Compute likelihood of data
    for (s_i in 1:3){
      kernal_sum_all_changepoints = rep(NA, T_data)
      # log_sigma_product = rep(NA, T_data)
      changepoint_range = changepoint_interval[s_i,1]:changepoint_interval[s_i,2]
      interval_ind_temp = interval_ind
      for (changepoint_temp in changepoint_range){
        interval_ind_temp[s_i, 2] = changepoint_temp-1
        interval_ind_temp[s_i+1, 1] = changepoint_temp
        kernal_sum_all_changepoints[changepoint_temp] =
          sum(kernel_all[1,interval_ind_temp[1,1]:interval_ind_temp[1,2]])/ sigma_epsilon[1]^ 2 +
          sum(kernel_all[2,interval_ind_temp[2,1]:interval_ind_temp[2,2]])/ sigma_epsilon[2]^ 2 +
          sum(kernel_all[3,interval_ind_temp[3,1]:interval_ind_temp[3,2]])/ sigma_epsilon[3]^ 2 +
          sum(kernel_all[4,interval_ind_temp[4,1]:interval_ind_temp[4,2]])/ sigma_epsilon[4]^ 2 
      }
      # plot(kernal_sum_all_changepoints)
      # changepoint = which.min(kernal_sum_all_changepoints)
      
      w1 = - 0.5 * kernal_sum_all_changepoints  
      # plot(exp(w1))
      w1_max = max(w1, na.rm = TRUE)
      # plot( w1 - w1_max)
      # min(w1 - w1_max, na.rm = TRUE)
      # plot(exp(w1 - w1_max))
      # sum(exp(w1 - w1_max), na.rm = TRUE)
      # w = exp( w1 - w1_max) / sum( exp( w1- w1_max ) )
      # plot(exp(w))
      w = exp(w1 - w1_max)
      w = w[changepoint_range]
      changepoint_vec[s_i] = sample(changepoint_range, 1, prob=w)
      interval_ind[s_i, 2] = changepoint_vec[s_i]-1
      interval_ind[s_i+1, 1] = changepoint_vec[s_i]
    }
    
    
    # Sample sigma epsilon -------------------------------
    X = compute_X(B, FLC, interval_ind, p)
    for (s_i in 1:4){
      time_index = interval_ind[s_i,1]:interval_ind[s_i,2]
      sigma_epsilon[s_i] = sample_sigma_epsilon(Y[,time_index,], X[,time_index,])
    }
    
    # Store output----------------------------------------------------------------------------
    if (iter > nburn){
      for (s_i in 1:4){
        # Sig_save[[s_i]][,,iter - nburn] = Sig[[s_i]]
        C_save[[s_i]][,,iter - nburn] = C[[s_i]]
        adj_save[[s_i]][,,iter - nburn] = adj[[s_i]]
        pii_block_save[[s_i]][,,iter - nburn] = pii_block[[s_i]]
        B_save[[s_i]][,,iter-nburn] = B[[s_i]]
        sigma_epsilon_save[[s_i]][iter-nburn] = sigma_epsilon[s_i]
      }
      
      changepoint_save[,iter - nburn] = changepoint_vec
      
    }
  }  # end iteration
  
  # Store the parameters and results
  mcmc_output = list()
  for (jj in 1:4){
    mcmc_output[[jj]] = list()
    mcmc_output[[jj]]$C_save = C_save[[jj]]
    # mcmc_output[[jj]]$Sig_save = Sig_save[[jj]]
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





init_B = function(FLC, Y, n, p, K, changepoint_vec, T_data){
  B = list(); 
  n_interval = length(changepoint_vec) + 1
  for (s_i in 1:n_interval){
    B[[s_i]] = array(rnorm(n*p*K, mean = 0, sd = 1), c(n, p*K))
  }

  interval_ind = matrix(NA, n_interval, 2)
  interval_ind[1:length(changepoint_vec), 2] = changepoint_vec -1
  interval_ind[2:n_interval, 1] = changepoint_vec
  interval_ind[1,1] = 1
  interval_ind[n_interval,2] = T_data
  
  non_zero_basis = list()
  temp = 1:K
  for (s_i in 1:n_interval){
    temp1 = colSums(FLC[interval_ind[s_i,1]:interval_ind[s_i,2], ])
    non_zero_basis[[s_i]] = (temp[abs(temp1 )> 1])
  }
  for (i in 1:n){
    for (j in 1:p){
      for (s_i in 1:n_interval){
        time_index = interval_ind[s_i,1]:interval_ind[s_i,2]
        FLC_temp = FLC[time_index, non_zero_basis[[s_i]]]
        temp1 = ((j-1)*K+1):(j*K); temp2 = temp1[non_zero_basis[[s_i]]]
        B[[s_i]][ i, temp2 ] = solve( t(FLC_temp) %*% FLC_temp) %*% 
          t(FLC_temp) %*% Y[i, time_index, j]
      }
    }
  }
  return(B)
}

compute_X = function(B, FLC, interval_ind, p){
  n = dim(B[[1]])[1]
  T_data = dim(FLC)[1]
  K = dim(FLC)[2]
  X = array(NA, c(n, T_data, p))
  for (i in 1:p){
    for (s_i in 1:dim(interval_ind)[1]){
      b = B[[s_i]][, ((i-1)*K+1):(i*K)]  # B is n*pk, b is of size n * K
      temp = interval_ind[s_i,1]:interval_ind[s_i,2]
      X[, temp, i] = b %*% t(FLC[temp, ])
    }
  }
  return(X)
}


compute_mcmc_values = function(FLC, Y, K){
  
  n = dim(Y)[1];  T_data = dim(Y)[2]; p = dim(Y)[3];  
  p_all = p * K   # dimension for Omega
  
  tFF = array(NA, c(T_data, K, K))
  tFy = array(NA, c(T_data, p_all, n))
  for (t in 1:T_data){
    tFF[t,,] = FLC[t,] %*% t(FLC[t,])
    temp1 = t(FLC[t,])
    temp2 = kronecker(diag(p), temp1)
    Yt = Y[,t,]
    tFy[t,,] =t(Yt %*% temp2)
  }
  
  tFF[abs(tFF)< 1e-7] = 0
  
  output = list()
  output$tFF = tFF
  output$tFy = tFy
  return (output);
}

compute_ind_noi = function(p_all_input){
  ind_noi_all = matrix(0, p_all_input, p_all_input-1);
  for(i in 1:p_all_input){
    if(i==1){
      ind_noi = 2:p_all_input
    }else if(i==p_all_input){
      ind_noi = 1:(p_all_input-1) 
    }else{
      ind_noi = c(1:(i-1), (i+1):p_all_input)
    }
    ind_noi_all[i,] = ind_noi
  }
  ind_noi_all = t(ind_noi_all)
  return(ind_noi_all)
}

sample_C = function(Scov, C_input, Sig_input,adj_input, tau_input, 
                    pii_block_expand_input, V0, V1, lambda, ind_noi_all, n){
  
  p_all_input = dim(Scov)[1]
  
  
  for (i in 1:p_all_input){
    
    # Sample the concentration matrix -------------------------------------------------------
    ind_noi = ind_noi_all[,i]
    
    tau_temp = tau_input[ind_noi,i]
    
    Sig11 = Sig_input[ind_noi,ind_noi]; Sig12 = Sig_input[ind_noi,i]
    
    invC11 = Sig11 - Sig12 %*% t(Sig12)/Sig_input[i,i]
    
    Ci = (Scov[i,i] + lambda) * invC11 + diag(1/tau_temp)
    Ci = (Ci + t(Ci))/2
    
    Ci_chol = chol(Ci)
    beta = backsolve(Ci_chol, forwardsolve(t(Ci_chol), - Scov[ind_noi,i])
                     + rnorm(p_all_input-1))
    
    C_input[ind_noi,i] = beta
    C_input[i,ind_noi] = beta
    
    a_gam = 0.5 * n + 1
    b_gam = (Scov[i,i] + lambda) * 0.5
    gam = rgamma( n = 1, shape = a_gam, rate = b_gam)
    
    c = beta %*% invC11 %*% beta
    C_input[i,i] = gam + c
    
    # Below updating Covariance matrix according to one-column change of precision matrix
    invC11beta = invC11 %*% beta
    Sig_input[ind_noi,ind_noi] = invC11 + invC11beta %*% t(invC11beta)/gam
    Sig12 = -invC11beta/gam
    Sig_input[ind_noi,i] = Sig12
    Sig_input[i,ind_noi] = t(Sig12)
    Sig_input[i,i] = 1/gam
    
    #----------------------------------------------------------------------------
    # Bernoulli update of the adjacency matrix
    v0 = V0[ind_noi, i]
    v1 = V1[ind_noi, i]
    pii = pii_block_expand_input[ind_noi, i]
    
    w1 = - 0.5 * log(v0) - 0.5 * beta^2 / v0 + log(1-pii)
    w2 = - 0.5 * log(v1) - 0.5 * beta^2 / v1 + log(pii)
    w_max = apply( cbind(w1,w2), 1, max)
    w = exp( w2 - w_max) /
      rowSums( exp( cbind(w1,w2)- cbind(w_max, w_max) ) )
    
    z = runif(p_all_input-1) < w
    
    v = v0
    v[z] = v1[z]
    tau_input[ind_noi,i] = v
    tau_input[i,ind_noi] = v
    
    adj_input[ind_noi,i] = z
    adj_input[i,ind_noi] = z
    
  }
  
  output = list()
  output$C = C_input
  output$Sig = Sig_input
  output$adj = adj_input
  output$tau = tau_input
  return(output)
}


sample_pii_block = function(pii_block_input, adj_input, ind_upper_block, idx_upper, K, p, a_pi, b_pi){
  adj_upper_block = array(adj_input[ind_upper_block], c(K^2, (p-1)*p/2))
  temp = colSums(adj_upper_block)
  temp1 = temp + a_pi
  temp2 = K^2 - temp + b_pi
  pii_sample = rbeta(n = (p - 1) * p/2, shape1 =  temp1, shape2 = temp2)
  pii_block_temp = matrix(0, nrow = p, ncol = p)
  pii_block_temp[idx_upper] = pii_sample
  pii_block_temp = pii_block_temp + t(pii_block_temp)
  # On-diagonal blocks are set
  diag(pii_block_temp) = diag(pii_block_input)
  return(pii_block_temp)
}


sample_B = function(tFF, tFy, FLC, interval_ind, sigma_epsilon, C, p_all, n, T_data){
  
  B = list()
  K = dim(FLC)[2]; p = p_all/K
  
  for (s_i in 1:dim(interval_ind)[1]){
    time_index = interval_ind[s_i, 1]:interval_ind[s_i, 2]
    tFF.sum_temp = colSums(tFF[time_index,,], dim = 1)
    tFF.sum_temp = kronecker(diag(p), tFF.sum_temp)
    tFy.sum_temp = colSums(tFy[time_index,,], dim = 1)
    Q = 1 / sigma_epsilon[s_i]^2 * tFF.sum_temp + C[[s_i]]
    l = 1 / sigma_epsilon[s_i]^2 * tFy.sum_temp
    temp = array(rnorm(p_all*n), c(p_all,n))
    Q_chol = chol(Q)
    B_temp = backsolve(Q_chol, forwardsolve(t(Q_chol),l) +
                         temp)
    B[[s_i]] = t(B_temp)
  }
  return(B);
}

sample_sigma_epsilon = function(Y_input, X_input){
  diff_temp = c(Y_input - X_input)
  gamma_shape = length(diff_temp) / 2 + 1
  gamma_rate = sum(diff_temp^2)/2
  temp = rgamma(n = 1, shape = gamma_shape, rate = gamma_rate)
  sigma_epsilon_output = 1 / sqrt(temp)
  return(sigma_epsilon_output)
}


sample_C_from_prior = function(C_input, Sig_input, adj_input, 
                               tau_input, pii_block_expand_input, 
                               V0, V1, lambda, ind_noi_all){
  
  p_all_input = dim(C_input)[1]
  Scov =  matrix(0, nrow = p_all_input, ncol = p_all_input) 
  n = 0
  
  for (i in 1:p_all_input){
    
    # Sample the concentration matrix -------------------------------------------------------
    ind_noi = ind_noi_all[,i]
    
    tau_temp = tau_input[ind_noi,i]
    
    Sig11 = Sig_input[ind_noi,ind_noi]; Sig12 = Sig_input[ind_noi,i]
    
    invC11 = Sig11 - Sig12 %*% t(Sig12)/Sig_input[i,i]
    
    Ci = (Scov[i,i] + lambda) * invC11 + diag(1/tau_temp)
    Ci = (Ci + t(Ci))/2
    
    Ci_chol = chol(Ci)
    beta = backsolve(Ci_chol, forwardsolve(t(Ci_chol), - Scov[ind_noi,i])
                     + rnorm(p_all_input-1))
    
    C_input[ind_noi,i] = beta
    C_input[i,ind_noi] = beta
    
    a_gam = 0.5 * n + 1
    b_gam = (Scov[i,i] + lambda) * 0.5
    gam = rgamma( n = 1, shape = a_gam, rate = b_gam)
    
    c = beta %*% invC11 %*% beta
    C_input[i,i] = gam + c
    
    # Below updating Covariance matrix according to one-column change of precision matrix
    invC11beta = invC11 %*% beta
    Sig_input[ind_noi,ind_noi] = invC11 + invC11beta %*% t(invC11beta)/gam
    Sig12 = -invC11beta/gam
    Sig_input[ind_noi,i] = Sig12
    Sig_input[i,ind_noi] = t(Sig12)
    Sig_input[i,i] = 1/gam
    
    #----------------------------------------------------------------------------
    # Bernoulli update of the adjacency matrix
    v0 = V0[ind_noi, i]
    v1 = V1[ind_noi, i]
    pii = pii_block_expand_input[ind_noi, i]
    
    w1 = - 0.5 * log(v0) - 0.5 * beta^2 / v0 + log(1-pii)
    w2 = - 0.5 * log(v1) - 0.5 * beta^2 / v1 + log(pii)
    w_max = apply( cbind(w1,w2), 1, max)
    w = exp( w2 - w_max) /
      rowSums( exp( cbind(w1,w2)- cbind(w_max, w_max) ) )
    
    z = runif(p_all_input-1) < w
    
    v = v0
    v[z] = v1[z]
    tau_input[ind_noi,i] = v
    tau_input[i,ind_noi] = v
    
    adj_input[ind_noi,i] = z
    adj_input[i,ind_noi] = z
    
  }
  
  output = list()
  output$C = C_input
  output$Sig = Sig_input
  output$adj = adj_input
  output$tau = tau_input
  return(output)
}



call_SSSL = function(Y,
                     burnin, nmc, 
                     v0, v1, pii, lambda,
                     disp){
  #' @param Y 
  #' @param v0
  #' @param v1
  #' @param pii
  #' @param lambda
  #' @param disp
  #'
  #' @return results: list containing{ 
  #' @item C_save : p x p x nm
  #' c array of precision matrices across iterations
  #' @item adj_save : p x p x nmc array of adjacency matrices across iterations
  #' @item ppi_edges : p x p matrix of posterior probability of each edge
  #' }   

  
  # Compute the dimensions:
  p = ncol(Y);  n = nrow(Y)
  
  # Rescale data:
  meanY = colMeans(Y)
  sdY = apply(Y,2,sd)
  Y = t(Y) - meanY
  Y = t(Y/sdY)
  
  #----------------------------------------------------------------------------
  # Initialize the main terms:
  #----------------------------------------------------------------------------
  ## Sample values
  # Sig = cov(Y)
  # C = solve(Sig)
  ## Diag
  Sig = diag(p)
  C = Sig
  ## True values
  # C = data$Omega_b
  # Sig = solve(C)
  
  ## Based on C
  # adj = abs(C) > 1e-5  # initial adjacency matrices as defined by initial precision matrices 
  # C = C * adj
  # Dense
  # adj = matrix(T, nrow = p, ncol = p)
  # Diag
  adj = diag(T, p)
  # True values
  # adj = data$G_b

  #----------------------------------------------------------------------------
  # Run the MCMC algorithm
  #---------------------------------------------------------------------------
  mcmc_output = SSSL(burnin, nmc, Y, 
                         v0, v1, pii,   # parameters 
                         Sig, C, adj,   # initializations
                         disp)
  
  return(mcmc_output);
}


call_blocked_SSSL_v2 = function(Y, 
                                K,
                                burnin, nmc, 
                                v0, v1, #  lambda,
                                disp){
  
  #' @param Y 
  #' @param v0
  #' @param v1
  #' @param pii
  #' @param lambda
  #' @param disp
  
  #' @return results: list containing{ 
  #' @item C_save : p x p x nm
  #' c array of precision matrices across iterations
  #' @item adj_save : p x p x nmc array of adjacency matrices across iterations
  #' @item ppi_edges : p x p matrix of posterior probability of each edge
  #' }   
  
  # Compute the dimensions:
  p = ncol(Y);  n = nrow(Y);  nblock = p/K
  
  # Rescale data:
  # meanY = colMeans(Y)
  # sdY = apply(Y,2,sd)
  # Y = t(Y) - meanY
  # Y = t(Y/sdY)
  
  #----------------------------------------------------------------------------
  # Initialize the main terms:
  #----------------------------------------------------------------------------
  ## Sample values
  # Sig = cov(Y)
  # C = solve(Sig)
  ## Diag
  Sig = diag(p)
  C = Sig
  # True values
  # C = data$Omega_b
  # Sig = solve(C)
  
  ## Based on C
  # adj = abs(C) > 1e-5  # initial adjacency matrices as defined by initial precision matrices 
  # C = C * adj
  # Dense
  # adj = matrix(T, nrow = p, ncol = p)
  # Diag
  adj = diag(T, p)
  # True values
  # adj = data$G_b
  
  # Initialize pii_block (p-by-p)
  pii_block = matrix(1/2, nrow = nblock, ncol = nblock)   # block-wise edge inclusion probabilities of edges in the coefficient space
  
  
  #----------------------------------------------------------------------------
  # Run the MCMC algorithm
  #---------------------------------------------------------------------------
  mcmc_output = blocked_SSSL_v2(burnin, nmc, Y, 
                                K,
                                v0, v1,   # lambda      # parameters 
                                pii_block, 
                                Sig, C, adj,   # initializations
                                disp)
  mcmc_output$adj = adj
  mcmc_output$Sig = Sig
  mcmc_output$C = C
  #mcmc_output$SDmean = meanY
  #mcmc_output$SDsd = sdY
  
  mcmc_output$nburn = nburn
  mcmc_output$v0 = v0
  mcmc_output$v1 = v1
  mcmc_output$pii_block = pii_block
  
  return(mcmc_output);
}


call_DBFGM_static = function(data,
                             K,
                             nburn, nsave,
                             v0, v1, a_pi, b_pi,#   lambda,
                             basis_type,
                             disp){
  
  ## Inputs
  #'@input data: data$Y is the n x T x p array of observed time series
  #'@param K: number of basis functions
  #'@param nburn @param nsave: number of burn-in and saved MCMC iterations
  #'@param v0 @param v1 : spike and slab standard deviations in the precision matrix prior
  #'@param a_pi @param b_pi:  shape and rate in the gamma prior of edge inclusion probability  
  #'@param basis_type: spline, polynomial, fourier
  #'@param disp: true or false, whether to display the number of iterations
  
  ## Outputs: a list containing
  #'@item Sig_save: pK x pK x S x nmc array of covariance matrices across iterations
  #'@item C_save: pK x pK x S x nmc array of precision matrices across iterations
  #'@item adj_save: pK x pK x S x nmc array of adjacency matrices across iterations
  #'@item pii_block_save: p x p x S x nmc array of edge inclusion probability matrices across iterations
  #'@item B_save: n x pK x nmc array of basis coefficients across iterations
  #'@item sigma_epsilon_save: nmc array of standard deviation across iterations
  #'@item running_time: in minutes
  
  ## Initial Values:
  #'@param 
  
  
  
  Y = data$Y
  # Compute the dimensions:
  n = dim(Y)[1];  T_data = dim(Y)[2]; p = dim(Y)[3];  
  p_all = p * K   # dimension for Omega
  
  # Rescale data:
  # temp = Y
  # for (i in 1:n){
  #   for(j in 1:p){
  #     temp[i,,j] = (Y[i,,j] - mean(Y[i,,j]))/sd(Y[i,,j])
  #   }
  # }
  # Y = temp
  
  if (basis_type == 'True Basis'){
    FLC = data$F_true
  } else if (basis_type == 'spline'){
    U = seq(0, 1, length.out = T_data)
    knots = U[seq(0, length(U), length.out = K-1)]
    b = create.bspline.basis(rangeval = c(0,1), breaks = knots, norder = 4)
    FLC = eval.basis(U, b)   # T_data * K
  } else if ( basis_type == 'polynomial'){
    U = seq(0, 1, length.out = T_data)
    FLC = cbind(1/sqrt(T_data),
                poly(U, K - 1)) 
  } else if ( basis_type == 'fpca'){
    output = generate_fpca(Y, M = K)
    B = output$B
    FLC = output$FLC
    X = array(NA, c(n, T_data, p))
    for (i in 1:p){
      Y_hat = FLC %*% t(B[,((i-1)*K+1):(i*K)])
      X[,,i] = t(Y_hat)
    }
  } else if ( basis_type == 'fourier'){
    U = seq(0, 1, length.out = T_data)
    b <- create.fourier.basis(c(0,1),nbasis = K  )
    FLC = eval.basis(U, b)   # T_data * K
  }
  #matplot(FLC, type = 'l')
  
  #----------------------------------------------------------------------------
  # Initialize the main terms:
  #----------------------------------------------------------------------------
  # B = array(rnorm(n*p*K), c(n, p*K)) 
  # lm fit
  if (basis_type != 'fpca'){
    B = array(NA, c(n, p*K))
    X = array(NA, c(n, T_data, p))
    H_1 = solve(t(FLC)%*%FLC)%*%t(FLC)
    for (i in 1:p){
      B_hat = H_1 %*% t(Y[,,i])
      B[,((i-1)*K+1):(i*K)] = t(B_hat)
      Y_hat = FLC %*% B_hat
      X[,,i] = t(Y_hat)
    }
  }
  
  # B = data$B_true
  # i = 1; j = 1
  # plot(data$Y[i,,j])
  # lines(FLC%*% B[i,((j-1)*K + 1): (j*K)], type = 'l', col = 'red')
  # lines(data$F_true %*% data$B_true[i,((j-1)*K_true + 1): (j*K_true)], type = 'l')
  
  sigma_epsilon = sd(Y-X)
  # sigma_epsilon = data$sigma_epsilon_true
  
  
  Sig = diag(p_all); C = Sig
  # Sig = cov(B)
  # C = solve(Sig)
  # C = data$Omega_b_true; Sig = solve(C)
  
  adj = diag(TRUE, p_all) 
  
  pii_block = matrix(1/2, nrow = p, ncol = p)   # block-wise edge inclusion probabilities of edges in the coefficient space
  
  #----------------------------------------------------------------------------
  # Run the MCMC algorithm
  #---------------------------------------------------------------------------
  # mcmc_output = MCMC_known_FLC_v2(nburn, nsave,
  #                                 Y,
  #                                 K,
  #                                 v0, v1,
  #                                 FLC,
  #                                 B, Sig, C, adj, pii_block,
  #                                 sigma_epsilon = data$sigma_epsilon_true,
  #                                 disp)
  # mcmc_output = MCMC_known_FLC_v2_random_sigma(nburn, nsave,
  mcmc_output = MCMC_known_FLC_v2_random_sigma_packed(nburn, nsave,
                                               Y,
                                               K,
                                               v0, v1, a_pi, b_pi,
                                               FLC,
                                               B, Sig, C, adj, pii_block,
                                               sigma_epsilon,
                                               disp)
  
  
  mcmc_output$adj = adj
  mcmc_output$Sig = Sig
  mcmc_output$C = C
  # mcmc_output$meanY = meanY  # TODO
  # mcmc_output$sdY = sdY
  mcmc_output$nburn = nburn
  mcmc_output$v0 = v0
  mcmc_output$v1 = v1
  mcmc_output$pii_block = pii_block
  mcmc_output$B = B
  mcmc_output$FLC = FLC
  mcmc_output$sigma_epsilon = sigma_epsilon
  mcmc_output$a_pi = a_pi
  mcmc_output$b_pi = b_pi
  
  return(mcmc_output);
}

call_DBFGM_fourseasons = function(data,
                            K,
                            nburn, nsave,
                            v0, v1, a_pi, b_pi,
                            basis_type,
                            changepoint_interval,
                            disp){
  
  Y = data$Y
  # Compute the dimensions:
  n = dim(Y)[1]; T_data = dim(Y)[2]; p = dim(Y)[3];  
  p_all = p * K
  
  # Rescale data:
  # temp = Y
  # for (i in 1:n){
  #   for(j in 1:p){
  #     temp[i,,j] = (Y[i,,j] - mean(Y[i,,j]))/sd(Y[i,,j])
  #   }
  # }
  # Y = temp
  
  #----------------------------------------------------------------------------
  # Initialize the main terms:
  #----------------------------------------------------------------------------
  if (basis_type == 'True Basis'){
    FLC = data$param_true$F_true
  } else if (basis_type == 'spline'){
    U = seq(0, 1, length.out = T_data)
    knots = U[seq(0, length(U), length.out = K-1)]
    b = create.bspline.basis(rangeval = c(0,1), breaks = knots, norder = 4)
    FLC = eval.basis(U, b)   # T_data * K
  } else if ( basis_type == 'polynomial'){
    U = seq(0, 1, length.out = T_data)
    FLC = cbind(1/sqrt(T_data),
                poly(U, K - 1)) 
  } else if (basis_type == 'fpca'){
    FLC = generate_fpca(Y, M = K)
  } else if ( basis_type == 'fourier'){
    U = seq(0, 1, length.out = T_data)
    b <- create.fourier.basis(c(0,1),nbasis = K  )
    FLC = eval.basis(U, b)   # T_data * K
  }
  
  
  ### Initialize change point
  changepoint_vec = c(round(T_data/4),round(T_data/2), round(T_data/4*3))
  
  ### Initialize basis coefficients using LSE
  # Initialize B, the non-identifiable coefficients are set to LSE on the whole time interval
  B = list(); for (s_i in 1:4){B[[s_i]] = array(0, c(n, p*K))}
  for (i in 1:n){
    for (j in 1:p){
      temp = solve( t(FLC) %*% FLC) %*% t(FLC) %*% Y[i,, j]
      temp1 = ((j-1)*K+1):(j*K)
      for (s_i in 1:2){
        B[[s_i]][ i, temp1 ] = temp
      }
    }
  }
  # The identifiable coefficients are set to LSE
  interval_ind = matrix(NA, 4, 2)
  interval_ind[1,] = c(1, changepoint_vec[1] - 1)
  interval_ind[2,] = c(changepoint_vec[1], changepoint_vec[2] - 1)
  interval_ind[3,] = c(changepoint_vec[2], changepoint_vec[3] - 1)
  interval_ind[4,] = c(changepoint_vec[3], T_data)
  non_zero_basis = list()
  # for (s_i in 1:2){
  #   non_zero_basis[[s_i]] = 1:K
  # }
  temp = 1:K
  for (s_i in 1:4){
    temp1 = colSums(FLC[interval_ind[s_i, 1]:interval_ind[s_i, 2], ])
    # non_zero_basis[[s_i]] = (temp[temp1 != 0])
    non_zero_basis[[s_i]] = (temp[abs(temp1 )> 1])
  }
  for (i in 1:n){
    for (j in 1:p){
      for (s_i in 1:4){
        time_index = interval_ind[s_i, 1]:interval_ind[s_i, 2]
        FLC_temp = FLC[time_index, non_zero_basis[[s_i]]]
        temp1 = ((j-1)*K+1):(j*K); temp2 = temp1[non_zero_basis[[s_i]]]
        B[[s_i]][ i, temp2 ] = solve( t(FLC_temp) %*% FLC_temp) %*% 
          t(FLC_temp) %*% Y[i, time_index, j]
      }
    }
  }
  
  
  ### Initialize sigma_epsilon
  X = compute_X(B, FLC, interval_ind, p)
  sigma_epsilon = c()
  temp = Y - X
  for (s_i in 1:4){
    sigma_epsilon[s_i] = sd(temp[,interval_ind[s_i, 1]:interval_ind[s_i, 2],])
  }
  
  ### Initialize pii_block
  temp = matrix(1/2, nrow = p, ncol = p)   # block-wise edge inclusion probabilities of edges in the coefficient space
  pii_block = list()
  for (s_i in 1:4){pii_block[[s_i]] = temp}
  
  ### Initialize precision matrices and graphs
  # Use diag
  C = list(); adj = list()
  for (s_i in 1:4){C[[s_i]] = diag(p_all); adj[[s_i]] = diag(TRUE, p_all)}
  Sig = C
  #----------------------------------------------------------------------------
  # Run the MCMC algorithm
  #---------------------------------------------------------------------------
  #mcmc_output = MCMC_changepoint_newB(nburn, nsave,
  mcmc_output = MCMC_DBFGM_fourseasons(nburn, nsave, 
                                       Y, 
                                       K, 
                                       v0, v1, a_pi, b_pi,
                                       FLC, 
                                       changepoint_interval,
                                       changepoint_vec,
                                       B, Sig, C, adj, pii_block,
                                       sigma_epsilon,
                                       disp);
  
  # Save parameters and initializations
  mcmc_output$B = B
  mcmc_output$sigma_epsilon = sigma_epsilon
  mcmc_output$adj = adj; mcmc_output$Sig = Sig; mcmc_output$C = C
  mcmc_output$pii_block = pii_block
  mcmc_output$FLC = FLC
  mcmc_output$nburn = nburn; mcmc_output$nsave = nsave
  mcmc_output$v0 = v0; mcmc_output$v1 = v1
  # mcmc_output$meanY = meanY  # TODO
  # mcmc_output$sdY = sdY
  mcmc_output$a_pi = a_pi
  mcmc_output$b_pi = b_pi
  
  return(mcmc_output);
}




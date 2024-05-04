call_DBFGM = function(data_input,
                       K,
                       nburn, nsave,
                       v0, v1, a_pi, b_pi,
                       basis_type,
                       changepoint_interval,
                       disp){
  # Run MCMC for the proposed Dynamic Bayesian Function Graphical Model
  #' @param data_input : data_input$Y is of size n x T x p
  #' @param K : number of basis
  #' @basis_type : type of basis
  #' @param v0  : variance of the spike component in the prior of the precision matrix
  #' @param v1  : variance of the slab component
  #' @param a_pi, @param b_pi  : hyperparameters in the beta prior of block-wise edge inclusion probabilities
  #' @param changepoint_interval : a matrix of 2 columns that defines the range of change points
  #' @param nburn, @param nsave : number of burn-ins and number of saved iterations after burn-in
  #' @param disp : logical, whether to display the MCMC progress
  #'
  #' @return results: list containing{ 
  #' @item C_save : p x p x nmc array of precision matrices across iterations
  #' @item adj_save : p x p x nmc array of adjacency matrices across iterations
  #' @item ppi_edges : p x p matrix of posterior probability of each edge
  #' }   
  Y = data_input$Y
  # Compute the dimensions:
  n = dim(Y)[1]; T_data = dim(Y)[2]; p = dim(Y)[3];  
  p_all = p * K
  num_changepoints = dim(changepoint_interval)[1]
  num_interval = num_changepoints + 1
  
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
    FLC = data_input$param_true$F_true
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
  changepoint_vec = c()
  for (point_i in 1:num_changepoints){
    changepoint_vec = c(changepoint_vec, round(mean(changepoint_interval[point_i,])))
  }
  
  ### Initialize basis coefficients using LSE
  # Initialize B, the non-identifiable coefficients are set to LSE on the whole time interval
  B = init_B(FLC, Y, n, p, K, changepoint_vec, T_data)
  
  ### Initialize sigma_epsilon
  interval_ind = matrix(NA, num_interval, 2)
  interval_ind[1:length(changepoint_vec), 2] = changepoint_vec -1
  interval_ind[2:num_interval, 1] = changepoint_vec
  interval_ind[1,1] = 1
  interval_ind[num_interval,2] = T_data
  X = compute_X(B, FLC, interval_ind, p)
  sigma_epsilon = c()
  temp = Y - X
  for (s_i in 1:num_interval){
    sigma_epsilon[s_i] = sd(temp[,interval_ind[s_i, 1]:interval_ind[s_i, 2],])
  }
  
  ### Initialize pii_block
  temp = matrix(1/2, nrow = p, ncol = p)   # block-wise edge inclusion probabilities of edges in the coefficient space
  diag(temp) = 1
  pii_block = list()
  for (s_i in 1:num_interval){pii_block[[s_i]] = temp}
  
  ### Initialize precision matrices and graphs
  # Use diag
  C = list(); adj = list()
  for (s_i in 1:num_interval){C[[s_i]] = diag(p_all); adj[[s_i]] = diag(TRUE, p_all)}
  Sig = C
  
  #----------------------------------------------------------------------------
  # Run the MCMC algorithm
  #---------------------------------------------------------------------------
  mcmc_output = MCMC_DBFGM(nburn, nsave, 
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
  mcmc_output$changepoint_vec = changepoint_vec
  mcmc_output$sigma_epsilon = sigma_epsilon
  mcmc_output$adj = adj; mcmc_output$Sig = Sig; mcmc_output$C = C
  mcmc_output$pii_block = pii_block
  mcmc_output$FLC = FLC
  mcmc_output$nburn = nburn; mcmc_output$nsave = nsave
  mcmc_output$v0 = v0; mcmc_output$v1 = v1
  mcmc_output$a_pi = a_pi
  mcmc_output$b_pi = b_pi
  
  return(mcmc_output);
}


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
  Sig = diag(p)
  C = Sig
  adj = diag(T, p)

  #----------------------------------------------------------------------------
  # Run the MCMC algorithm
  #---------------------------------------------------------------------------
  mcmc_output = SSSL(burnin, nmc, Y, 
                     v0, v1, pii,   # parameters 
                     Sig, C, adj,   # initializations
                     disp)
  
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

  Y = data$Y
  # Compute the dimensions:
  n = dim(Y)[1];  T_data = dim(Y)[2]; p = dim(Y)[3];  
  p_all = p * K   # dimension for Omega
  
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
  
  #----------------------------------------------------------------------------
  # Initialize the main terms:
  #----------------------------------------------------------------------------
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
  
  sigma_epsilon = sd(Y-X)
  Sig = diag(p_all); C = Sig
  adj = diag(TRUE, p_all) 
  
  pii_block = matrix(1/2, nrow = p, ncol = p)   # block-wise edge inclusion probabilities of edges in the coefficient space
  
  #----------------------------------------------------------------------------
  # Run the MCMC algorithm
  #---------------------------------------------------------------------------
  mcmc_output = MCMC_DBFGM_static(nburn, nsave,
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


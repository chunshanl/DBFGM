call_DBFGM = function(data_input,
                      K, basis_type,
                      v0, v1, a_pi, b_pi,
                      changepoint_interval,
                      nburn, nsave,
                      disp){
  #' @param data_input : data_input$Y is of size n x T x p
  #' @param K : number of basis
  #' @basis_type : type of basis, spline, polynomial, fpca, or user define
  #' @param v0  : variance of the spike component in the prior of the precision matrix
  #' @param v1  : variance of the slab component
  #' @param a_pi, @param b_pi  : hyperparameters in the beta prior of block-wise edge inclusion probabilities
  #' @param changepoint_interval : a vector of length 2 that defines the range of change point
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
  
  # Rescale data_input:
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
  # matplot(FLC, type = 'l')
  
  ### Initialize change point
  changepoint = round(T_data/2)
  
  ### Initialize basis coefficients
  B = init_B(FLC, Y, n, p, K, changepoint, T_data)
  #B = data_input$param_true$B_true
  # X = compute_X(B, FLC, interval_ind, p)
  # plot(Y[10,,1])
  # lines(X[10,,1], col = 2)
  
  # s_i = 1
  # temp = var(B[[s_i]])
  # diag(temp)
  # plot(diag(temp))
  
  ### Initialize sigma_epsilon
  interval_ind = matrix(NA, 2, 2)
  interval_ind[1,] = c(1, changepoint - 1)
  interval_ind[2,] = c(changepoint, T_data)
  X = compute_X(B, FLC, interval_ind, p)
  sigma_epsilon = c()
  temp = Y - X
  for (s_i in 1:2){
    sigma_epsilon[s_i] = sd(temp[,interval_ind[s_i,1]:interval_ind[s_i,2],])
  }
  
  # sigma_epsilon = data_input$param_true$sigma_epsilon_true
  
  ### Initialize block-wise edge inclusion probabilities
  temp = matrix(1/2, nrow = p, ncol = p)   
  diag(temp) = 1
  pii_block = list()
  for (s_i in 1:2){pii_block[[s_i]] = temp}
  
  ### Initialize precision matrices and graphs
  # # True values
  # C = data_input$param_true$Omega_b_true
  # Sig = list()
  # for (s_i in 1:2){Sig[[s_i]] = solve(C[[s_i]])}
  # adj = data_input$param_true$G_b_true
  
  # Use diag
  C = list(); adj = list()
  for (s_i in 1:2){C[[s_i]] = diag(p_all); adj[[s_i]] = diag(TRUE, p_all)}
  Sig = C
  
  #----------------------------------------------------------------------------
  # Run the MCMC algorithm
  #---------------------------------------------------------------------------
  #mcmc_output = MCMC_changepoint_newB(nburn, nsave,
  mcmc_output = MCMC_changepoint_DBFGM(Y, 
                                       K, FLC, 
                                       v0, v1, a_pi, b_pi,
                                       changepoint_interval,
                                       changepoint,
                                       B, 
                                       Sig, C, adj, pii_block,
                                       sigma_epsilon,
                                       nburn, nsave, 
                                       disp);
  
  # Save parameters and initialization
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

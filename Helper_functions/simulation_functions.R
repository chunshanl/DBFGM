### Functions to simulate data
simulate_blocked_Omega = function(p, K, pii_local){
  
  ## Simulate blocked G and Omega -----------------------------
  # Omega is p by p blocks; each block is size K by K
  # p: number of nodes in the global graph
  # K: the block size is K by K
  # pii_local: between 0-1, block-wise edge inclusion probability
  
  ## Outputs:
  # G_x: the p by p global graph
  # G_b: the blocked graph
  # Omega_b: the blocked precision matrix, 1s on the diagonal
  
  # Generate G_x, the p-by-p global graph
  prob = 2/(p - 1)
  G_x = matrix(0, p, p);
  G_x[upper.tri(G_x)] <- rbinom(p * (p - 1)/2, 
                                1, prob)  
  diag(G_x) = 1

  # Generate G_b, the graph of size p*K by p*K 
  # pii_local = 1/2
  G_b = matrix(0, p * K, p * K)
  G_b[upper.tri(G_b)] <- rbinom(p * K * (p * K - 1)/2, 
                                1, pii_local)  
  diag(G_b) = 1
  
  temp = kronecker(G_x, matrix(1, K, K))  # Add block structure according to G_x
  G_b = temp * G_b
  # pheatmap(G_b, cluster_cols = FALSE, cluster_rows = FALSE)
  
  G_x = matrix(as.logical(G_x), p, p)
  G_x = G_x | t(G_x)
  # pheatmap(G_x+0, cluster_cols = FALSE, cluster_rows = FALSE)
  G_b = matrix(as.logical(G_b), p*K, p*K)
  G_b = G_b | t(G_b)
  # pheatmap(G_b+0, cluster_cols = FALSE, cluster_rows = FALSE)
  
  # Generate the sparse precision matrix 
  Omega_b = rgwish(adj = G_b)
  Omega_b = cov2cor(Omega_b)
  Omega_b[!G_b] = 0  # p * p blocks
  # pheatmap(Omega_b, cluster_rows=FALSE, cluster_cols=FALSE)
  
  # Output 
  output = list(G_x, G_b, Omega_b, pii_local)
  names(output) <- c("G_x_true", "G_b_true", "Omega_b_true", "pii_local")
  return(output)
}


simulate_mvn_blocked =  function(p = 10,
                        K,
                        n = 100,
                        pii_local){
  # Simulate blocked Omega and simulate MVN data given Omega--------------
  
  out = simulate_blocked_Omega(p, K, pii_local) 
  
  Y = rmvnorm(n = n, sigma = solve(out$Omega_b_true))
  
  output = list(out$G_x_true, out$G_b_true, out$Omega_b_true, Y, pii_local, p, K)
  
  names(output) <- c("G_x_true", "G_b_true", "Omega_b_true", "Y", "pii_local", p, K)
  
  return(output)
  
}


simulate_static_data = function(p = 10,
                                K,
                                n = 100,
                                T_data = 100,
                                pii_local,
                                sigma_epsilon_true = 0.05,
                                basis_type_true,
                                data = NULL,
                                disp){
  
  # Simulate multivariate functional data  -----------------------------------------------------------------------
  # if (basis_type_true == 'spline'){
  #   selfcor = T
  # }else{
  #   selfcor = FALSE
  # }
  
  if (!is.null(data)){
    G_x = data$G_x_true
    G_b = data$G_b_true
    Omega_b = data$Omega_b_true
    B_true = data$Y
  }else{
    # Generate blocked precision matrix ------------------------------------
    # out = simulate_blocked_Omega(p, K, selfcor) 
    out = simulate_blocked_Omega(p, K, pii_local) 
    G_x = out$G_x_true; G_b = out$G_b_true; Omega_b = out$Omega_b_true
    # if (disp){
    #   pheatmap(Omega_b, cluster_rows=FALSE, cluster_cols=FALSE)
    # }
    
    # Generate factors ------------------------------------
    B_true = rmvnorm(n = n, sigma = solve(Omega_b)) # size n * (pK)
  }
  
  # Generate Factor Loading Curves -------------------------------------------------------
  # Observation points:
  U = seq(0, 1, length.out = T_data)

  if (basis_type_true == 'spline'){
    knots = U[seq(0, length(U), length.out = K-1)]
    b = create.bspline.basis(rangeval = c(0,1), breaks = knots, norder = 4)
    F_true = eval.basis(U, b)   # T_data * K
    
    # jump point trying, to be deleted
    # U = seq(0, 1, length.out = 200)
    # knots = U[seq(0, 200, length.out = K-1)]
    # b = create.bspline.basis(rangeval = c(0,1), breaks = knots, norder = 4)
    # F_true = eval.basis(U, b)
    # F_true = F_true[1:99, ]
    
  }else if (basis_type_true == 'polynomial'){
    # FLCs: orthonormalized polynomials
    F_true = cbind(1/sqrt(T_data),
                   poly(U, K - 1))  # size T_data*K
    # F_true = cbind(1/sqrt(200),
    #                poly(U, K - 1))  # size T_data*K
    # F_true = F_true[1:100,]
  }
  
  # Compute p functions --------------------------------------------
  X = array(NA, c(n, T_data, p))
  for (i in 1:p){
    b = B_true[, ((i-1)*K+1):(i*K)]  # b is of size n * K
    X[,,i] = b %*% t(F_true)
  }

  # Add noise -------------------------------------------------
  noise = array( rnorm(n*T_data*p, mean = 0, sd = sigma_epsilon_true), 
                 c(n,T_data,p))
  Y = X + noise
  if (disp){
    matplot(F_true, type = 'l', xlab = 'True Basis Functions')
    matplot(X[1,,], type = 'l', xlab = 'One sample of p functions (without noise)')
    matplot(Y[1,,], type = 'l', xlab = 'One sample of p functions (with noise)')
  }
  
  
  # Output -----------------------------
  output = list(G_x, G_b, Omega_b, B_true, 
                F_true, Y, sigma_epsilon_true, pii_local)
  
  names(output) <- c("G_x_true", "G_b_true", "Omega_b_true", "B_true", 
                     "F_true", "Y", "sigma_epsilon_true", 'pii_local')
  
  return(output)
  
}


simulate_changepoint_data = function(p,
                              K_true,
                              n,
                              T_data,
                              pii_local,
                              sigma_epsilon_true,
                              basis_type_true,
                              continuous = FALSE,
                              changepoint_true){
  ## Simulate multivariate functional data with change of functional graph -------------
  
  # Inputs:
  # p: p random functions in one observation
  # K_true: true number of basis used in data simulation
  # n: sample size; T_data: number of observed points of each function
  # pii_local: probability of edge inclusion probability in each block when generating true precision matrices
  # sigma_epsilon_true: sigma_epsilon (std of random noise) used for data generation
  # basis_type_true: spline or polynomial
  # continuous: True or False: whether to adjust the functions before and after the change point to be continuous at the change point
  # p by p blocks, each block is size K by K
  # changepoint_true: true change point
  
  ## Outputs:
  # param_true: list of set-ups used to generate the data
  # Y: simulated data, size n by T by p
  
  # Save set-up
  param_true = list()
  param_true$changepoint_true = changepoint_true
  param_true$pii_local = pii_local
  param_true$sigma_epsilon_true = sigma_epsilon_true
  param_true$K_true = K_true
  
  # Simulate precision matrices and coefficients
  n_interval = length(changepoint_true) + 1
  for (i_interval in 1:n_interval){
    out = simulate_blocked_Omega(p, K_true, pii_local)  # G_x_true, G_b_true, Omega_b_true, pii_local
    param_true$G_x_true[[i_interval]] = out$G_x_true 
    param_true$G_b_true[[i_interval]] = out$G_b_true
    param_true$Omega_b_true[[i_interval]] = out$Omega_b_true
    B_true = rmvnorm(n = n, sigma = solve(out$Omega_b_true)) # size n * (pK)
    param_true$B_true[[i_interval]] = B_true
  }

  # Generate Factor Loading Curves 
  U = seq(0, 1, length.out = T_data) 
  if (basis_type_true == 'spline'){
    knots = U[seq(0, length(U), length.out = K_true-1)]
    b = create.bspline.basis(rangeval = c(0,1), breaks = knots, norder = 4)
    F_true = eval.basis(U, b)   # T_data * K
  }else if (basis_type_true == 'polynomial'){
    # FLCs: orthonormalized polynomials
    F_true = cbind(1/sqrt(T_data),
                   poly(U, K_true - 1))  # size T_data*K
  }else if (basis_type_true == 'spline_change_point'){
    F_true = c()
    # before
      U_temp = seq(0, 1, length.out = changepoint_true - 1) 
      knots = U_temp[seq(0, length(U_temp), length.out = K_true-1)]
      b = create.bspline.basis(rangeval = c(0,1), breaks = knots, norder = 4)
      F_true_temp = eval.basis(U_temp, b)  
      F_true = rbind(F_true, F_true_temp)
    # after
      U_temp = seq(0, 1, length.out = T_data - changepoint_true + 1) 
      knots = U_temp[seq(0, length(U_temp), length.out = K_true-1)]
      b = create.bspline.basis(rangeval = c(0,1), breaks = knots, norder = 4)
      F_true_temp = eval.basis(U_temp, b) 
      F_true = rbind(F_true, F_true_temp)
  } else if (basis_type_true == 'polynomial_change_point'){
    F_true = c()
    # before
    U_temp = seq(0, 1, length.out = changepoint_true - 1) 
    F_true_temp = cbind(1/sqrt(changepoint_true - 1),
                   poly(U_temp, K_true - 1))  
    F_true = rbind(F_true, F_true_temp)
    # after
    U_temp = seq(0, 1, length.out = T_data - changepoint_true + 1)
    F_true_temp = cbind(1/sqrt(T_data - changepoint_true + 1),
                        poly(U_temp, K_true - 1))  
    F_true = rbind(F_true, F_true_temp)
  }
  # matplot(F_true, type = 'l') 
  
  # Multiply factor loading curves with coefficients 
  interval_ind = matrix(NA, 2, 2)
  interval_ind[1,] = c(1, changepoint_true - 1)
  interval_ind[2,] = c(changepoint_true, T_data)
  X = compute_X(param_true$B_true, F_true, interval_ind, p)
  # plot(X[1,,1])
  # matplot(X[1,,], type = 'l')
  
  # Continuous at jump point constraint
  if (continuous){
    X_noadjust = X
    for (i in 1:n){
      temp = (X[i, changepoint_true-1,] - X[i, changepoint_true,])/2
      nrows = T_data - changepoint_true + 1
      temp1 =  matrix(rep(temp,each = nrows),nrow = nrows)
      X[i, changepoint_true:T_data,] = X[i, changepoint_true:T_data,] + temp1
      nrows = changepoint_true - 1
      temp1 =  matrix(rep(temp,each = nrows),nrow = nrows)
      X[i, 1:(changepoint_true-1),] = X[i, 1:(changepoint_true-1),] - temp1
    }
  }
  
  param_true$F_true = F_true
  
  # Add noise -------------------------------------------------
  noise = array( rnorm(n*T_data*p, mean = 0, sd = sigma_epsilon_true), 
                 c(n, T_data, p))
  Y = X + noise
  # matplot(Y[1,,], type = 'l')
  # matplot(X[1,,], type = 'l')
  
  # Output -----------------------------
  output = list(param_true, Y)
  names(output) <- c("param_true", "Y")
  # if (continuous){
  #   output$B_true_before_noadjust = B_true_before_noadjust
  #   output$B_true_after_noadjust = B_true_after_noadjust
  # }
  # 
  return(output)
  
}


simulate_changepoint_data_replications = function(data_changepoint,
                                                  continuous = FALSE,
                                                  random_seed){
  
  # data_changepoint: data generated from the simulate_changepoint_data function
  # random_seed: random seed for this replication
  
  # Save parameter
  param_true = data_changepoint$param_true
  
  n = dim(data_changepoint$Y)[1]; T_data = dim(data_changepoint$Y)[2]
  p = dim(data_changepoint$Y)[3]; K_true = data_changepoint$param_true$K_true
  n_interval = 2
  F_true = param_true$F_true; changepoint_true = param_true$changepoint_true
  sigma_epsilon_true = param_true$sigma_epsilon_true
  
  set.seed(random_seed)
  # Simulate coefficients
  for (i_interval in 1:n_interval){
    B_true = rmvnorm(n = n, sigma = solve(param_true$Omega_b_true[[i_interval]])) # size n * (pK)
    param_true$B_true[[i_interval]] = B_true
  }
  
  # Multiply factor loading curves with coefficients 
  interval_ind = matrix(NA, 2, 2)
  interval_ind[1,] = c(1, changepoint_true - 1)
  interval_ind[2,] = c(changepoint_true, T_data)
  X = compute_X(param_true$B_true, F_true, interval_ind, p)
  # X = array(NA, c(n, T_data, p))
  # for (i in 1:p){
  #   b = param_true$B_true[[1]][, ((i-1)*K_true+1):(i*K_true)]  # B_true is n*pk, b is of size n * K
  #   X[, 1:(changepoint_true-1), i] = b %*% t(F_true[1:(changepoint_true-1), ])
  #   
  #   b = param_true$B_true[[2]][, ((i-1)*K_true+1):(i*K_true)]  
  #   X[, changepoint_true:T_data, i] = b %*% t(F_true[changepoint_true:T_data, ])
  # }
  
  # plot(X[1,,1])
  # matplot(X[1,,], type = 'l')
  
  # Continuous at jump point constraint
  if (continuous){
    X_noadjust = X
    for (i in 1:n){
      temp = (X[i, changepoint_true-1,] - X[i, changepoint_true,])/2
      nrows = T_data - changepoint_true + 1
      temp1 =  matrix(rep(temp,each = nrows),nrow = nrows)
      X[i, changepoint_true:T_data,] = X[i, changepoint_true:T_data,] + temp1
      nrows = changepoint_true - 1
      temp1 =  matrix(rep(temp,each = nrows),nrow = nrows)
      X[i, 1:(changepoint_true-1),] = X[i, 1:(changepoint_true-1),] - temp1
    }
    # adjusted_B = continuous_adjust(F_true, changepoint_true,
    #                                   solve(Omega_b_true_before),
    #                                   solve(Omega_b_true_after),
    #                                   B_true_before_noadjust, B_true_after_noadjust)
  }
  
  
  # Add noise -------------------------------------------------
  noise = array( rnorm(n*T_data*p, mean = 0, sd = sigma_epsilon_true), 
                 c(n, T_data, p))
  Y = X + noise
  # matplot(Y[1,,], type = 'l')
  # matplot(X[2,,], type = 'l')
  
  # Output -----------------------------
  output = list(param_true, Y)
  names(output) <- c("param_true", "Y")
  return(output)
  
}

simulate_cluster_data = function(p,
                              K_true,
                              n,
                              T_data,
                              pii_local,
                              sigma_epsilon_true,
                              basis_type_true,
                              continuous = FALSE,
                              changepoint_true,
                              n_cluster_true){
  
  param_true = list()
  param_true$K_true = K_true
  param_true$pii_local = pii_local
  param_true$sigma_epsilon_true = sigma_epsilon_true
  param_true$basis_type_true = basis_type_true
  param_true$changepoint_true = changepoint_true
  param_true$n_cluster_true = n_cluster_true

  n_interval = length(changepoint_true) + 1

  # Simulate random cluster assignments
  cluster_ind_true = list()  
  for (ii in 1:n_interval){
    cluster_ind_true[[ii]] = sample(x=1:n_cluster_true[ii], size = n, prob = rep(1, n_cluster_true[ii]), replace =TRUE)
  }
  cluster_ind_true[[2]] = cluster_ind_true[[2]]+2  ################### CHANGE
  param_true$cluster_ind_true = cluster_ind_true
  
  # Simulate Omega and coefficients
  for (ii in 1:sum(n_cluster_true)){
    out = simulate_blocked_Omega(p, K_true, pii_local) 
    param_true$G_x_true[[ii]] = out$G_x_true  # G_x_true, G_b_true, Omega_b_true, pii_local
    param_true$G_b_true[[ii]] = out$G_b_true
    param_true$Omega_b_true[[ii]] = out$Omega_b_true
  }
  B_true = list()  
  for (ii in 1:n_interval){
    B_true[[ii]] = array(NA, dim = c(n, p*K_true))
    for (jj in 1:sum(n_cluster_true)){
      n_temp = sum(cluster_ind_true[[ii]] == jj)
      B_true_temp = rmvnorm(n = n_temp, sigma = solve(param_true$Omega_b_true[[jj]])) # size n * (pK)
      B_true[[ii]][cluster_ind_true[[ii]]==jj, ] = B_true_temp
    }
  }
  param_true$B_true = B_true
  
  # Generate Factor Loading Curves 
  U = seq(0, 1, length.out = T_data) # Observation points:
  if (basis_type_true == 'spline'){
    knots = U[seq(0, length(U), length.out = K_true-1)]
    b = create.bspline.basis(rangeval = c(0,1), breaks = knots, norder = 4)
    F_true = eval.basis(U, b)   # T_data * K
    
    # jump point trying, to be deleted
    # U = seq(0, 1, length.out = 200)
    # knots = U[seq(0, 200, length.out = K-1)]
    # b = create.bspline.basis(rangeval = c(0,1), breaks = knots, norder = 4)
    # F_true = eval.basis(U, b)
    # F_true = F_true[1:99, ]
    
  }else if (basis_type_true == 'polynomial'){
    # FLCs: orthonormalized polynomials
    F_true = cbind(1/sqrt(T_data),
                   poly(U, K_true - 1))  # size T_data*K
    # F_true = cbind(1/sqrt(200),
    #                poly(U, K - 1))  # size T_data*K
    # F_true = F_true[1:100,]
  }
  # matplot(F_true, type = 'l') 
  
  # Multiply factor loading curves with coefficients 
  X = array(NA, c(n, T_data, p))
  for (i in 1:p){
    b = param_true$B_true[[1]][, ((i-1)*K_true+1):(i*K_true)]  # B_true is n*pk, b is of size n * K
    X[, 1:(changepoint_true-1), i] = b %*% t(F_true[1:(changepoint_true-1), ])
    
    b = param_true$B_true[[2]][, ((i-1)*K_true+1):(i*K_true)]  
    X[, changepoint_true:T_data, i] = b %*% t(F_true[changepoint_true:T_data, ])
  }
  
  # plot(X[1,,1])
  # matplot(X[1,,], type = 'l')
  
  # Continuous at jump point constraint
  if (continuous){
    X_noadjust = X
    for (i in 1:n){
      temp = (X[i, changepoint_true-1,] - X[i, changepoint_true,])/2
      nrows = T_data - changepoint_true + 1
      temp1 =  matrix(rep(temp,each = nrows),nrow = nrows)
      X[i, changepoint_true:T_data,] = X[i, changepoint_true:T_data,] + temp1
      nrows = changepoint_true - 1
      temp1 =  matrix(rep(temp,each = nrows),nrow = nrows)
      X[i, 1:(changepoint_true-1),] = X[i, 1:(changepoint_true-1),] - temp1
    }
  }
  param_true$F_true = F_true
  
  # Add noise -------------------------------------------------
  noise = array( rnorm(n*T_data*p, mean = 0, sd = sigma_epsilon_true), 
                 c(n, T_data, p))
  Y = X + noise
  # matplot(Y[1,,], type = 'l')
  # matplot(X[2,,], type = 'l')
  
  # Output -----------------------------
  output = list(param_true, Y)
  names(output) <- c("param_true", "Y")
  # if (continuous){
  #   output$B_true_before_noadjust = B_true_before_noadjust
  #   output$B_true_after_noadjust = B_true_after_noadjust
  # }
  return(output)
  
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



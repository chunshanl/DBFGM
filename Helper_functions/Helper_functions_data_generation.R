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
  G_b = matrix(0, p * K, p * K)
  G_b[upper.tri(G_b)] <- rbinom(p * K * (p * K - 1)/2, 
                                1, pii_local)  
  diag(G_b) = 1
  
  temp = kronecker(G_x, matrix(1, K, K))  # Add block structure according to G_x
  G_b = temp * G_b
  
  G_x = matrix(as.logical(G_x), p, p)
  G_x = G_x | t(G_x)
  G_b = matrix(as.logical(G_b), p*K, p*K)
  G_b = G_b | t(G_b)
  
  # Generate the sparse precision matrix 
  Omega_b = rgwish(adj = G_b)
  Omega_b = cov2cor(Omega_b)
  Omega_b[!G_b] = 0  # p * p blocks
  
  # Output 
  output = list(G_x, G_b, Omega_b, pii_local)
  names(output) <- c("G_x_true", "G_b_true", "Omega_b_true", "pii_local")
  return(output)
}


simulate_data = function(p,
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
  }
  
  # Multiply factor loading curves with coefficients 
  interval_ind = matrix(NA, n_interval, 2)
  interval_ind[1:length(changepoint_true), 2] = changepoint_true -1
  interval_ind[2:n_interval, 1] = changepoint_true
  interval_ind[1,1] = 1
  interval_ind[n_interval,2] = T_data
  
  X = compute_X(param_true$B_true, F_true, interval_ind, p)  # n by T by p
  
  # Continuous at jump point constraint
  if (continuous){
    X_noadjust = X
    for (i in 1:n){
      for (changepoint_i in 1:length(changepoint_true)){
        changepoint = changepoint_true[changepoint_i]
        adjust_val = (X[i, changepoint-1,] - X[i, changepoint,])/2  # vec of length p
        
        ind_1 = 1:(changepoint-1)
        nrows = length(ind_1)
        adjust_val_matrix =  matrix(rep(adjust_val,each = nrows),nrow = nrows)  # number of time points by p, each col is the same value
        X[i, ind_1,] = X[i, ind_1,] - adjust_val_matrix
        
        ind_2 = interval_ind[changepoint_i+1,1]:interval_ind[changepoint_i+1,2] # For the first change point, want the second interval
        nrows = length(ind_2)
        adjust_val_matrix =  matrix(rep(adjust_val,each = nrows),nrow = nrows)
        X[i, ind_2,] = X[i, ind_2,] + adjust_val_matrix
        
      }
    }
  }
  
  param_true$F_true = F_true
  
  # Add noise -------------------------------------------------
  noise = array( rnorm(n*T_data*p, mean = 0, sd = sigma_epsilon_true), 
                 c(n, T_data, p))
  Y = X + noise
  
  # Output -----------------------------
  output = list(param_true, Y)
  names(output) <- c("param_true", "Y")
  
  return(output)
  
}


simulate_data_replications = function(data_changepoint, 
                                      continuous = FALSE,
                                      random_seed){
  
  # data_changepoint: data generated from the simulate_changepoint_data function
  # random_seed: random seed for this replication
  
  # Save parameter
  param_true = data_changepoint$param_true
  
  n = dim(data_changepoint$Y)[1]; T_data = dim(data_changepoint$Y)[2]
  p = dim(data_changepoint$Y)[3]; K_true = data_changepoint$param_true$K_true
  n_interval = length(param_true$changepoint_true)+1
  F_true = param_true$F_true; changepoint_true = param_true$changepoint_true
  sigma_epsilon_true = param_true$sigma_epsilon_true
  
  set.seed(random_seed)
  # Simulate coefficients
  for (i_interval in 1:n_interval){
    B_true = rmvnorm(n = n, sigma = solve(param_true$Omega_b_true[[i_interval]])) # size n * (pK)
    param_true$B_true[[i_interval]] = B_true
  }
  
  # Multiply factor loading curves with coefficients 
  interval_ind = matrix(NA, n_interval, 2)
  interval_ind[1:length(changepoint_true), 2] = changepoint_true -1
  interval_ind[2:n_interval, 1] = changepoint_true
  interval_ind[1,1] = 1
  interval_ind[n_interval,2] = T_data
  
  X = compute_X(param_true$B_true, F_true, interval_ind, p)  # n by T by p
  
  # Continuous at jump point constraint
  if (continuous){
    X_noadjust = X
    for (i in 1:n){
      for (changepoint_i in 1:length(changepoint_true)){
        changepoint = changepoint_true[changepoint_i]
        adjust_val = (X[i, changepoint-1,] - X[i, changepoint,])/2  # vec of length p
        
        ind_1 = 1:(changepoint-1)
        nrows = length(ind_1)
        adjust_val_matrix =  matrix(rep(adjust_val,each = nrows),nrow = nrows)  # number of time points by p, each col is the same value
        X[i, ind_1,] = X[i, ind_1,] - adjust_val_matrix
        
        ind_2 = interval_ind[changepoint_i+1,1]:interval_ind[changepoint_i+1,2] # For the first change point, want the second interval
        nrows = length(ind_2)
        adjust_val_matrix =  matrix(rep(adjust_val,each = nrows),nrow = nrows)
        X[i, ind_2,] = X[i, ind_2,] + adjust_val_matrix
        
      }
    }
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

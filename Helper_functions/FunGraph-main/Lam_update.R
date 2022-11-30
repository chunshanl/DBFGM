lam.update = function(C, lam_prior){
    # funcion to update lasso parameter lambda
    # input: C - current correlation matrix with diagonals 1
    #        lam_prior - gamma parameter for lambda prior distirubtion
    
    p <- nrow(C)
    rho <- sapply(C[upper.tri(C)], C2rho)
    
    lam = rgamma(1,lam_prior$a+p*(p-1)/2, rate=lam_prior$b+sum(abs(rho)));

    return(lam)
}

## define function to transform C_ij to rho_ij
C2rho <- function(C_ij){
    rho_ij  <- -C_ij/(1-C_ij^2)
    return(rho_ij)
}

# ## test/debug
# Sig <- toeplitz(0.7^(0:9))
# O <- solve(Sig)
# D <- diag(diag(O)^.5)
# C.true <- round(solve(D) %*% O %*% solve(D),4)
# 
# p=10; B=1000; 
# rho.samp <- array(NA,B)
# for(it in 1:B){
#     rho.samp[it] <- lam.update(C.true,lam_prior=list(a=0.1,b=0.1))
# }
# mean(rho.samp)
# plot(rho.samp,type='l')

D.update = function(n,S,C,D.diag,D_prior) {
    # function to update correlation matrix C
    # input: n - sample size
    #        C - correlation matrix with diagonals 1
    #        S - sample covariance 
    #        D.diag - current diagonal entries of conditional sds
    #        D_prior - gamma parameters for D_jj prior distribution

    p <- nrow(C)
    a <- D_prior$a
    b <- D_prior$b
    
    # instead of updating D_jj separately, which moves very slowly due to high correlations
    # propose all D diagnal entries together within Omega.hat as arising from Glasso
    D <- diag(D.diag)
    O <- D %*% C %*% D
    for(j in 1:p) {
        noj <- setdiff(1:p,j)
        O.pm <- O[c(noj,j),c(noj,j)]
        
        R <- chol(O.pm)
        c <- sum(R[-p,p]^2)
        
        O[j,j] <- rgamma(1,n/2+a,S[j,j]/2+b) + c
    }
    D.star <- diag(diag(O)^.5)
    
    log.alpha <- D.post(D.star,S,C,n,a,b) - D.post(D,S,C,n,a,b)
    if(log(runif(1)) < log.alpha) D.diag <- diag(D.star)
    
    return(D.diag)
}

## define the function to calculate log posterior of D
D.post <- function(D,S,C,n,a,b){
    tt <- (n+2*a-1)*log(det(D)) - sum(diag(S %*% D %*% C %*% D))/2 -2*b
    return(tt)
}

# ## test/debug
# Sig <- toeplitz(0.7^(0:9))
# O.true <- solve(Sig)
# D.true <- diag(O.true)^.5
# D <- diag(D.true)
# C.true <- round(solve(D) %*% O.true %*% solve(D),4)
# 
# chol(O.true)
# chol(C.true) %*% D
# 
# n = 500; p = 10; B=5000;
# y <- matrix(rnorm(n*p),n,p) %*% chol(Sig)
# S <- t(y)%*%y
# 
# D.samp <- array(NA,c(p,B))
# D.cur <- rep(NA,p)
# for(j in 1:p) D.cur[j] <- summary(lm(y[,j] ~ y[,-j]))$sigma^-1
# C.samp <- array(NA,c(p,p,B))
# O.cur <- solve(cov(y)+diag(0.1,p))
# C.cur <- diag(diag(O.cur)^-0.5) %*% O.cur %*% diag(diag(O.cur)^-0.5)
# lam.samp <- rep(NA,B)
# lam.cur <- 0.1
# for(it in 1:B){
#     # D.samp[,it] <- D.update(n,S,C=C.true,D.diag=D.cur,D_prior=list(a=0.1,b=0.1))
#     D.cur <- D.update(n,S,C=C.cur,D.diag=D.cur,D_prior=list(a=0.1,b=0.1))
#     C.cur <- C.update(n,S,C=C.cur,D=D.cur,lam=lam.cur)
#     lam.cur <- lam.update(C=C.cur,lam_prior=list(a=1,b=11))
#     D.samp[,it] <- D.cur
#     C.samp[,,it] <- C.cur
#     lam.samp[it] <- lam.cur
# }
# apply(D.samp[,1001:5000],1,mean)
# D.true
# plot(D.samp[1,],type='l')
# # apply(C.samp[,,1001:5000],1:2,mean) * (apply(C.samp[,,1001:5000],1:2,quantile,0.025) > 0
# # | apply(C.samp[,,1001:5000],1:2,quantile,0.975) < 0)
# # C.true
# # plot(lam.samp,type='l')

# tmp <- apply(C.samp,1:2,mean)
# tmp2 <- apply(D.samp,1,mean)
# O.est <- diag(tmp2) %*% tmp %*% diag(tmp2)
# Sig.est <- solve(O.est)


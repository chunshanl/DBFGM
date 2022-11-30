C.update = function(n,S,C,D.diag,lam) {
    # function to update correlation matrix C
    # input: n - sample size
    #        C - current correlation matrix with diagonals 1
    #        S - sample covariance 
    #        D.diag - diagonal entries of conditional sds
    #        lam - lasso parameter

    p <- nrow(C)
    D <- diag(D.diag)
    DSD <- D %*% S %*% D
    
    for(i in 1:(p-1)) {
        for(j in (i+1):p) {
            
            # move i,j to last two rows/columns
            noij <- setdiff(1:p,c(i,j))
            C.pm <- C[c(noij,i,j),c(noij,i,j)]
            DSD.pm <- DSD[c(noij,i,j),c(noij,i,j)]
            
            R <- chol(C.pm)
            a = sum(R[1:(p-2),p-1]*R[1:(p-2),p])
            bc = (R[p-1,p-1]^2)*(R[p,p]^2 + R[p-1,p]^2)
            DSDij = DSD.pm[p-1,p]
            
            Cij.grid <- seq(a-sqrt(bc)*0.99,a+sqrt(bc)*0.99,length.out=100)
            post.grid <- Cij.post(Cij.grid,a,bc,DSDij,n,lam)
            Cij.star <- sample(Cij.grid, size=1, prob=exp(post.grid - max(post.grid))) + runif(1,-1,1)*sqrt(bc)/100
            
            log.alpha <- Cij.post(Cij.star,a,bc,DSDij,n,lam) - Cij.post(C.pm[p-1,p],a,bc,DSDij,n,lam)
            if(log(runif(1)) < log.alpha) C[i,j] <- C[j,i] <- Cij.star
        }
    }
        
    return(C)
}

## define function to calculate log posterior density of C_ij
Cij.post <- function(Cij,a,bc,DSDij,n,lam){
    t1 <- n/2*log(1 - (Cij-a)^2/(bc))
    t2 <- -DSDij*Cij
    t3 <- -lam*abs(Cij)/(1-Cij^2) + log(1+Cij^2) - 2*log(1-Cij^2)
    return(t1+t2+t3)
}

# ## test/debug
# Sig <- toeplitz(0.7^(0:9))
# O <- solve(Sig)
# D <- diag(diag(O)^.5)
# C.true <- round(solve(D) %*% O %*% solve(D),4)
# 
# chol(O)
# chol(C.true) %*% D
# 
# n = 100; p = 10; B=5000;
# y <- matrix(rnorm(n*p),n,p) %*% chol(Sig)
# 
# C.samp <- array(NA,c(p,p,B))
# O.cur <- solve(cov(y))
# C.cur <- diag(diag(O.cur)^-0.5) %*% O.cur %*% diag(diag(O.cur)^-0.5)
# for(it in 1:B){
#     C.cur <- C.update(n=n,S=t(y)%*%y,C=C.cur,D=D,lam=2)
#     C.samp[,,it] <- C.cur
# }
# apply(C.samp[,,3001:5000],1:2,mean)
# C.true
# plot(C.samp[1,2,],type='l')


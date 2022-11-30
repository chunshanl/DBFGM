## load packages and R functions
library(wavelets)
library(doParallel)
library(foreach)
#registerDoParallel(24)
source('fgraph_cclasso.R')
source('Cmat_update.R')
source('Dmat_update.R')
source('Lam_update.R')


#### specify graphical model
p = 10
G0 <- diag(0,p)
G0[1,2] <- G0[1,3] <- G0[2,3] <- G0[7,8] <- G0[7,9] <- G0[8,9] <- 1 
G1 <- G2 <- G0
G1[4,5] <- G1[5,6] <- G1[4,6] <- 1
G2[7,10] <- G2[8,10] <- G2[9,10] <- 1
G0 <- G0 + t(G0); G1 <- G1 + t(G1); G2 <- G2 + t(G2)

## specify precision matrix of error and auto-regressive matrix A
T = 128
tt = 1:T
r.max <- -0.99
r1 <- sapply(tt,function(t) min(max(0,t-T/4),T/2)*(r.max/(T/2))) 
r2 <- r.max - r1
             
Sigma.err <- rep(list(diag(p)),T)
for(t in 1:T) {
    Omega.tmp <- diag(c(rep(2,6),rep(3,4)))
    Omega.tmp[which(G1==1)] <- r1[t]
    Omega.tmp[which(G2==1)] <- Omega.tmp[which(G2==1)] + r2[t]
    Sigma.tmp <- solve(Omega.tmp)
    Sigma.err[[t]] <- diag(diag(Sigma.tmp)^-0.5) %*% Sigma.tmp %*% diag(diag(Sigma.tmp)^-0.5)*2
}
A <- diag(1,p)

## true conditional cross-covariance graph
Sigmas1.true <- G1.true <- G1.neg <- rep(list(diag(p)),T)
for(t in 1:T){
    if(t == 1) {
        Sigmas1.true[[t]] <- Sigma1 <- Sigma.err[[1]] 
    } else {
        Sigmas1.true[[t]] <- Sigma1 <- A %*% Sigma1 %*% t(A) * 0.5 + Sigma.err[[t]] * 0.5
    }
    G1.true[[t]] <- (abs(solve(Sigma1)) > 0.1)*1
    G1.neg[[t]] <- (abs(solve(Sigma1)) < 0.001)*1
}
G1.true <- abind(G1.true,along=3)
G1.neg <- abind(G1.neg,along=3)


#### generate basis functions
tmp = diag(T)
W = matrix(NA,T,T)
for(t in 1:T) {
    a = dwt(tmp[t,],filter='haar',n.levels=5)
    W[t,] = c(a@V[[5]],unlist(a@W[5:1]))
}
Phi <- solve(W)

    
####generate data
set.seed(9999)

n <- 50; p=10; T=128
y <- array(NA,c(n,p,T))

for(t in 1:T){
    err <- matrix(rnorm(n*p),n) %*% chol(Sigma.err[[t]])
    if(t == 1) y[,,t] <- err else
        y[,,t] <- t(apply(y[,,t-1],1,function(x) A%*%x))*sqrt(0.5) + err*sqrt(0.5)
}

    
#### run Bayesian functional graph
  
MCMCspecs = list(B=1000,thin=5,burnin=1000,update=1000);
D_prior = list(a=0.1,b=0.1); lam_prior = list(a=0.1,b=1);
fgbay.samp <- fgraph_ccc(y,Phi,MCMCspecs,D_prior,lam_prior)

fgbay.est <- array(diag(p),c(p,p,T))
for(t in 1:T) fgbay.est[,,t][upper.tri(diag(p))] <- apply(fgbay.samp$ccc[,t,],1,function(x) (quantile(x,0.025)>0 | quantile(x,0.975)<0)*1)


#### calculate TPs and FPs

tmp <- matrix(FALSE,p,p); tmp[upper.tri(tmp)] <- TRUE
uptri.ind <- array(tmp,c(p,p,T))

## TPs
sum((fgbay.est+G1.true==2)[uptri.ind])/sum(G1.true[uptri.ind]) 
## FPs
sum((fgbay.est+G1.neg==2)[uptri.ind])/sum(G1.neg[uptri.ind])
    

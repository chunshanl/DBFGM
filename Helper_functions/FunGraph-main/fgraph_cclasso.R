fgraph_ccc = function(Y,Phi,MCMCspecs,D_prior,lam_prior) {
# This is the non-parallel computing version
# Note: functions to conduct functional network inference using lasso on conditional cross-covariance
# Input: Y - (nXpxT) data matrix, p: # of variables; T: # of functional location; n: sample size
#        Phi - (KxT) basis function matrix
#        MCMCspecs - parameters for MCMC iterations
#        D_prior - list of a&b giving gamma parameters for D_jj^2 prior distribution
#        lam_prior - list of a&b giving gamma parameters for lambda prior distribution
#
# Output:res - data structure with elements: 

library(abind)
library(plyr)
library(MASS)
library(glasso)
# library(doParallel)
# library(foreach)
# if(!is.null(ncor)) registerDoParallel(ncor)

### Step 1: Apply DWT to project observed functions into basis dual space
n = dim(Y)[1]; p = dim(Y)[2]; T = dim(Y)[3];
K = nrow(Phi); W = ginv(Phi) #t(Phi) %*% solve(Phi %*% t(Phi));

y.dual = alply(abind(sapply(alply(Y,1),function(x) x%*%W, simplify=FALSE),along=0),3)
S.dual <- sapply(y.dual, function(x) t(x)%*%x,simplify=FALSE)


### Step 2: initial value
Dd.dual <- rep(list(rep(NA,p)),K)
C.dual <- rep(list(diag(p)),K)
for(k in 1:K) {
    # for(j in 1:p) Dd.dual[[k]][j] <- summary(lm(y.dual[[k]][,j] ~ y.dual[[k]][,-j]))$sigma^-1
    # O.hat <- solve(S.dual[[k]]/n+0.01*diag(p))
    # C.dual[[k]] <- diag(diag(O.hat)^-0.5) %*% O.hat %*% diag(diag(O.hat)^-0.5)

    tmp <- glasso(S.dual[[k]]/n,rho=0.1)$wi
    O.hat <- (tmp + t(tmp))/2
    Dd.dual[[k]] <- diag(O.hat)^.5
    C.dual[[k]] <- diag(diag(O.hat)^-0.5) %*% O.hat %*% diag(diag(O.hat)^-0.5)
}
lam.dual <- rep(0.1,K)


### Step 3: MCMC sampling
B = MCMCspecs$B;
thin = MCMCspecs$thin;
burnin = MCMCspecs$burnin;
update = MCMCspecs$update;

# arrays to store MCMC samples
C.samp = array(NA,c(p*(p-1)/2,K,B));
D.samp = array(NA,c(p,K,B));
lam.samp = matrix(NA,K,B);

is=0;
for(t in 1:(B*thin+burnin)) {
    
    for(k in 1:K){
        
    # (1) Update C
    C.dual[[k]] = C.update(n,S.dual[[k]],C.dual[[k]],Dd.dual[[k]],lam.dual[k])
    
    # (2) Update D.diag
    for(stp in 1:5)
    Dd.dual[[k]] = D.update(n,S.dual[[k]],C.dual[[k]],Dd.dual[[k]],D_prior)
    
    # (3) Update lambda
    lam.dual[k] = lam.update(C.dual[[k]],lam_prior)
    
    }

    # (4) Save MCMC sample
    if(t>burnin && (t-burnin)%%thin==0) { 
        is = is+1; # this is the real row number among the B samples.
        C.samp[,,is] <- sapply(C.dual, function(x) x[upper.tri(diag(p))])
        D.samp[,,is] <- abind(Dd.dual, along=2)
        lam.samp[,is] <- lam.dual
    }

    if(t%%update ==0) cat(t,'\n');          

}


### Step 4: Project MCMC results back to conditional cross-covariance in data space 
rho.samp <- -C.samp/(1-C.samp^2)
csd.samp <- apply(D.samp,2:3,function(x) (diag(1/x)%*%matrix(1,p,p)%*%diag(1/x))[upper.tri(diag(p))])

rhosd.arr <- rho.samp*csd.samp
ccc.samp <- abind(sapply(alply(rhosd.arr,3), function(x) x%*%(Phi^2), simplify=FALSE), along=3)

# ccc.samp <- array(NA,c(p*(p-1)/2,T,B))
# for(t in 1:T) ccc.samp[,t,] <- apply(sapply(1:K, function(k) Phi[k,t]^2*csd.samp[,k,]*rho.samp[,k,],simplify='array'),1:2,sum)

return(list(C=C.samp,D=D.samp,lam=lam.samp,ccc=ccc.samp)) 

}

    

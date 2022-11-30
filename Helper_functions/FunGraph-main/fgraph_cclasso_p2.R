fgraph_ccc = function(Y,Phi,MCMCspecs,D_prior,lam_prior,ncor=NULL) {

# Note: functions to conduct functional network inference using lasso on conditional cross-covariance
# This is the parallel compuating version
# Input: Y - (nXpxT) data matrix, p: # of variables; T: # of functionallocation; n: sample size
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
library(doParallel)
library(foreach)
if(!is.null(ncor)) registerDoParallel(ncor)

### Step 1: Apply DWT to project observed functions into basis dual space
n = dim(Y)[1]; p = dim(Y)[2]; T = dim(Y)[3];
K = nrow(Phi); W = ginv(Phi);

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

dir.create('tmp_output')

foreach(k=1:K,.export=c('C.update','Cij.post','D.update','D.post','lam.update','C2rho')) %dopar% {
    
    # arrays to store MCMC samples
    C.samp.k = array(NA,c(p*(p-1)/2,B));
    D.samp.k = array(NA,c(p,B));
    lam.samp.k = rep(NA,B);

    S <- S.dual[[k]]
    C.cur <- C.dual[[k]]
    Dd.cur <- Dd.dual[[k]]
    lam.cur <- lam.dual[k]

    is=0;
    for(t in 1:(B*thin+burnin)) {
    
        # (1) Update C
        C.cur <- C.update(n,S,C.cur,Dd.cur,lam.cur)
    
        # (2) Update D.diag
        for(stp in 1:5) # repeat this step five times for better mixing
        Dd.cur <- D.update(n,S,C.cur,Dd.cur,D_prior)
    
        # (3) Update lambda
        lam.cur <- lam.update(C.cur,lam_prior)
    
        # (4) Save MCMC sample
        if(t>burnin && (t-burnin)%%thin==0) { 
            is = is+1; # this is the real row number among the B samples.
            C.samp.k[,is] <- C.cur[upper.tri(diag(p))]
            D.samp.k[,is] <- Dd.cur
            lam.samp.k[is] <- lam.cur
        }

        if(t%%MCMCspecs$update ==0) cat(t,'\n');     
    }
    save(C.samp.k,D.samp.k,lam.samp.k,file=paste0('tmp_output/mcmc_k',k,'.rda'))
    
    cat('Basis', k, 'is done. \n')
} 


# arrays to pool MCMC samples
C.samp = array(NA,c(p*(p-1)/2,K,B));
D.samp = array(NA,c(p,K,B));
lam.samp = matrix(NA,K,B);
for(k in 1:K){
    load(paste0('tmp_output/mcmc_k',k,'.rda'))
    C.samp[,k,] <- C.samp.k
    D.samp[,k,] <- D.samp.k
    lam.samp[k,] <- lam.samp.k
}
unlink('tmp_output',recursive=TRUE)


### Step 4: Project MCMC results back to conditional cross-covariance in data space 
rho.samp <- -C.samp/(1-C.samp^2)
csd.samp <- apply(D.samp,2:3,function(x) (diag(1/x)%*%matrix(1,p,p)%*%diag(1/x))[upper.tri(diag(p))])
ccc.samp <- array(NA,c(p*(p-1)/2,T,B))
for(t in 1:T) ccc.samp[,t,] <- apply(sapply(1:K, function(k) 
    Phi[k,t]^2*csd.samp[,k,]*rho.samp[,k,],simplify='array'),1:2,sum)

return(list(C=C.samp,D=D.samp,lam=lam.samp,ccc=ccc.samp)) 

}

    

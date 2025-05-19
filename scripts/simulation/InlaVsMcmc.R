library(INLA)
library(mvtnorm)
library(parallel)

load("nbhd.RData")

D = diag(colSums(W))
n = nrow(W)
W_n = sum(W)/2
W_sparse = matrix(0, nrow = W_n, ncol =2)
counter = 1
for(i in 1:(n-1)){
  for(j in (i+1):n){
    if(W[i,j] == 1){
      W_sparse[counter,1] = i
      W_sparse[counter,2] = j
      counter = counter + 1
    }
  }
}
D_sparse = colSums(W)
lambda = eigen(D-W, T)$values

## STAN
sp_d <- list(n = nrow(W),
             N = 2*nrow(W),
             y = numeric(),
             W = W,
             W_n = sum(W)/2,
             W_sparse = W_sparse,
             D_sparse = D_sparse,
             lambda = lambda)

clusters = parallel::detectCores()
NSIM = 10000
indexList = split(1:NSIM, ceiling((1:NSIM)/(NSIM/clusters)))

InlaVsMCMC = function(nsim, W, stanData, eta0 = 0.5, eta1 = -0.25, alpha = 0.2, tau1 = 1, tau2 = 1){
  Q = (alpha * diag(rep(1, nrow(W))) + (1-alpha) * (diag(colSums(W))-W))
  Sigma = solve(Q)
  
  resInla = data.frame(tau = rep(0, length(nsim)),
                       eta0 = rep(0, length(nsim)),
                       eta1 = rep(0, length(nsim)),
                       alpha = rep(0, length(nsim)))
  resStan = resInla
  
  for(i in 1:length(nsim)){
    psi1 = c(rmvnorm(1,sigma = tau1*Sigma))
    mu = c((diag(rep(eta0, nrow(W))) + eta1 * W) %*% psi1)
    psi2 = c(rmvnorm(1, mu ,sigma = tau2*Sigma))
    y1 = rpois(length(psi1), exp(psi1))
    y2 = rpois(length(psi1), exp(psi2))
    
   
    start.time.inla <- Sys.time()
    m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
    mdat = data.frame(y = y1, idx = 1:length(y1))
    l.car <- inla(y ~ f(idx, model = m) - 1, data = mdat, family = "poisson")
    
    mc <- inla.CCAR.model(W = W, phi = l.car$summary.random$idx$mean, alpha.min = 0, alpha.max = 1)
    mdat = data.frame(y = y2 - 1, idx = 1:length(y1))
    c.car <- inla(y2 ~ f(idx, model = mc), data = mdat, family = "poisson")
    end.time.inla <- Sys.time()
    time.taken.inla <- end.time.inla - start.time.inla
    
    start.time.stan <- Sys.time()
    stanData$y = c(y1, y2)
    sp_fit <- stan('BCAR.stan', data = stanData, 
                   iter = 2e4, chains = 1, verbose = FALSE)
    end.time.stan <- Sys.time()
    time.taken.stan <- end.time.stan - start.time.stan
    
    # save results
    resInla$alpha[i] = 1/(1+exp(-c.car$summary.hyperpar$mean[[1]]))
    resInla$tau[i] = exp(c.car$summary.hyperpar$mean[[2]])
    resInla$eta0[i] = c.car$summary.hyperpar$mean[3]
    resInla$eta1[i] = c.car$summary.hyperpar$mean[4]
    resInla$time[i] = time.taken.inla
    
    resStan$alpha[i] = mean(sp_fit@sim$samples[[1]]$alpha2[1e4:2e4])
    resStan$tau[i] = mean(sp_fit@sim$samples[[1]]$tau2[1e4:2e4])
    resStan$eta0[i] = mean(sp_fit@sim$samples[[1]]$eta0[1e4:2e4])
    resStan$eta1[i] = mean(sp_fit@sim$samples[[1]]$eta1[1e4:2e4])
    resStan$time[i] = time.taken.stan
  }
  resInla$estimator = "INLA"
  resStan$estimator = "MCMC"
  
  return(rbind(resInla, resStan))
}

cl = makeCluster(clusters, outfile = "log.txt")
clusterEvalQ(cl, {
  source('parlibs.R')
})
sim = parSapply(cl, indexList, InlaVsMCMC, W, sp_d)
stopCluster(cl)


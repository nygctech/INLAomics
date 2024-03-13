## spotMCCAR implementation spotMCCAR(alpha, phi) -> proper conditonal CAR with multiple processes in conditioning set where some only have spot-to-spot effect
# W: Adjacency SPARSE matrix for spatial effect
# phi: Point estimates of the GMRF in the conditioning set, a n by k matrix
# spotid: vector with positions of assays that only have a spot-to-spot effect
# k: Number of processes(assays)
# alpha.min: Minimum value of the spatial convolution parameter
# alpha.max: Maximum value of the spatial convolution parameter

'inla.rgeneric.mixedMCCAR.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                    "log.prior", "quit"), theta = NULL){
  interpret.theta <- function()
  {
    alpha <- alpha.min + (alpha.max - alpha.min) / (1 + exp(-theta[1L]))
    prec <- exp(theta[2L])
    param = c(alpha, prec)
    return(list(alpha = alpha, prec = prec, param = param))
  }
  
  graph <- function()
  {
    G <- Matrix::Diagonal(nrow(W), 1) + W
    return (G)
  }
  
  #Precision matrix
  Q <- function()
  {
    #Parameters in model scale
    param <- interpret.theta()
    #Precision matrix
    Q <- param$prec*(param$alpha*Matrix::Diagonal(nrow(W), 1) + (1-param$alpha)*(Matrix::Diagonal(nrow(W), apply(W, 1, sum)) - W))
    return (Q)
  }
  
  mu <- function() {
    counter = 0
    loc = numeric(nrow(W))
    for(i in 1:k){
      if(i %in% spotid){
        counter = counter + 1
        loc = loc + (Matrix::Diagonal(nrow(W), theta[2+i])%*%phi[,i])[,1]
      } else {
        loc = loc + (((theta[2+k+i-counter]*W) + Matrix::Diagonal(nrow(W), theta[2+i]))%*%phi[,i])[,1]
      }
    }
    return(loc)
  }
  
  log.norm.const <- function() {
    val <- numeric(0)
    return (val)
  }
  
  log.prior <- function() {
    param <- interpret.theta()
    val <- log(1) + log(param$alpha) + log(1 - param$alpha)
    val <- val + dchisq(param$prec, 1, log=T) + log(param$prec) + sum(dnorm(theta[-(1:2)], 0, 5, log = T))
    return (val)
  }
  
  initial <- function() {
    return (rep(0, 2*(k+1)-length(spotid)))
  }
  
  quit <- function() {
    return (invisible())
  }
  
  val <- do.call(match.arg(cmd), args = list())
  return (val)
}

inla.mixedMCCAR.model <- function(...) {
  INLA::inla.rgeneric.define(inla.rgeneric.mixedMCCAR.model, ...)
}
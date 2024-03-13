## spotMCCAR implementation spotMCCAR(alpha, phi) -> proper conditonal CAR with multiple processes in conditioning set, only spot-to-spot effect
# W: Adjacency SPARSE matrix for spatial effect
# phi: Point estimates of the GMRF in the conditioning set, a n by k matrix
# k: Number of processes(assays)
# alpha.min: Minimum value of the spatial convolution parameter
# alpha.max: Maximum value of the spatial convolution parameter

'inla.rgeneric.spotMCCAR.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                                    "log.prior", "quit"), theta = NULL){
  interpret.theta <- function()
  {
    #alpha <- 1 / (1 + exp(-theta[1L]))
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
    loc = numeric(nrow(W))
    for(i in 1:k){
      loc = loc + (Matrix::Diagonal(nrow(W), theta[2+i])%*%phi[,i])[,1]
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
    return (rep(0, 2+k))
  }
  
  quit <- function() {
    return (invisible())
  }
  
  val <- do.call(match.arg(cmd), args = list())
  return (val)
}

inla.spotMCCAR.model <- function(...) {
  INLA::inla.rgeneric.define(inla.rgeneric.spotMCCAR.model, ...)
}
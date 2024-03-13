## spotCCAR implementation spotCCAR(alpha, phi) -> proper conditonal CAR with single process in conditioning, only with spot-to-spot effect
# W: Adjacency SPARSE matrix for spatial effect
# phi: Point estimates of the GMRF in the conditioning set, a n by k matrix
# alpha.min: Minimum value of the spatial convolution parameter
# alpha.max: Maximum value of the spatial convolution parameter


'inla.rgeneric.spotCCAR.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
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
    return((Matrix::Diagonal(nrow(W), theta[3L])%*%phi)[,1])
  }
  
  log.norm.const <- function() {
    val <- numeric(0)
    return (val)
  }
  
  log.prior <- function() {
    param <- interpret.theta()
    val <- log(1) + log(param$alpha) + log(1 - param$alpha)
    val <- val + dchisq(param$prec, 1, log=T) + log(param$prec) + dnorm(theta[3L], 0, 5, log = T)
    return (val)
  }
  
  initial <- function() {
    return (c(0,0,0))
  }
  
  quit <- function() {
    return (invisible())
  }
  
  val <- do.call(match.arg(cmd), args = list())
  return (val)
}

inla.spotCCAR.model <- function(...) {
  INLA::inla.rgeneric.define(inla.rgeneric.spotCCAR.model, ...)
}
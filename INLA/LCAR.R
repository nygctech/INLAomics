# W: Adjacency SPARSE matrix for spatial effect
# alpha.min: Minimum value of the spatial convolution parameter
# alpha.max: Maximum value of the spatial convolution parameter

'inla.rgeneric.LCAR.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
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
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    val <- numeric(0)
    return (val)
  }
  
  log.prior <- function() {
    param <- interpret.theta()
    
    val <- - theta[1L] - 2 * log(1 + exp(-theta[1L]))
    val <- val + dchisq(param$prec, 1, log=T) + theta[2L]
    return (val)
  }
  
  initial <- function() {
    ## return initial values
    return (c(0,0))
  }
  
  quit <- function() {
    return (invisible())
  }
  
  val <- do.call(match.arg(cmd), args = list())
  return (val)
}

inla.LCAR.model <- function(...) {
  INLA::inla.rgeneric.define(inla.rgeneric.LCAR.model, ...)
}
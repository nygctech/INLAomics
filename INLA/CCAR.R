'inla.rgeneric.CCAR.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
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
    #return((theta[3L]*W)%*%phi)
    return((((theta[4L]*W) + Matrix::Diagonal(nrow(W), theta[3L]))%*%phi)[,1])
  }
  
  log.norm.const <- function() {
    val <- numeric(0)
    return (val)
  }
  
  log.prior <- function() {
    param <- interpret.theta()
    val <- log(1) + log(param$alpha) + log(1 - param$alpha)
    val <- val + dchisq(param$prec, 1, log=T) + log(param$prec) + sum(dnorm(c(theta[3L], theta[4L]), 0, 5, log = T))
    return (val)
  }
  
  initial <- function() {
    return (c(0,0,0,0))
  }
  
  quit <- function() {
    return (invisible())
  }
  
  val <- do.call(match.arg(cmd), args = list())
  return (val)
}

inla.CCAR.model <- function(...) {
  INLA::inla.rgeneric.define(inla.rgeneric.CCAR.model, ...)
}
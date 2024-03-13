'inla.rgeneric.LCAR.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
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
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    ## return the log(normalising constant) for the model
    #param = interpret.theta()
    #
    #val = n * (- 0.5 * log(2*pi) + 0.5 * log(prec.innovation)) +
    #0.5 * log(1.0 - param$alpha^2)
    
    val <- numeric(0)
    return (val)
  }
  
  log.prior <- function() {
    ## return the log-prior for the hyperparameters.
    ## Uniform prior in (alpha.min, alpha.max) on model scale
    param <- interpret.theta()
    
    #val <- log(1) + log(param$alpha) + log(1 - param$alpha)
    val <- - theta[1L] - 2 * log(1 + exp(-theta[1L]))
    # Chisquare for precision
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
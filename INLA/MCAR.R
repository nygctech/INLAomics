## MCAR implementation MCAR(alpha, Lambda) -> PROPER CAR with non-diagonal Lambda and Leroux precision matrix
# W: Adjacency matrix for spatial effect
# k: Number of processes(assays)
# alpha.min: Minimum value of the spatial convolution parameter
# alpha.max: Maximum value of the spatial convolution parameter

'inla.rgeneric.MCAR.model' <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                   "log.prior", "quit"), theta = NULL)
  {
    ## MCAR implementation MCAR(alpha, Lambda) ->
    ## ->  PROPER CAR, common alpha, dense Lambda
    ## k: number of diseases/blocks
    ## W: adjacency matrix
  
    #theta: 1 common correlation parameter alpha,
    #        (k + 1) * k / 2  for lower-tri matrix by col.
    interpret.theta <- function()
    {
      #Function for changing from internal scale to external scale
      # also, build the inverse of the matrix used to model in the external scale
      # the between-disease variability and then this matrix is inverted.
      
      # First parameter is the common autocorrelation parameter
      #alpha <- 1 / (1 + exp(-theta[1L]))
      alpha <- alpha.min + (alpha.max - alpha.min) / (1 + exp(-theta[1L]))
      
      # The next k parameters are the marginal precisions,
      # the other parameters are the correlation parameters ordered by columns.
      mprec <- sapply(theta[as.integer(2:(k+1))], function(x) { exp(x) })
      corre <- sapply(theta[as.integer(-(1:(k+1)))], function(x) {
        (2 * exp(x))/(1 + exp(x)) - 1 })
      
      param <- c(alpha, mprec, corre)
      
      #length non-diagonal elements
      n <- (k - 1) * k / 2
      
      # intial matrix with 1s at the diagonal
      M <- diag(1, k)
      
      #Adding correlation parameters (lower.tri) and (upper.tri)
      M[lower.tri(M)] <- param[k + 2:(n+1)]
      M[upper.tri(M)] <- t(M)[upper.tri(M)]
      
      #Preparing the st. dev matrix
      st.dev <- 1 / sqrt(param[2:(k+1)])
      
      # Matrix of st. dev.
      st.dev.mat <- matrix(st.dev, ncol = 1) %*% matrix(st.dev, nrow = 1)
      
      # Final inversed matrix
      M <- M * st.dev.mat
      
      # Inverting matrix
      #PREC <- MASS::ginv(M)# Generalized inverse
      PREC <- solve(M)
      
      return (list(alpha = alpha, param = param, VACOV = M, PREC = PREC))
    }
    
    
    #Graph of precision function; i.e., a 0/1 representation of precision matrix
    graph <- function()
    {
      PREC <- matrix(1, ncol = k, nrow = k)
      G <- kronecker(PREC, Matrix::Diagonal(nrow(W), 1) + W)
      return (G)
    }
    
    #Precision matrix
    Q <- function()
    {
      #Parameters in model scale
      param <- interpret.theta()
      
      #Precision matrix
      Q <- kronecker(param$PREC,
                     param$alpha*Matrix::Diagonal(nrow(W), 1) + (1-param$alpha)*(Matrix::Diagonal(nrow(W), apply(W, 1, sum)) - W)
      )
      
      return (Q)
    }
    
    #Mean of model
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
      
      # log-Prior for the mixture parameter
      #val <- log(1) + log(param$alpha) + log(1 - param$alpha)
      val <- - theta[1L] - 2 * log(1 + exp(-theta[1L]))
      # Whishart prior for joint matrix of hyperparameters
      val <- val + 
        log(MCMCpack::dwish(W = param$PREC, v =  k, S = diag(rep(1, k)))) +
        sum(theta[as.integer(2:(k + 1))]) +  # This is for precisions
        sum(log(2) + theta[-as.integer(1:(k + 1))] - 2 * log(1 + exp(theta[-as.integer(1:(k + 1))])))  # This is for correlation terms
      
      return (val)
    }
    
    initial <- function() {
      ## return initial values
      
      # The Initial values form a diagonal matrix
      return ( c(0, rep(log(1), k), rep(0, (k * (k - 1) / 2))))
      
    }
    
    quit <- function() {
      return (invisible())
    }
    
    val <- do.call(match.arg(cmd), args = list())
    return (val)
  }

inla.MCAR.model <- function(...) {
  INLA::inla.rgeneric.define(inla.rgeneric.MCAR.model, ...)
}
## MCAR implementation MCAR(alpha, Lambda) -> PROPER CAR with Lambda diagonal and Leroux precision matrix
# W: Adjacency matrix for spatial effect
# k: Number of processes(assays)
# alpha.min: Minimum value of the spatial convolution parameter
# alpha.max: Maximum value of the spatial convolution parameter

'inla.rgeneric.indep.MCAR.model' <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                   "log.prior", "quit"), theta = NULL)
  {
    interpret.theta <- function()
    {
      # Function for changing from internal scale to external scale
      # also, build the precion matrix.
      
      # First parameter is the common autocorrelation parameter
      alpha <- alpha.min + (alpha.max - alpha.min) / (1 + exp(-theta[1L]))
      
      # Next k parameters are the marginal precisions.
      mprec <- sapply(theta[as.integer(2:(k+1))], function(x) { exp(x) })
      
      # Diagonal precion matrix
      PREC <- diag(mprec, k)
      
      return (list(alpha = alpha, mprec = mprec, PREC = PREC))
    }
    
    
    #Graph of precision function; i.e., a 0/1 representation of precision matrix
    graph <- function()
    {
      
      PREC <- diag(1, k)
      G <- kronecker(PREC, Matrix::Diagonal(nrow(W), 1) + W)
      return (G)
    }
    
    #Precision matrix
    Q <- function()
    {
      #Parameters in model scale
      param <- interpret.theta()
      
      # Precision matrix
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
      
      val <- numeric(0)
      return (val)
    }
    
    log.prior <- function() {
      ## return the log-prior for the hyperparameters.
      ## Uniform prior in (alpha.min, alpha.max) on model scale
      param <- interpret.theta()
      val <- - theta[1L] - 2 * log(1 + exp(-theta[1L]))
      val <- val + sum(dchisq(param$mprec, 1, log=T)) - sum(theta[as.integer(2:(k+1))]) / 2 - k * log(2)
      
      return (val)
    }
    
    initial <- function() {
      ## return initial values
      
      #Initial values
      return ( c(0, rep(log(1), k)) )
      
    }
    
    quit <- function() {
      return (invisible())
    }
    
    # FIX for rgeneric to work on R >= 4
    # Provided by E. T. Krainski
    if (as.integer(R.version$major) > 3) {
      if (!length(theta))
        theta = initial()
    } else {
      if (is.null(theta)) {
        theta <- initial()
      }
    }
    
    
    
    val <- do.call(match.arg(cmd), args = list())
    return (val)
  }

##' @rdname indepmcar
##' @param ...  Arguments to be passed to 'inla.rgeneric.define'.
##' @export
##' @usage inla.INDMCAR.model(...)

inla.INDMCAR.model <- function(...) {
  INLA::inla.rgeneric.define(inla.rgeneric.indep.MCAR.model, ...)
}
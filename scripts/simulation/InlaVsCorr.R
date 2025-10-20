library(INLA)
library(mvtnorm)
library(foreach)

load("./scripts/simulation/nbhd.RData")
D <- diag(colSums(W))
Q <- (0.1 * diag(rep(1, nrow(W))) + 0.9 * (D - W))
Sigma <- solve(Q)
n <- nrow(df)

cores <- parallel::detectCores()
NSIM <- 10000
my.cluster <- parallel::makeCluster(cores)
parallel::clusterEvalQ(my.cluster, {
  source("./scripts/simulation/libs.R")
})
doParallel::registerDoParallel(cl = my.cluster)

sim <- foreach(i = 1:NSIM) %dopar% {
  Lambda <- rWishart(1, 10, diag(rep(1, 3)))[, , 1]

  psi_1 <- c(rmvnorm(1, sigma = kronecker(solve(Lambda), Sigma)))
  eta <- matrix(c(
    1, -1.,
    0.5, 0,
    3, 1
  ), nrow = 3, byrow = T)

  A <- cbind((eta[1, 1] * diag(rep(1, n)) + eta[1, 2] * W), (eta[2, 1] * diag(rep(1, n)) + eta[2, 2] * W), (eta[3, 1] * diag(rep(1, n)) + eta[3, 2] * W))
  mu <- c(A %*% psi_1)

  y2 <- rpois(n, exp(1 + mu))
  y1 <- rpois(3 * n, exp(2 + psi_1))
  tryCatch(
    {
      k <- 3
      X <- kronecker(diag(rep(1, 3)), matrix(1, nrow = n))
      mdat <- data.frame(y = y1, idx = 1:(3 * n), x1 = X[, 1], x2 = X[, 2], x3 = X[, 3])
      m <- inla.MCAR.model(W = W, k = 3, alpha.min = 0, alpha.max = 1)
      m.car <- inla(y ~ x1 + x2 + x3 - 1 + f(idx, model = m), data = mdat, family = "poisson")

      m <- inla.MCCAR.model(
        W = W, phi = matrix(m.car$summary.random$idx$mean, ncol = k),
        k = k, alpha.min = 0, alpha.max = 1
      )
      mdat <- data.frame(y = y2, idx = 1:n)
      mc.car <- inla(y ~ f(idx, model = m),
        data = mdat, family = "poisson"
      )

      data.frame(
        "eta0" = mc.car$internal.summary.hyperpar$mean[c(3, 4, 5)],
        "eta1" = mc.car$internal.summary.hyperpar$mean[c(2 + k + 1, 2 + k + 2, 2 + k + 3)],
        "corr" = c(cor(y2, y1[1:n]), cor(y2, y1[(n + 1):(2 * n)]), cor(y2, y1[(2 * n + 1):(3 * n)]))
      )
    },
    error = function(cond) {
      data.frame("eta0" = NA, "eta1" = NA, "corr" = NA)
    }
  )
}
parallel::stopCluster(cl = my.cluster)

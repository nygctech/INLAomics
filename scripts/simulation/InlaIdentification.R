library(INLA)
library(parallel)
library(Rcpp)
sourceCpp("simLargeNormal.cpp") # Path to the C++ file

source("./scripts/utils.R")
source("./scripts/Highplex/helpers.R")
source("./INLA/LCAR.R")
source("./INLA/MCAR.R")
source("./INLA/MCCAR.R")


load("./scripts/simulation/nbhd.RData")
coords <- cbind(cn87$x, cn87$y)

createGrid <- function(coords, reps = 0) {
  if (reps > 0) {
    base <- coords
    for (i in 1:reps) {
      shift <- max(coords[, 1]) - min(coords[, 1]) + 1
      temp <- base
      temp[, 1] <- temp[, 1] + shift
      coords <- rbind(coords, temp)
    }
  }
  return(coords)
}

###############
#### G = 2 ####
###############

InlaVShetero <- function(nsim, W, D, params) {
  n <- nrow(W)
  Lambda <- matrix(c(
    1 / params$tau11,
    params$rho / sqrt(params$tau11 * params$tau12),
    params$rho / sqrt(params$tau11 * params$tau12),
    1 / params$tau12
  ), 2, byrow = T)
  Q <- (params$pi1 * diag(rep(1, nrow(W))) + (1 - params$pi1) * (D - W))

  res <- data.frame(
    tau11 = numeric(length(nsim)),
    tau12 = numeric(length(nsim)),
    tau2 = numeric(length(nsim)),
    rho = numeric(length(nsim)),
    pi1 = numeric(length(nsim)),
    pi2 = numeric(length(nsim)),
    eta01 = numeric(length(nsim)),
    eta02 = numeric(length(nsim)),
    eta11 = numeric(length(nsim)),
    eta12 = numeric(length(nsim)),
    time = numeric(length(nsim))
  )

  for (i in 1:length(nsim)) {
    start <- proc.time()
    Psi <- sample_from_inverse_kron(solve(Lambda), Q)
    rna <- rpois(2 * n, exp(Psi))
    mu <- c((diag(params$eta01, n, n) + params$eta11 * W) %*% Psi[1:n])
    mu <- mu + c((diag(params$eta02, n, n) + params$eta12 * W) %*% Psi[(n + 1):(2 * n)])
    Q <- (params$pi2 * diag(rep(1, nrow(W))) + (1 - params$pi2) * (D - W))
    psi <- c(rmvnorm(1, mean = mu, sigma = solve(Q) / 2))
    prot <- rpois(n, exp(psi))
    df <- data.frame(
      prot = prot,
      rna1 = rna[1:n],
      rna2 = rna[(n + 1):(2 * n)],
      size_prot = rep(1, n),
      size_rna = rep(1, n),
      idx = 1:n
    )
    k <- 2
    X <- kronecker(diag(rep(1, k)), rep(1, nrow(df)))
    mdat <- data.frame(
      "rna" = unname(unlist(as.vector(df[, names(df) %in% c("rna1", "rna2")]))),
      "idx" = 1:(k * nrow(df)),
      "size" = rep(df$size_rna, k)
    )
    mdat <- cbind(mdat, X)
    names(mdat)[4:ncol(mdat)] <- paste("x", 1:k, sep = "")
    rnaform <- as.formula(paste("rna ~", paste(names(mdat)[4:(ncol(mdat))], collapse = "+"), "-1"))

    m <- inla.MCAR.model(W = W, k = k, alpha.min = 0, alpha.max = 1)
    m.car <- inla(update(rnaform, . ~ . + f(idx, model = m)),
      data = mdat, family = "poisson", offset = log(size)
    )

    protform <- as.formula(paste(paste("prot", "~"), "1"))
    m <- inla.MCCAR.model(
      W = W, phi = matrix(m.car$summary.random$idx$mean, ncol = k),
      k = k, alpha.min = 0, alpha.max = 1
    )

    mc.car <- inla(update(protform, . ~ . + f(idx, model = m)),
      data = df, family = "poisson",
      offset = log(size_prot)
    )

    end <- proc.time() - start
    res$tau11[i] <- exp(m.car$summary.hyperpar$mean[2])
    res$tau12[i] <- exp(m.car$summary.hyperpar$mean[3])
    res$rho[i] <- 2 * exp(m.car$summary.hyperpar$mean[4]) / (1 + exp(m.car$summary.hyperpar$mean[4])) - 1
    res$pi1[i] <- 1 - plogis(m.car$summary.hyperpar$mean[1])
    res$tau2[i] <- exp(mc.car$summary.hyperpar$mean[2])
    res$pi2[i] <- 1 - plogis(mc.car$summary.hyperpar$mean[1])
    res$eta01[i] <- mc.car$summary.hyperpar$mean[3]
    res$eta02[i] <- mc.car$summary.hyperpar$mean[4]
    res$eta11[i] <- mc.car$summary.hyperpar$mean[5]
    res$eta12[i] <- mc.car$summary.hyperpar$mean[6]
    res$time[i] <- as.numeric(end[3])
  }

  return(res)
}

simgrid <- createGrid(coords, 2) # n = 3*nrow(coords) = 1023
n <- nrow(simgrid)
W <- matrix(0, nrow(simgrid), nrow(simgrid))
for (i in 1:nrow(W)) {
  for (j in i:nrow(W)) {
    x <- simgrid[i, ]
    if (crossprod(simgrid[i, ] - simgrid[j, ]) == 1) {
      W[i, j] <- 1
      W[j, i] <- 1
    }
  }
}
D <- diag(colSums(W))

clusters <- 5
NSIM <- 20 * clusters
indexList <- split(1:NSIM, ceiling((1:NSIM) / (NSIM / clusters)))
params <- list(
  tau11 = 2, tau12 = 0.5, tau2 = 0.5,
  rho = 0.8, pi1 = 0.9, pi2 = 0.9, eta01 = 1.,
  eta02 = -0.8, eta11 = 0.25, eta12 = -0.5
)

cl <- makeCluster(clusters)
clusterEvalQ(cl, {
  source("./scripts/simulation/libs.R")
})
a <- parSapply(cl, indexList, InlaVShetero, W, D, params)
stopCluster(cl)


###############
#### G = 3 ####
###############

InlaVShetero <- function(nsim, W, D, params) {
  n <- nrow(W)

  Lambda <- matrix(
    c(
      1 / params$tau11, params$rho12 / sqrt(params$tau11 * params$tau12), params$rho13 / sqrt(params$tau11 * params$tau13),
      params$rho12 / sqrt(params$tau11 * params$tau12), 1 / params$tau12, params$rho23 / sqrt(params$tau12 * params$tau13),
      params$rho13 / sqrt(params$tau11 * params$tau13), params$rho23 / sqrt(params$tau12 * params$tau13), 1 / params$tau13
    ),
    3, 3,
    byrow = T
  )
  k <- ncol(Lambda)

  res <- data.frame(
    tau11 = numeric(length(nsim)),
    tau12 = numeric(length(nsim)),
    tau13 = numeric(length(nsim)),
    tau2 = numeric(length(nsim)),
    rho12 = numeric(length(nsim)),
    rho13 = numeric(length(nsim)),
    rho23 = numeric(length(nsim)),
    pi1 = numeric(length(nsim)),
    pi2 = numeric(length(nsim)),
    eta01 = numeric(length(nsim)),
    eta02 = numeric(length(nsim)),
    eta03 = numeric(length(nsim)),
    eta11 = numeric(length(nsim)),
    eta12 = numeric(length(nsim)),
    eta13 = numeric(length(nsim)),
    time = numeric(length(nsim))
  )
  Q1 <- (params$pi1 * diag(rep(1, nrow(W))) + (1 - params$pi1) * (D - W))
  Q2 <- (params$pi2 * diag(rep(1, nrow(W))) + (1 - params$pi2) * (D - W))

  for (i in 1:length(nsim)) {
    start <- proc.time()
    Q <- (params$pi1 * diag(rep(1, nrow(W))) + (1 - params$pi1) * (D - W))
    Psi <- sample_from_inverse_kron(solve(Lambda), Q1)
    rna <- rpois(k * n, exp(Psi))
    mu <- c((diag(params$eta01, n, n) + params$eta11 * W) %*% Psi[1:n])
    mu <- mu + c((diag(params$eta02, n, n) + params$eta12 * W) %*% Psi[(n + 1):(2 * n)])
    mu <- mu + c((diag(params$eta03, n, n) + params$eta13 * W) %*% Psi[(2 * n + 1):(3 * n)])
    psi <- c(rmvnorm(1, mean = mu, sigma = solve(Q2) * params$tau2))
    prot <- rpois(n, exp(psi))
    df <- data.frame(
      prot = prot,
      rna1 = rna[1:n],
      rna2 = rna[(n + 1):(2 * n)],
      rna3 = rna[(2 * n + 1):(3 * n)],
      size_prot = rep(1, n),
      size_rna = rep(1, n),
      idx = 1:n
    )
    X <- kronecker(diag(rep(1, k)), rep(1, nrow(df)))
    mdat <- data.frame(
      "rna" = unname(unlist(as.vector(df[, names(df) %in% c("rna1", "rna2", "rna3")]))),
      "idx" = 1:(k * nrow(df)),
      "size" = rep(df$size_rna, k)
    )
    mdat <- cbind(mdat, X)
    names(mdat)[4:ncol(mdat)] <- paste("x", 1:k, sep = "")
    rnaform <- as.formula(paste("rna ~", paste(names(mdat)[4:(ncol(mdat))], collapse = "+"), "-1"))

    m <- inla.MCAR.model(W = W, k = k, alpha.min = 0, alpha.max = 1)
    m.car <- inla(update(rnaform, . ~ . + f(idx, model = m)),
      data = mdat, family = "poisson", offset = log(size)
    )

    protform <- as.formula(paste(paste("prot", "~"), "1"))
    m <- inla.MCCAR.model(
      W = W, phi = matrix(m.car$summary.random$idx$mean, ncol = k),
      k = k, alpha.min = 0, alpha.max = 1
    )

    mc.car <- inla(update(protform, . ~ . + f(idx, model = m)),
      data = df, family = "poisson",
      offset = log(size_prot)
    )

    end <- proc.time() - start

    res$tau11[i] <- exp(m.car$summary.hyperpar$mean[2])
    res$tau12[i] <- exp(m.car$summary.hyperpar$mean[3])
    res$tau13[i] <- exp(m.car$summary.hyperpar$mean[4])
    res$rho12[i] <- 2 * exp(m.car$summary.hyperpar$mean[5]) / (1 + exp(m.car$summary.hyperpar$mean[5])) - 1
    res$rho13[i] <- 2 * exp(m.car$summary.hyperpar$mean[6]) / (1 + exp(m.car$summary.hyperpar$mean[6])) - 1
    res$rho23[i] <- 2 * exp(m.car$summary.hyperpar$mean[7]) / (1 + exp(m.car$summary.hyperpar$mean[7])) - 1
    res$pi1[i] <- 1 - plogis(m.car$summary.hyperpar$mean[1])
    res$tau2[i] <- exp(mc.car$summary.hyperpar$mean[2])
    res$pi2[i] <- 1 - plogis(mc.car$summary.hyperpar$mean[1])
    res$eta01[i] <- mc.car$summary.hyperpar$mean[3]
    res$eta02[i] <- mc.car$summary.hyperpar$mean[4]
    res$eta03[i] <- mc.car$summary.hyperpar$mean[5]
    res$eta11[i] <- mc.car$summary.hyperpar$mean[6]
    res$eta12[i] <- mc.car$summary.hyperpar$mean[7]
    res$eta13[i] <- mc.car$summary.hyperpar$mean[8]
    res$time[i] <- as.numeric(end[3])
  }
  return(res)
}

simgrid <- createGrid(coords, 2) # n = 3*nrow(coords) = 1023
n <- nrow(simgrid)
W <- matrix(0, nrow(simgrid), nrow(simgrid))
for (i in 1:nrow(W)) {
  for (j in i:nrow(W)) {
    x <- simgrid[i, ]
    if (crossprod(simgrid[i, ] - simgrid[j, ]) == 1) {
      W[i, j] <- 1
      W[j, i] <- 1
    }
  }
}
D <- diag(colSums(W))

clusters <- 5
NSIM <- 20 * clusters
indexList <- split(1:NSIM, ceiling((1:NSIM) / (NSIM / clusters)))

params <- list(
  tau11 = 2, tau12 = 0.5, tau13 = 0.7, tau2 = 0.5,
  rho12 = 0.9, rho13 = -0.9, rho23 = -0.8,
  pi1 = 0.9, pi2 = 0.9,
  eta01 = 1., eta02 = 0.5, eta03 = -0.4, eta11 = 0.25, eta12 = 0.2, eta13 = -0.2
)

cl <- makeCluster(clusters)
clusterEvalQ(cl, {
  source("./scripts/simulation/libs.R")
})
a <- parSapply(cl, indexList, InlaVShetero, W, D, params)
stopCluster(cl)

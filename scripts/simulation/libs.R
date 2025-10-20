library(INLA)
library(mvtnorm)
library(rstan)
library(Rcpp)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
inla.setOption(inla.timeout = 5 * 60)

source("./simulation/CCAR.R")
source("./INLA/MCAR.R")
source("./INLA/LCAR.R")
source("./INLA/MCCAR.R")
source("./INLA/CCAR.R")

sourceCpp("./simulation/simLargeNormal.cpp")

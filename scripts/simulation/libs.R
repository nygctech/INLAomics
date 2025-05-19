library(INLA)
library(mvtnorm)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("CCAR.R")

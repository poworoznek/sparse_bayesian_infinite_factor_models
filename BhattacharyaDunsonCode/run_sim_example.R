## MCMC for Latent Factor models as in Bhattacharya and Dunson (2011)

# Load in the data (or simulate)
rm(list = ls())
library(mvtnorm)
library(tidyverse)
library(abind)
library(MSFA)
setwd("~/documents/factorR/")
source("./simulate_data_fxns.R")
source("./BhattacharyaDunsonCode/post_process_functions.R")
source("./BhattacharyaDunsonCode/mcmc.R")

sigmasq = 1
p = 1000
n = 400
k = 20

output = simulate_x(n, p, k, sigmasq, simulate_lambda, return_lambda = TRUE)
Y = output$x
#Y = matrix(rnorm(10000), ncol = 50)
n = nrow(Y)
p = ncol(Y)

## Set MCMC constants
iter = 3000
burn = 500

rm(list = c("k", "n", "p"))

# Run the Gibbs Sampler from Bhattacharya and Dunson (2011)
stores = mcmc(Y, iter, burn)

# Post processing of the samples
Omega = get_omega(stores)
round(Omega - output$sigma, 1)
round(Omega, 1)
round(output$sigma, 1)

# Using the Orthogonal Procrustes approach to recover the Lambda factor loadings matrix
# Just using the De Vito and Parmigiani function at the moment
Lambda = get_lambda(stores)
head(round(Lambda, 1), 10)
head(round(output$lambda, 1), 10)


# These are... okay? Not great for sure, doesn't look as positive as the Bayesian MSFA paper
lambda_cor(Lambda, output$lambda)

Sigma = get_sigma(stores)

tcrossprod(Lambda)
tcrossprod(output$lambda)
tcrossprod(Lambda) + Sigma
Omega

# Compare to the Lambda derived directly from the Bayesian MSFA package


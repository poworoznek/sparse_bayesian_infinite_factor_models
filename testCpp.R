library(tidyverse)
library(mvtnorm)
library(abind)
source("simulate_data_fxns.R")

sigmasq = 1
p = 100
n = 100
k = 10

outp = simulate_x(n, p, k, sigmasq, simulate_lambda, 
                    return_lambda = TRUE)
Y = outp$x

prop = 1
epsilon = 1e-3
nrun=20000
burn=1
thin = 1
output = "numFactors"

p = ncol(Y)
n = nrow(Y)

as = 1                          # gamma hyperparameters for residual precision
bs = 0.3                        
df = 3                          # gamma hyperparameters for t_{ij}
ad1 = 2.1
bd1 = 1                         # gamma hyperparameters for delta_1
ad2 = 3.1
bd2 = 1                         # gamma hyperparameters delta_h, h >= 2
adf = 1
bdf = 1 
b0 = 1
b1 = 0.0005

k = floor(log(p)*3)

sp = floor((nrun - burn)/thin)        # number of posterior samples

VY= apply(Y, 2, var)                  # explicitly preserve scale
scaleMat = sqrt((VY) %*% t(VY))
Y = scale(Y)
num = 0

# --- Initial values --- #
ps = rgamma(p, as, bs)
Sigma = diag(1/ps)                             # Sigma = diagonal residual covariance
Lambda = matrix(1, nrow = p, ncol = k)
ta = matrix(rnorm(n*k), nrow = n, ncol = k)    # factor loadings & latent factors
meta = matrix(0,nrow = n, ncol = k)
veta = diag(k)                                 # latent factor distribution = standard normal

psijh = matrix(rgamma(p*k, df/2, df/2), nrow = p, ncol = k)     # local shrinkage coefficients
theta = c(1 / rgamma(1,ad1,bd1), 1 / rgamma(k-1,ad2,bd2))       # gobal shrinkage coefficients multilpliers
tauh = 1 / cumprod(theta)                                       # global shrinkage coefficients
Plam = t(t(psijh) * (tauh))                                     # precision of loadings rows

Rcpp::sourceCpp('MGSPc.cpp')

out = MGSPsamp(p, n, k, as, bs, df, ad1, bd1, ad2, bd2, adf, bdf, 
               b0, b1, sp, nrun, burn, thin, prop, epsilon, ps, Sigma, Lambda, 
               meta, veta, psijh, theta, tauh, Plam, Y, scaleMat, output) 

source("safefact.R")
out2 = safefact(Y = outp$x, nrun = nrun, burn = burn, output=output)

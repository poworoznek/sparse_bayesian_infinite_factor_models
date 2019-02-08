# test rotations
library(mvtnorm)
source("fixedfact.R")
source("mcrotfact.R")
source("spclone.R")

set.seed(13)

vals = rep(c(rep(1, 15), rep(0, 97)), 10)[1:1000]

lam = matrix(vals, nrow = 100)
lam[,10] = runif(100) * (-1)^(round(runif(100)))

image(t(lam))

cov = lam %*% t(lam) + diag(100)

dat = rmvnorm(1000, mean = rep(0, 100), sigma = cov)

fixedfact(Y = dat, nrun = 10000, burn = 2000, k = 10, 
          output = "factSamples", factfilename = "Lambda2.rds")

mcrotated = mcrotfact(lambdafile = "Lambda2.rds", method = "BADFM", maxiter = 100, ncores = 4, tol = 1e-5)

image(t(mcrotated$mean))

fact = readRDS("Lambda2.rds")
fmean = Reduce("+", fact)

image(t(fmean))

sprotated = spclone(lambdafile = "Lambda2.rds", method = "BADFM", maxiter = 100, ncores = 4, tol = 1e-5)

image(t(sprotated$mean))

image(t(mcrotated$mean))
image(t(sprotated$mean))

norm(mcrotated$mean - sprotated$mean, type = "F")

vmrotated = mcrotfact(lambdafile = "Lambda2.rds", method = "varimax", maxiter = 100, ncores = 4, tol = 1e-5)

image(t(vmrotated$mean))
image(t(fmean))

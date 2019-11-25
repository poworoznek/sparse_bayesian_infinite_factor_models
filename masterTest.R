library(tidyverse)
Rcpp::sourceCpp('samplerBits.cpp')
Rcpp::sourceCpp('msf.cpp')
source("plotmat.R")
source("safefact.R")
source('lmean.R')
source('matchsignfact.R')

library(mvtnorm)
library(abind)
source("simulate_data_fxns.R")

set.seed(10)
sigmasq = 1
p = 10
n = 100
k = 3

true = simulate_x(n, p, k, sigmasq, simulate_lambda, return_lambda = T)
true$Sigma = with(true, sigma - lambda %*% t(lambda))

X = true$x

#################################
########## MGSP LINEAR ##########
#################################

source('linearMGSP.R')

linmgsp = linearMGSP(X, nrun = 100000, burn = 10000, adapt = F, kinit = 20)

plotmat((true$sigma))

plotmat((linmgsp$covMean))

plotmat(cov(X))

rot = lapply(linmgsp$lambdaSamps, varimax, FALSE, 1e-15)
rot = lapply(rot, `[[`, 1)
al = lapply(rot, matchsignfact, rot[[2500]])

plotmat(al[[1]])
plotmat(al[[10000]])


plotmat(t(t(summat(al)[,c(4, 3, 2, 1)]) * c(-1, -1, 1, 1)))
plotmat(varimax(true$lambda, FALSE, 1e-15)[[1]])

###############################
########## DL LINEAR ##########
###############################

source("linearDL.R")
lindl = linearDL(X, nrun = 10000, burn = 5000, kinit = 10)

plotmat(true$sigma)
plotmat(lindl$covMean)

rot2 = lapply(lindl$lambdaSamps, varimax)
rot2 = lapply(rot2, `[[`, 1)
al2 = lapply(rot2, matchsignfact, rot2[[2500]])

plotmat(al2[[1]])
plotmat(al2[[100]])

plotmat(lmean(al))
plotmat(lmean(al2))
plotmat(true$lambda)

plotmat(linmgsp$lambdaSamps[[1000]])
plotmat(lindl$lambdaSamps[[1000]])


####################################
########## DL INTERACTION ##########
####################################

source("interactionDL.R")

####### Simulate #######

set.seed(1)
k0 = 5
p = 20
n = 500

lambda = matrix(rnorm(p*k0, 0, 0.01), ncol = k0)
lambda[sample.int(p, 40, replace = TRUE) +
         p*(sample.int(k0, 40, replace = TRUE)-1)] = rnorm(40, 0, 1)
lambda[1:7, 1] = rnorm(7, 2, 0.5)
lambda[8:14, 2] = rnorm(7, -2, 0.5)
lambda[15:20, 3] = rnorm(6, 2, 0.5)
lambda[,4] = rnorm(p, 0, 0.5)
lambda[,5] = rnorm(p, 0, 0.5)
plotmat(varimax(lambda)[[1]])

X = matrix(rnorm(n*k0),n,k0)%*%t(lambda) + bayesSurv::rMVNorm(n,Sigma = 3*diag(p))

beta_true = numeric(p); beta_true[c(1,3,6,8,10,11)] =c(1,1,0.5,-1,-2,-0.5)
Omega_true = matrix(0,p,p)
Omega_true[1,2] = 1; Omega_true[5,2] = -1; Omega_true[10,8] = 1; 
Omega_true[11,5] = -2; Omega_true[1,1] = 0.5; 
Omega_true[2,3] = 0.5; 
Omega_true = Omega_true + t(Omega_true)
y = X%*%beta_true + diag(X%*%Omega_true%*%t(X)) +  rnorm(n,0.5)

###### RUN ######

intdl = interactionDL(y, X, 10000, 5000, kinit = 10)

plotmat(intdl$covMean)
plotmat(tcrossprod(lambda) + diag(rep(3, 20)))

rot3 = lapply(intdl$lambdaSamps, varimax, TRUE, 1e-15)
rot3 = lapply(rot3, `[[`, 1)
al3 = lapply(rot3, msf, rot3[[2500]])

plotmat(al3[[1]])
plotmat(al3[[100]])
plotmat(rot3[[1]])
plotmat(rot3[[100]])

plotmat(lmean(al))
plotmat(lmean(al3))
plotmat(varimax(lambda)[[1]])


######################################
########## MGSP INTERACTION ##########
######################################

source("interactionMGSP.R")

intmgsp = interactionMGSP(y, X, 10000, 5000, kinit = 10)

plotmat(intmgsp$covMean)
plotmat(tcrossprod(lambda) + diag(rep(3, 20)))

rot4 = lapply(intmgsp$lambdaSamps, varimax)
rot4 = lapply(rot4, `[[`, 1)
al4 = lapply(rot4, msf, rot4[[2500]])

plotmat(al4[[1]])
plotmat(al4[[100]])

plotmat(lmean(al))
plotmat(lmean(al4))
plotmat(varimax(lambda)[[1]])

plotmat(intmgsp$lambdaSamps[[1000]])
plotmat(intdl$lambdaSamps[[1000]])

plotmat(amean(intmgsp$interactionSamps))
plotmat(amean)




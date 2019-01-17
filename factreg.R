# regression with factorized covariance representation
# see fastfact for full prior outline
# priors are applied to factor loadings, not coefficients

# ARGUMENTS: Y: Matrix of responses;
#            X: Matrix of predictors;
#            intercept: whether or not to add intercept column to X
#            nrun: number of mcmc iterations;
#            burn: number of burn in iterations;
#            thin: thinning interval;
#            prop: proportion of elements in each column less than epsilon in magnitude cutoff;
#            epsilon: tolerance;
#            kinit: initial value for the number of factors;
#            covfilename: optional filename for covariance matrix samples;
#            factfilename: optional filename for factor matrix samples;
#            sigfilename: optional filename for sigma matrix samples;

source("fastfact.R")
source("SMWreg.R")

factreg = function(Y, X, intercept = FALSE, nrun, burn, thin = 1, prop = 1, epsilon = 1e-3,
                   kinit = NULL, covfilename = "OmegaReg.rds", factfilename = "LambdaReg.rds", 
                   sigfilename = "SigmaReg.rds"){
  
  if(nrow(Y) != nrow(X)) stop("Mismatching response and predictor matrix dimensions")
  
  responseLength = ncol(Y)
  if(intercept) X = cbind(1, X)
  predictorLength = ncol(X)
  data = cbind(Y, X)
  
  cat("---Sampling---\n")
  sample = fact2(Y = data, prop = prop, epsilon = epsilon, 
                 kinit = kinit, nrun = nrun, burn = burn, 
                 thin = thin, output = c("covSamples","factSamples", "sigSamples"),
                 covfilename = covfilename, factfilename = factfilename, 
                 sigfilename = sigfilename)
  
  cat("---Coefficient Computation---\n")
  betaSample = SMWreg(response = 1:responseLength, covfilename, factfilename, sigfilename)
  
  return(betaSample)
}



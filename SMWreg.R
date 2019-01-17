# underlying regression function
# develops coefficient mcmc samples
# relies on Lambda and Sigma RDS objects
# based on the Sherman–Morrison–Woodbury formula

# ARGUMENTS: response: column indices of response variables;
#            covfile: file path to covariance matrix sample array;
#            lambdafile: file path to factor matrix sample list;
#            sigmafile: file path to sigma matrix sample array;

SMWreg = function(response, covfile, lambdafile, sigmafile){
  
  lambda = readRDS(lambdafile) # take in output from fact2
  sigma = readRDS(sigmafile)
  cov = readRDS(covfile)
  n = length(lambda)           # same chain length
  
  omegaXXinv = lapply(1:n, function(ind){
    lam = lambda[[ind]][-response,]      # lambda has variable column count
    sig = sigma[-response,-response,ind] # sigma is always p x p
    
    sigsolve = solve(sig)                   # avoid inverting the sig matrix 4 times
    omega = sigsolve - sigsolve %*% lam %*% # implement SMW identity
      solve(diag(ncol(lam)) + t(lam) %*% sigsolve %*% lam) %*%
      t(lam) %*% sigsolve
    return(omega)
  })
  
  omegaZX = lapply(1:n, function(ind){
    omega = cov[-response, response, ind]
    return(omega)
  })
  
  beta = lapply(1:n, function(ind){
    betasamp = omegaXXinv[[ind]] %*% omegaZX[[ind]] # compute regression coefficients for each sample
    return(betasamp)
  })
  
  return(beta)
}


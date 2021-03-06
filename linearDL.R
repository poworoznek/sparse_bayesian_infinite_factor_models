# Gibbs sampler for factorized covariance estimation
# using dirichlet-laplace prior on factor loadings from Battacharya Dunson (2014)
# same eta, lambda, and sigma updates as MGSP linear

# ARGUMENTS: Y: Data matrix (n x p); 
#            nrun: number of iterations;
#            burn: burn-in period;
#            thin: thinning interval;
#            prop: proportion of elements in each column less than epsilon in magnitude cutoff;
#            epsilon: tolerance;
#            kinit: initial value for the number of factors;
#            output: output type, a vector including some of:
#            c("covMean", "covSamples", "factSamples", "sigSamples", "numFactors");
#            covfilename: optional filename for covariance matrix samples;
#            factfilename: optional filename for factor matrix samples;
#            sigfilename: optional filename for sigma matrix samples;

linearDL = function(Y, nrun, burn, thin = 1, prop = 1, epsilon = 1e-3,
                      kinit = NULL, output = c("covMean", "covSamples", 
                                               "factSamples", "sigSamples", 
                                               "numFactors"), 
                      covfilename = "Omega.rds", factfilename = "Lambda.rds", 
                      sigfilename = "Sigma.rds", verbose = TRUE,
                      dump = FALSE, buffer = 10000){
  
  cm = any(output %in% "covMean")
  cs = any(output %in% "covSamples")
  fs = any(output %in% "factSamples")
  ss = any(output %in% "sigSamples")
  nf = any(output %in% "numFactors")
  
  p = ncol(Y)
  n = nrow(Y)

  a = 1/2
  as = 1
  bs = 0.3
  
  if(is.null(kinit)) kinit = floor(log(p)*3)
  
  sp = floor((nrun - burn)/thin)        # number of posterior samples
  
  VY= apply(Y, 2, var)                  # explicitly preserve scale
  scaleMat = sqrt((VY) %*% t(VY))
  Y = scale(Y)
  k=kinit                               # no. of factors to start with
  
  # --- Initial values --- #
  ps = rgamma(p,as,bs)
  lambda = matrix(0,p,k)
  eta = matrix(rnorm(n*k),n,k)
  
  tau = rgamma(p,p*a,1/2)
  psi = matrix(rexp(p*k,1/2),p,k)
  phi = matrix(0,p,k)
  for(j in 1:p){
    gam = rgamma(k, a)
    phi[j,] = gam /sum(gam)
  }
  Plam = phi*(phi^2)*matrix(rep(tau^2,k),p,k,byrow=F)
  
  if(cm) COVMEAN = matrix(0, nrow = p, ncol = p)
  if(cs) OMEGA = array(dim = c(p, p, sp))
  if(fs) LAMBDA = list()
  if(ss) SIGMA = array(dim = c(p, sp))
  if(nf) K = rep(NA, sp)
  ind = 1
  
  if(verbose) {
    pb = txtProgressBar(style = 3)
    at = ceiling(nrun/100)
  }
  
  for(i in 1:nrun){
    eta = eta_lin(lambda, ps, k, n, Y)
    Plam = plm_dl(psi, phi, tau)
    lambda = lam_lin(eta, Plam, ps, k, p, Y)
    psi = psi_dl(lambda, phi, tau)
    tau = tau_dl(lambda, phi, k, p)
    phi = phi_dl(lambda, a, k, p)
    ps = sig_lin(lambda, eta, k, p, n, Y, as, bs)
    
    
    if((i %% thin == 0) & (i > burn)) {
      if(cm | cs) Omega = (tcrossprod(lambda) + diag(1/c(ps))) * scaleMat
      if(cm) COVMEAN = COVMEAN + Omega / sp
      if(cs) OMEGA[,,ind] = Omega
      if(fs) LAMBDA[[ind]] = lambda
      if(ss) SIGMA[,ind] = 1/ps
      if(nf) K[ind] = k
      ind = ind + 1
    }
    
    if(verbose & (i %% at == 0)) setTxtProgressBar(pb, i / nrun)
    
    if(dump & (i %% buffer == 0)){
      if(cs) saveRDS(OMEGA, "Omega.rds", compress = FALSE)
      if(fs) saveRDS(LAMBDA, "Lambda.rds", compress = FALSE)
      if(ss) saveRDS(SIGMA, "Sigma.rds", compress = FALSE)
    }
  }
  
  if(verbose) close(pb)
  
  out = list()
  if(cm) out = c(out, list(covMean = COVMEAN))
  if(cs) out = c(out, list(omegaSamps = OMEGA))
  if(fs) out = c(out, list(lambdaSamps = LAMBDA))
  if(ss) out = c(out, list(sigmaSamps = SIGMA))
  if(nf) out = c(out, list(numFacts = K))
  
  return(out)
}

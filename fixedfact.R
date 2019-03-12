# Gibbs sampler for factorized covariance estimation
# using mgps prior on factor loadings from Battacharya Dunson (2011)
# prior edits informed by Durante (2017)
# version for fixed number of factors

# ARGUMENTS: Y: Data matrix (n x p); 
#            nrun: number of iterations;
#            burn: burn-in period;
#            thin: thinning interval;
#            k: value for the number of factors (Must be passed);
#            output: output type, a vector including some of:
#            c("covMean", "covSamples", "factSamples", "sigSamples");
#            covfilename: optional filename for covariance matrix samples;
#            factfilename: optional filename for factor matrix samples;
#            sigfilename: optional filename for sigma matrix samples;

fixedfact = function(Y, nrun, burn, thin = 1, 
                 k, output = "covMean", 
                 covfilename = "Omega.rds", factfilename = "Lambda.rds", 
                 sigfilename = "Sigma.rds"){
  
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
  
  # --- Allocate output object memory --- #
  if(any(output %in% "covMean")) COVMEAN = matrix(0, nrow = p, ncol = p)
  if(any(output %in% "covSamples")) OMEGA = array(dim = c(p, p, sp))
  if(any(output %in% "factSamples")) LAMBDA = list()
  if(any(output %in% "sigSamples")) SIGMA = array(dim = c(p, p, sp))
  ind = 1
  
  #------start gibbs sampling-----#
  
  cat("Start\n")
  
  for(i in 1:nrun) {
    
    # -- Update eta -- #
    Lmsg = Lambda * ps
    Veta1 = diag(k) + t(Lmsg) %*% Lambda
    Tmat = chol(Veta1)
    R = qr.R(qr(Tmat))
    S = solve(R)
    Veta = S %*% t(S)                                               # Veta = inv(Veta1)
    Meta = Y %*% Lmsg %*% Veta                                      # n x k 
    eta = Meta + matrix(rnorm(n*k), nrow = n, ncol = k) %*% t(S)    # update eta in a block
    
    # -- update Lambda (rue & held) -- #
    eta2 = t(eta) %*% eta    # prepare eta crossproduct before the loop
    zlams = rnorm(k*p)       # generate normal draws all at once 
    
    for(j in 1:p) {
      Llamt = chol(diag(Plam[j,]) + ps[j]*eta2)
      Lambda[j,] = t(solve(Llamt,
                           zlams[1:k + (j-1)*k]) + 
                       solve(Llamt,
                             solve(t(Llamt),
                                   ps[j] * t(eta) %*% Y[,j])))
    }
    
    #------Update psi_{jh}'s------#
    psijh = matrix(rgamma(p*k,
                          df/2 + 0.5,
                          df/2 + t(t(Lambda)^2 * (tauh))),
                   nrow = p, ncol = k)
    
    #------Update theta & tauh------#
    mat = psijh * Lambda^2
    ad = ad1 + 0.5*p*k
    bd = bd1 + 0.5 * theta[1] * sum(tauh*colSums(mat))
    theta[1] = 1 / rgamma(1,ad,bd)           
    tauh = 1 / cumprod(theta)
    
    
    for(h in 2:k) {
      ad = ad2 + 0.5*p*(k-h+1)
      bd = bd2 + 0.5 * theta[h] * sum(tauh[h:k]*colSums(mat[,h:k, drop = F]))
      theta[h] = 1 / rgamma(1,ad,bd)
      tauh = 1 / cumprod(theta)
    }
    
    # -- Update Sigma -- #
    Ytil = Y - eta %*% t(Lambda)
    ps= rgamma(p, as + 0.5*n, bs+0.5*colSums(Ytil^2))
    Sigma=diag(1/ps)
    
    #---update precision parameters----#
    Plam = t(t(psijh) * tauh)
    
    # -- save sampled values (after thinning) -- #
    if((i %% thin == 0) & (i > burn)) {
      Omega = (Lambda %*%  t(Lambda) + Sigma) * scaleMat
      if(any(output %in% "covMean")) COVMEAN = COVMEAN + Omega / sp
      if(any(output %in% "covSamples")) OMEGA[,,ind] = Omega
      if(any(output %in% "factSamples")) LAMBDA[[ind]] = Lambda
      if(any(output %in% "sigSamples")) SIGMA[,,ind] = Sigma
      ind = ind + 1
    }
    
    if((i %% 1000) == 0) {
      cat(i,"\n")
    }
  }
  out = lapply(output, function(x) {
    if(x == "covMean") return(COVMEAN)
    if(x == "covSamples") {
      saveRDS(OMEGA, file = covfilename)
      return(paste("see", covfilename))
    }
    if(x == "factSamples") {
      saveRDS(LAMBDA, file = factfilename)
      return(paste("see", factfilename))
    }
    if(x == "sigSamples") {
      saveRDS(SIGMA, file = sigfilename)
      return(paste("see", sigfilename))
    }
    if(x == "numFactors") return(K)
  })
  names(out) = output
  return(out)
}

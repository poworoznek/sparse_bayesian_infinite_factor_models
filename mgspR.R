# Gibbs sampler for factorized covariance estimation
# using mgps prior on factor loadings from Battacharya Dunson (2011)
# prior hyperparameters informed by Durante (2017)
# version for adaptive number of factors weighted by iteration
# see https://github.com/poworoznek/sparse_bayesian_infinite_factor_models for license

# ARGUMENTS: Y: Data matrix (n x p); 
#            prop: proportion of elements in each column less than epsilon in magnitude cutoff;
#            epsilon: tolerance;
#            nrun: number of iterations;
#            burn: burn-in period;
#            thin: thinning interval;
#            kinit: initial value for the number of factors;
#            output: output type, a vector including some of:
#            c("covMean", "covSamples", "factSamples", "sigSamples", "numFactors");
#            covfilename: optional filename for covariance matrix samples;
#            factfilename: optional filename for factor matrix samples;
#            sigfilename: optional filename for sigma matrix samples;
#            verbose: logical. print sample checkpoints?
#            adapt: logical. Adapt k across samples?

mgspR = function(Y, prop = 1, epsilon = 1e-3, nrun, burn, thin = 1, 
                    kinit = NULL, output = c("covMean", "covSamples", 
                                             "factSamples", "sigSamples",
                                             "numFactors"), 
                    covfilename = "Omega.rds", factfilename = "Lambda.rds", 
                    sigfilename = "Sigma.rds", verbose = TRUE, adapt = TRUE){
  
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
  
  if(is.null(kinit)) kinit = floor(log(p)*3)
  
  sp = floor((nrun - burn)/thin)        # number of posterior samples
  
  VY= apply(Y, 2, var)                  # explicitly preserve scale
  scaleMat = sqrt((VY) %*% t(VY))
  Y = scale(Y)
  num = 0
  k=kinit                               # no. of factors to start with
  
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
  if(any(output %in% "numFactors")) K = rep(NA, sp)
  ind = 1
  
  #------start gibbs sampling-----#
  
  if(verbose) cat("Start\n")
  
  for(i in 1:nrun) {
    
    # -- Update eta -- #
    Lmsg = Lambda * ps
    Veta1 = diag(k) + crossprod(Lmsg, Lambda)
    Tmat = chol(Veta1)
    S = backsolve(Tmat, diag(k))
    Veta = tcrossprod(S)                                              # Veta = inv(Veta1)
    Meta = Y %*% Lmsg %*% Veta                                      # n x k 
    eta = Meta + tcrossprod(matrix(rnorm(n*k), 
                                   nrow = n, 
                                   ncol = k), S)                    # update eta in a block
    
    # -- update Lambda (rue & held) -- #
    eta2 = crossprod(eta)    # prepare eta crossproduct before the loop
    zlams = rnorm(k*p)       # generate normal draws all at once 
    
    for(j in 1:p) {
      Llamt = chol(diag(Plam[j,]) + ps[j]*eta2)
      Lambda[j,] = t(backsolve(Llamt, zlams[1:k + (j-1)*k]) + 
                       backsolve(Llamt,
                             forwardsolve(t(Llamt),
                                   ps[j] * crossprod(eta, Y[,j]))))
    }
    
    #------Update psi_{jh}'s------#
    psijh = matrix(rgamma(p*k,
                          df/2 + 0.5,
                          df/2 + t(t(Lambda)^2 * (tauh))/2),
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
    Ytil = Y - tcrossprod(eta, Lambda)
    ps= rgamma(p, as + 0.5*n, bs+0.5*colSums(Ytil^2))
    Sigma = diag(1/ps)
    
    #---update precision parameters----#
    Plam = t(t(psijh) * tauh)
    
    # ----- make adaptations ----#
    prob = 1/exp(b0 + b1*i)                    # probability of adapting
    uu = runif(1)
    lind = colSums(abs(Lambda) < epsilon)/p    # proportion of elements in each column less than eps in magnitude
    vec = lind >= prop
    num = sum(vec)                             # number of redundant columns
    
    if((uu < prob) & adapt) {
      if((i > 20) & (num == 0) & all(lind < 0.995)) {
        k = k + 1
        Lambda = cbind(Lambda, rep(0,p))
        eta = cbind(eta,rnorm(n))
        psijh = cbind(psijh, rgamma(p,df/2,df/2))
        theta[k] = 1 / rgamma(1, ad2,bd2)
        tauh = 1 / cumprod(theta)
        Plam = t(t(psijh) * tauh)
      } else {
        if (num > 0) {
          k = max(k - num,1)
          Lambda = Lambda[,!vec, drop = F]
          psijh = psijh[,!vec, drop = F]
          eta = eta[,!vec, drop = F]
          theta = theta[!vec]
          tauh = 1 / cumprod(theta)
          Plam = t(t(psijh) * tauh)
        }
      }
    }
    
    # -- save sampled values (after thinning) -- #
    if((i %% thin == 0) & (i > burn)) {
      if(any(output %in% c("covMean", "covSamples"))) Omega = (tcrossprod(Lambda) + Sigma) * scaleMat
      if(any(output %in% "covMean")) COVMEAN = COVMEAN + Omega / sp
      if(any(output %in% "covSamples")) OMEGA[,,ind] = Omega
      if(any(output %in% "factSamples")) LAMBDA[[ind]] = Lambda
      if(any(output %in% "sigSamples")) SIGMA[,,ind] = Sigma
      if(any(output %in% "numFactors")) K[ind] = k
      ind = ind + 1
    }
    
    if(((i %% 1000) == 0) & verbose) {
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

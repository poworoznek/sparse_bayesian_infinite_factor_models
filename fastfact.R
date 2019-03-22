# Gibbs sampler for factorized covariance estimation
# using mgps prior on factor loadings from Battacharya Dunson (2011)
# prior edits informed by Durante (2017)
# version for adaptive number of factors weighted by iteration

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

fastfact = function(Y, prop = 1, epsilon = 1e-3, nrun, burn, thin = 1, 
                 kinit = NULL, output = "covMean", 
                 covfilename = "Omega.rds", factfilename = "Lambda.rds", 
                 sigfilename = "Sigma.rds",
                 dump = FALSE, buffer = 10000){
  
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
  meta = matrix(0,nrow = n, ncol = k)
  veta = diag(k)                                 # latent factor distribution = standard normal
  
  psijh = matrix(rgamma(p*k, df/2, df/2), nrow = p, ncol = k)     # local shrinkage coefficients
  theta = c(rgamma(1,ad1,bd1), rgamma(k-1,ad2,bd2))       # gobal shrinkage coefficients multilpliers
  tauh = cumprod(theta)                                       # global shrinkage coefficients
  Plam = t(t(psijh) * (tauh))                                     # precision of loadings rows
  start = 0
  
  if(!dump){
    coutput = MGSPsamp(p, n, k, as, bs, df, ad1, bd1, ad2, bd2, adf, bdf, 
                       b0, b1, sp, nrun, burn, thin, prop, epsilon, ps, Sigma, Lambda, 
                       meta, veta, psijh, theta, tauh, Plam, Y, scaleMat, output, start) 
    COVMEAN = coutput$covMean
    OMEGA = coutput$covSamps
    LAMBDA = coutput$factSamps
    SIGMA = coutput$sigSamps
    K = coutput$numFact
    
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
  }else{
    
    runs = nrun %/% buffer + as.logical(nrun %% buffer)
    coutlist = list()
    oldnrun = nrun
    nrun = min(buffer, nrun)
    sp = max(floor((nrun - burn)/thin), 1)
    
    for(run in 1:runs){
      coutput = MGSPsamp(p, n, k, as, bs, df, ad1, bd1, ad2, bd2, adf, bdf, 
                         b0, b1, sp, nrun, burn, thin, prop, epsilon, ps, Sigma, Lambda, 
                         meta, veta, psijh, theta, tauh, Plam, Y, scaleMat, output, start)
      start = start + buffer
      burned = burn - buffer
      burn = max(1, burned)
      nrun = min(buffer, oldnrun - buffer*run)
      sp = max(floor((nrun - burn)/thin), 0)
      
      ps = coutput[["lastState"]][[1]]
      Sigma = coutput[["lastState"]][[2]]
      Lambda = coutput[["lastState"]][[3]]
      meta = coutput[["lastState"]][[4]]
      veta = coutput[["lastState"]][[5]]
      psijh = coutput[["lastState"]][[6]]
      theta = coutput[["lastState"]][[7]]
      tauh = coutput[["lastState"]][[8]]
      Plam = coutput[["lastState"]][[9]]
      k = coutput[["lastState"]][[10]]
      out = lapply(output, function(x) {
        if(x == "covMean") return(coutput$covMean)
        if(x == "covSamples") {
          dump(coutput$covSamps, file = covfilename, append = TRUE)
          return(paste("see", covfilename))
        }
        if(x == "factSamples") {
          dump(coutput$factSamps, file = factfilename, append = TRUE)
          return(paste("see", factfilename))
        }
        if(x == "sigSamples") {
          dump(coutput$sigSamps, file = sigfilename, append = TRUE)
          return(paste("see", sigfilename))
        }
        if(x == "numFactors") return(coutput$numFact)
      })
      names(out) = output
      coutlist[[run]] = out
    }
    return(coutlist)
  }
  

}

# Gibbs sampler for factorized regression with all 2 way interactions
# using MGSP prior on factor loadings from Battacharya Dunson (2014)
# same eta, lambda, and sigma updates as interactionDL

# ARGUMENTS: Y: Data matrix (n x p); 
#            nrun: number of iterations;
#            burn: burn-in period;
#            thin: thinning interval;
#            prop: proportion of elements in each column less than epsilon in magnitude cutoff;
#            epsilon: tolerance;
#            kinit: initial value for the number of factors;
#            output: output type, a vector including some of:
#            c("covMean", "covSamples", "factSamples", "sigSamples", "coefSamples", "numFactors", "errSamples");
#            covfilename: optional filename for covariance matrix samples;
#            factfilename: optional filename for factor matrix samples;
#            sigfilename: optional filename for sigma matrix samples;
#            adapt: logical or "burn". Adapt proposal variance in metropolis hastings step? if "burn", 
#                   will adapt during burn in and not after

interactionMGSP = function(y, X, nrun, burn, thin = 1, delta_rw = 0.0526749, epsilon_rw = 0.5,
                         a = 1/2, kinit = NULL, output = c("covMean", "covSamples", "factSamples",
                                                           "sigSamples", "coefSamples", "numFactors",
                                                           "errSamples"),  covfilename = "Omega.rds",
                         factfilename = "Lambda.rds",  sigfilename = "Sigma.rds", verbose = TRUE,
                         dump = FALSE, buffer = 10000, adapt = "burn"){
  
  cm = any(output %in% "covMean")
  cs = any(output %in% "covSamples")
  fs = any(output %in% "factSamples")
  ss = any(output %in% "sigSamples")
  cf = any(output %in% "coefSamples")
  nf = any(output %in% "numFactors")
  es = any(output %in% "errSamples")
  
  p = ncol(X)
  n = nrow(y)
  
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
  
  VX= apply(X, 2, var)                  # explicitly preserve scale
  scaleMat = sqrt((VX) %*% t(VX))
  X = scale(X)
  k=kinit                               # no. of factors to start with
  
  ssy = 1                     # initial values
  phi = numeric(k)
  ps = rgamma(p,as,bs)
  lambda = matrix(0,p,k)
  eta = matrix(rnorm(n*k),n,k)
  Psi = matrix(0,k,k)
  
  psiMG = matrix(rgamma(p*k, df/2, df/2), nrow = p, ncol = k)     # local shrinkage coefficients
  delta = c(rgamma(1,ad1,bd1), rgamma(k-1,ad2,bd2))       # gobal shrinkage coefficients multilpliers
  tau = cumprod(delta)                                       # global shrinkage coefficients
  Plam = t(t(psiMG) * (tau))                                     # precision of loadings rows

  
  if(cm) COVMEAN = matrix(0, nrow = p, ncol = p)
  if(cs) OMEGA = array(dim = c(p, p, sp))
  if(fs) LAMBDA = list()
  if(ss) SIGMA = array(dim = c(p, sp))
  if(cf) {
    PHI = array(dim = c(k, sp))
    PSI = array(dim = c(k, k, sp))
    INTERCEPT = numeric(sp)
    MAIN = array(dim = c(p, sp))
    INTERACTION = array(dim = c(p, p, sp))
  }
  if(nf) K = numeric(sp)
  if(es) SSY = numeric(sp)
  ind = 1
  acp = numeric(n)
  
  if(verbose) {
    pb = txtProgressBar(style = 3)
    at = ceiling(nrun/100)
  }
  
  for(i in 1:nrun){
    eta = eta_int(lambda, eta, ps, phi, Psi, k, n, y, X, ssy, delta_rw, acp)
    Psi = psi_int(eta, y, phi, ssy, k, n)
    phi = phi_int(eta, y, ssy, Psi, k)
    ssy = ssy_int(eta, phi, Psi, y, n)
    psiMG = psi_mg(lambda, tau, ps, k, p, df)
    delta = del_mg(lambda, psiMG, tau, delta, k, p, ad1, bd1, ad2, bd2)
    tau = cumprod(delta)
    Plam = plm_mg(psiMG, tau)
    lambda = lam_lin(eta, Plam, ps, k, p, X)
    ps = c(sig_lin(lambda, eta, k, p, n, X, as, bs))
    
    if((i %% thin == 0) & (i > burn)) {
      if(cm | cs) Omega = (tcrossprod(lambda) + diag(1/ps)) * scaleMat
      if(cm) COVMEAN = COVMEAN + Omega / sp
      if(cs) OMEGA[,,ind] = Omega
      if(fs) LAMBDA[[ind]] = lambda
      if(ss) SIGMA[,ind] = 1/ps
      if(cf) {
        V_n = solve(t(lambda/ps)%*% lambda+diag(k))
        a_n = V_n %*% t(lambda/ps)
        dsVX = diag(sqrt(VX))
        dsVX_inv = diag(1/sqrt(VX))
        PHI[,ind] = phi
        PSI[,,ind] = Psi
        INTERCEPT[ind] = sum(diag(Psi%*%V_n))
        MAIN[,ind] = as.vector(t(phi)%*%a_n)
        INTERACTION[,,ind] = dsVX_inv%*%t(a_n)%*%Psi%*%a_n%*%dsVX_inv
      }
      if(nf) K[ind] = k
      if(es) SSY[ind] = ssy
      ind = ind + 1
    }
    
    if(verbose & (i %% at == 0)) setTxtProgressBar(pb, i / nrun)
    
    if (i%%100==0) {if(adapt == "burn") {if(i <= burn) {
      acp_mean = mean(acp)/100
      if(acp_mean > 0.3){
        delta_rw = delta_rw*2
      }else if(acp_mean < 0.2){
        delta_rw = delta_rw*2/3}
      acp = numeric(n)
    }} else if(adapt) {
      acp_mean = mean(acp)/100
      if(acp_mean > 0.3){
        delta_rw = delta_rw*2
      }else if(acp_mean < 0.2){
        delta_rw = delta_rw*2/3
      }
      acp = numeric(n)
    }}
    
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
  if(cf) out = c(out, list(phiSamps = PHI, PsiSamps = PSI, 
                           interceptSamps = INTERCEPT,
                           mainEffectSamps = MAIN,
                           interactionSamps = INTERACTION))
  if(nf) out = c(out, list(numFacts = K))
  if(es) out = c(out, list(ssySamps = SSY))
  
  return(out)
}

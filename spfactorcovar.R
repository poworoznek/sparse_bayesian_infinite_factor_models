# MATLAB code by Anirban Bhattacharya
# Ported to R by Evan Poworoznek

# Gibbs sampler for covariance estimation
# using mgps prior on factor loadings

rm(list = ls())

library(R.matlab)
list2env(readMat("~/documents/factorR/datagenp_30ktr_5rep_1.mat"), envir = .GlobalEnv)
p = as.numeric(p)
n = as.numeric(n)

# --- define global constants --- #
nrun=20000
burn=5000
thin=1

Rcpp::sourceCpp("MGSPc.cpp")
source("fastfact.R")
source("safefact.R")

t1 = system.time({output = safefact(Y = dat, nrun = nrun, burn = burn, verbose = FALSE)})
msesafe = mean((output$covMean - Ot)^2)

t2 = system.time({output2 = fastfact(Y = dat, nrun = nrun, burn = burn, verbose = FALSE)})
msecpp = mean((output2$covMean - Ot)^2)

t1;t2

sp =(nrun - burn)/thin # number of posterior samples

kinit = rep(floor(log(p)*3),rep)                 # number of factors to start with
b0 = 1
b1 = 0.0005
epsilon = 1e-3                                        # threshold limit
prop = 1.00                                           # proportion of redundant elements within columns


#---- define output files across replicates-----#

mserep = matrix(0,nrow = rep,ncol = 3)                          # mse,absolute bias(avg and maximum) in estimating cov matrix
mse1rep = matrix(0,nrow = rep,ncol = 3)                         # same as above in original scale in estimating cov matrix
nofrep = matrix(0,nrow = rep,ncol = sp)                         # evolution of factors across replicates
adrep = matrix(0,nrow = rep,ncol = 1)                           # number of adaptations across replicates
OMEGA = array(dim = c(p, p, sp))

#library(profvis)
t3 = system.time({
  for(g in 1:rep) {
    
    cat("start replicate",g, "\n")
    cat("--------------------\n")
    
    # ------read data--------#
    Y = dat[(g-1)*n + 1:(g*n),]               # n x p data matrix ; scale the training data
    VY= apply(Y, 2, var)
    M = colMeans(Y)
    Y = scale(Y)
    Ot1 = Ot / sqrt((VY) %*% t(VY))                # true dispersion matrix of the transformed data
    num = 0
    k=kinit[g]                                 # no. of factors to start with
    
    # ------end read data--------#
    
    # --- Define hyperparameter values --- #
    
    as = 1                          # gamma hyperparameters for residual precision
    bs = 0.3                        
    df = 3                                    # gamma hyperparameters for t_{ij}
    ad1 = 2.1
    bd1 = 1                         # gamma hyperparameters for delta_1
    ad2 = 3.1
    bd2 = 1                         # gamma hyperparameters delta_h, h >= 2
    adf = 1
    bdf = 1                          # gamma hyperparameters for ad1 and ad2 or df
    
    # --- Initial values --- #
    ps = rgamma(p, as, bs)
    Sigma=diag(1/ps)                  # Sigma = diagonal residual covariance
    Lambda = matrix(1, nrow = p, ncol = k)
    ta = matrix(rnorm(n*k), nrow = n, ncol = k)         # factor loadings & latent factors
    meta = matrix(0,nrow = n, ncol = k)
    veta = diag(k)                           # latent factor distribution = standrad normal
    
    psijh = matrix(rgamma(p*k, df/2, df/2), nrow = p, ncol = k)                            # local shrinkage coefficients
    delta = c(rgamma(1,ad1,bd1), rgamma(k-1,ad2,bd2)) ### bd1, bd2 as scale parameters ???  # gobal shrinkage coefficients multilpliers
    tauh = cumprod(delta)                                      # global shrinkage coefficients
    Plam = t(t(psijh) * (tauh))                         # precision of loadings rows
    
    # --- Define output files specific to replicate --- #
    nofout = rep(0, nrun+1)                  # number of factors across iteartions
    nofout[1] = k
    nof1out = rep(0,sp)
    mseout = matrix(0, nrow = sp, ncol = 3)                      # within a replicate, stores mse across mcmc iterations
    mse1out = matrix(0, nrow = sp, ncol = 3)                     # mse in original scale
    Omegaout = rep(0, p^2)
    Omega1out = rep(0,p^2)
    
    
    
    
    #------start gibbs sampling-----#
    
    for(i in 1:nrun) {
      
      # -- Update eta -- #
      Lmsg = Lambda * ps
      Veta1 = diag(k) + t(Lmsg) %*% Lambda
      eigs = eigen(Veta1)
      if(all(eigs$values > 1e-6)) {
        Tmat = sqrt(eigs$values) * t(eigs$vectors)
      } else {
        Tmat = chol(Veta1)
      }
      R = qr.R(qr(Tmat))
      S = solve(R)
      Veta = S %*% t(S)                   # Veta = inv(Veta1)
      Meta = Y %*% Lmsg %*% Veta                        # n x k 
      eta = Meta + matrix(rnorm(n*k), nrow = n, ncol = k) %*% t(S) # update eta in a block
      
      # -- update Lambda (rue & held) -- #
      eta2 = t(eta) %*% eta
      
      ### translate to mclapply() or openCL after profvis ???
      for(j in 1:p) {
        Qlam = diag(Plam[j,]) + ps[j]*eta2
        blam = ps[j]*(t(eta) %*% Y[,j])
        Llam = t(chol(Qlam))
        zlam = rnorm(k)
        vlam = solve(Llam,blam)
        mlam = solve(t(Llam),vlam)
        ylam = solve(t(Llam),zlam)
        Lambda[j,] = t(ylam + mlam)
      }
      
      #------Update psi_{jh}'s------#
      psijh = matrix(nrow = p, ncol = k,
                     data = rgamma(p*k,df/2 + 0.5,df/2 + t(t(Lambda)^2 * (tauh))))
      
      #------Update delta & tauh------#
      mat = psijh * Lambda^2
      ad = ad1 + 0.5*p*k
      bd = bd1 + 0.5 / delta[1] * sum(tauh*colSums(mat))
      delta[1] = rgamma(1,ad,bd)
      tauh = cumprod(delta)
      
      
      for(h in 2:k) {
        ad = ad2 + 0.5*p*(k-h+1)
        bd = bd2 + 0.5 / delta[h] * sum(tauh[h:k]*colSums(mat[,h:k, drop = F]))
        delta[h] = rgamma(1,ad,bd)
        tauh = cumprod(delta)
      }
      
      # -- Update Sigma -- #
      Ytil = Y - eta %*% t(Lambda) 
      ps= rgamma(p, as + 0.5*n, bs+0.5*colSums(Ytil^2))
      Sigma=diag(1/ps)
      
      #---update precision parameters----#
      Plam = t(t(psijh) * tauh)
      
      # ----- make adaptations ----#
      prob = 1/exp(b0 + b1*i)                # probability of adapting
      uu = runif(1)
      lind = colSums(abs(Lambda) < epsilon)/p    # proportion of elements in each column less than eps in magnitude
      vec = lind >= prop
      num = sum(vec)       # number of redundant columns
      
      if(uu < prob) {
        if((i > 20) & (num == 0) & all(lind < 0.995)) {
          k = k + 1
          Lambda = cbind(Lambda, rep(0,p))
          eta = cbind(eta,rnorm(n))
          psijh = cbind(psijh, rgamma(p,df/2,df/2))
          delta[k] = rgamma(1, ad2,bd2)
          tauh = cumprod(delta)
          Plam = t(t(psijh) * tauh)
        } else {
          if (num > 0) {
            k = max(k - num,1)
            Lambda = Lambda[,!vec, drop = F]
            psijh = psijh[,!vec, drop = F]
            eta = eta[,!vec, drop = F]
            delta = delta[!vec]
            tauh = cumprod(delta)
            Plam = t(t(psijh) * tauh)
          }
        }
      }
      nofout[i+1] = k
      
      # -- save sampled values (after thinning) -- #
      if((i %% thin == 0) & (i > burn)) {
        Omega = Lambda %*% t(Lambda) + Sigma
        Omega1 = Omega * sqrt((VY) %*% t(VY))
        Omegaout = Omegaout + as.vector(Omega)/sp
        Omega1out = Omega1out + as.vector(Omega1)/sp
        nof1out[(i-burn)/thin] = (nofout[(i-burn)/thin] - num)*(num > 0)
        OMEGA[,,(i-burn)/thin] = Omega1
      }
      
      if((i %% 1000) == 0) {
        cat(i,"\n")
      }
    }
    
    
    #---- summary measures specifi to replicate ---#
    #1. covariance matrix estimation
    errcov = Omegaout - as.vector(Ot1)
    err1cov = Omega1out - as.vector(Ot)
    mserep[g,] = c(mean(errcov^2), mean(abs(errcov)), max(abs(errcov)))
    mse1rep[g,] = c(mean(err1cov^2), mean(abs(err1cov)), max(abs(err1cov)))
    
    #w. evolution of factors
    nofrep[g,] = t(nof1out)
    
    cat("end replicate",g, "\n")
    cat("--------------------\n")

  }
  
})

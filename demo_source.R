# Source All Necessary Functions for the DEMO

generate_DL_reg = function(S = 10000, k = 15, p = 20,
                           a_DL = 0.5, var_omega = 5,
                           as = 1, bs = 0.3,
                           trace = FALSE, every = 500){
  
  # Storage
  # third order interactions in the model
  beta1_X_DL =  matrix(0,p,S)
  beta2_X_DL = array(0,c(p,p,S))
  beta3_X_DL = array(0,c(1,2,3,S))
  
  
  for (i in 1:S){
    
    Lambda_DL = matrix(0,p,k)
    ps = rgamma(p,as,bs)
    for(j in 1:p){
      phi_DL = rdirichlet(1, rep(a_DL,p))
      for(h in 1:k){
        
        #Dirichlet-Laplace
        psi = rexp(1,0.5)
        tau_DL = rgamma(1,k*a_DL,0.5)
        Lambda_DL[j,h] = rnorm(1,0,tau_DL*phi_DL[j]*sqrt(psi))
      }
    }
    
    omega1 = rnorm(k,0,var_omega)
    omega2 = rnorm(k,0,var_omega)
    omega3 = rnorm(k,0,var_omega)
    A_DL = solve(t(Lambda_DL)%*%diag(ps)%*%Lambda_DL + 
                   diag(k))%*%t(Lambda_DL)%*%diag(ps)
    
    
    
    # up to 3rd order interactions
    beta1_X_DL[,i] = t(A_DL)%*%omega1 + 3*diag(ps)%*%t(A_DL)%*%omega3
    
    beta2_X_DL[,,i] = 2*(2*t(A_DL)%*%diag(omega2)%*%A_DL + 
                           6*diag(1/ps)%*%t(A_DL)%*%diag(omega2)%*%A_DL)
    
    beta3_X_DL[1,2,3,i] = 6*sum(omega3*A_DL[,1]*A_DL[,2]*A_DL[,3])
    
    if(trace){
      if (i%%every==0){
        print(i)
      }
    }
    
  }
  
  
  xlim = c(-3,3)
  par(mfrow = c(1,3))
  hist(beta1_X_DL[1,],freq = F,breaks = seq(-1000,1000,by =0.05),xlim = xlim,
       xlab = "main effect",main = "")
  abline(v = quantile(beta1_X_DL[1,],probs = c(0.05,0.95)),col="red",lty="dotted")
  abline(v = quantile(beta1_X_DL[1,],probs = c(0.25,0.75)),col="green",lty="dotted")
  hist(beta2_X_DL[1,2,],freq = F,breaks = seq(-1000,1000,by =0.05),xlim = xlim,
       xlab = "2nd order interaction",main =  paste("k =",k,",","p =",p))
  abline(v = quantile(beta2_X_DL[1,2,],probs = c(0.05,0.95)),col="red",lty="dotted")
  abline(v = quantile(beta2_X_DL[1,2,],probs = c(0.25,0.75)),col="green",lty="dotted")
  hist(beta2_X_DL[1,2,],freq = F,breaks = seq(-1000,1000,by =0.05),xlim = xlim,
       xlab = "3rd order interaction",main = "")
  abline(v = quantile(beta3_X_DL[1,2,3,],probs = c(0.05,0.95)),col="red",lty="dotted")
  abline(v = quantile(beta3_X_DL[1,2,3,],probs = c(0.25,0.75)),col="green",lty="dotted")
  
  
}

# function to perform factorized regression with interactions
# See Federico's paper (2019?)
# Dirichlet Laplace prior with random walk met-hastings

# ARGUMENTS: y: vector of responses;
#            X: Matrix of predictors
#            delta_rw: MH tuning parameter
#            epsilon_rw:  MH tuning parameter
#            a: Dirichlet Laplace prior hyperparameter
#            k: number of factors

gibbs_DL = function(y, X ,nrun, burn, thin = 1, 
                    delta_rw = 0.0526749, epsilon_rw = 0.5,
                    a = 1/2, k = NULL){
  n = nrow(X)
  p = ncol(X)                    # collect data attributes
  
  if(is.null(k)) k = floor(log(p)*3)
  if(length(y) != n) stop("Mismatching input lengths")
  
  VX = apply(X, 2, var)          # prepare data matrix
  #VX = rep(1,p)
  X = as.matrix(scale(X))
  
  sp = floor((nrun - burn)/thin)
  
  b0 = 1            # factorization hyperparameters
  b1 = 0.0005
  epsilon = 1e-3
  prop = 1.00
  as = 1
  bs = 0.3
  df = 150
  ad1 = 3
  bd1 = 1
  ad2 = 4
  bd2 = 1
  adf = 1
  bdf = 1
  
  Omega_bayes = array(0,c(sp,p,p))  # sample storage memory allocation
  tau_st = array(0,c(sp,p))
  sigmasq_st = numeric(sp)
  alpha_bayes = numeric(sp)
  beta_bayes = matrix(0,sp,p)
  Lambda_st = array(0,c(sp,p,k))
  acp = numeric(n)
  
  sigmasq_y = 1                     # initial values
  phi = numeric(k)
  ps = rgamma(p,as,bs)
  Sigma = diag(1/ps)
  Lambda = matrix(0,p,k)
  Lambda.T = t(Lambda)
  eta = matrix(rnorm(n*k),n,k)
  meta = matrix(0,n,k)
  veta = diag(rep(1,k))
  Psi = matrix(0,k,k)
  
  tau = rgamma(p,p*a,1/2)
  psijh = matrix(rexp(p*k,1/2),p,k)
  phijh = matrix(0,p,k)
  for(j in 1:p){
    phijh[j,] = rdirichlet(1,rep(a,k))
  }
  Plam = psijh*(phijh^2)*matrix(rep(tau^2,k),p,k,byrow=F)
  
  t = t0 = Sys.time()               # begin sample timing
  count = 1 
  
  for(i in 1:nrun){
    
    # --- Update eta --- #
    aMH = phi%*%t(phi)/sigmasq_y + Lambda.T%*%diag(ps)%*%Lambda + diag(k)
    
    for (h in 1:n){                # Metropolis hastings step 
      
      eta_star = bayesSurv::rMVNorm(1,eta[h,],diag(k)*delta_rw)
      eta_star.T = t(eta_star)     # avoid repeated transpose calls
      eta.T = t(eta[h,])
      
      logr = eta_star.T%*%(aMH - 2*Psi*y[h]/sigmasq_y)%*%eta_star -
        2*eta_star.T%*%(Lambda.T%*%diag(ps)%*%X[h,] + phi*y[h]/sigmasq_y) +
        2*eta_star.T%*%phi*(eta_star%*%Psi%*%eta_star)/sigmasq_y + 
        (1/sigmasq_y)*(eta_star%*%Psi%*%eta_star)^2 - 
        (eta.T%*%(aMH - 2*Psi*y[h]/sigmasq_y)%*%eta[h,] -
           2*eta.T%*%(Lambda.T%*%diag(ps)%*%X[h,] + phi*y[h]/sigmasq_y) +
           2*eta.T%*%phi*(eta[h,]%*%Psi%*%eta[h,])/sigmasq_y +
           (1/sigmasq_y)*(eta[h,]%*%Psi%*%eta[h,])^2)
      logr = logr*(-0.5)
      
      logu = log(runif(1))
      
      if (logr > logu){
        eta[h,] = eta_star
        acp[h] = acp[h] + 1
      }
    }
    
    # --- Update Psi --- #
    MM = model.matrix(y~.^2 - 1,as.data.frame(eta))   # perform factorized regression
    X_reg = cbind(eta^2,MM[,(k+1):ncol(MM)])
    X_reg.T = t(X_reg)
    Lambda_n = X_reg.T%*%X_reg/sigmasq_y + diag(rep(1,ncol(X_reg)))
    Vcsi = solve(Lambda_n)
    Mcsi = Vcsi%*%X_reg.T%*%(y-eta%*%phi)/sigmasq_y
    csi = bayesSurv::rMVNorm(n=1,mean=Mcsi,Sigma=Vcsi)
    lambda_diag = csi[1:k]
    lambda = csi[(k+1):length(csi)]
    Psi[lower.tri(Psi)] = lambda/2
    Psi[upper.tri(Psi)] = 0
    Psi = Psi + t(Psi)
    diag(Psi) = lambda_diag
    
    # --- Update phi --- #
    eta.T = t(eta)
    Lambda_n = eta.T%*%eta/sigmasq_y + diag(rep(1,ncol(eta)))/2
    Vcsi = solve(Lambda_n)
    Mcsi = Vcsi%*%eta.T%*%(y-diag(eta%*%Psi%*%eta.T))/sigmasq_y     # using updated psi
    phi = bayesSurv::rMVNorm(n = 1, mean = Mcsi, Sigma = Vcsi)
    
    # --- Update sigmasq_y --- #
    an = 0.5 + n/2
    bn = 0.5 + 0.5*t(y-eta%*%phi-diag(eta%*%Psi%*%eta.T))%*%
      (y-eta%*%phi-diag(eta%*%Psi%*%eta.T))
    sigmasq_y = 1/rgamma(1,an,bn)
    
    # --- Update Lambda --- #
    Plam = 1/(psijh*(phijh^2)*matrix(rep(tau^2,k),p,k,byrow=F))
    eta2 = eta.T%*%eta
    zlams = rnorm(k*p, sd = 1)       # generate normal draws all at once 
    
    for(j in 1:p) {
      Llamt = chol(diag(Plam[j,]) + ps[j]*eta2)
      Lambda[j,] = t(solve(Llamt,
                           zlams[1:k + (j-1)*k]) + 
                       solve(Llamt,
                             solve(t(Llamt),
                                   ps[j] * eta.T %*% X[,j])))
    }
    Lambda.T = t(Lambda)
    
    # --- Update psijh --- #
    mujh = phijh*matrix(rep(tau,k),p,k,byrow = F)/abs(Lambda)
    psijh = matrix(rinvgauss(k*p, mujh),p,k,byrow = T)
    
    # --- Update tau --- #
    for(j in 1:p){
      tau[j] = GIGrvg::rgig(n=1,lambda = 1-k, psi = 1, 
                            chi = 2*sum(abs(Lambda[j,])/phijh[j,]))
    }
    
    # --- Update phijh --- #
    Tjh = numeric(k)
    for(j in 1:p){
      for(h in 1:k){
        Tjh[h] = GIGrvg::rgig(n=1,lambda = a-1,psi = 1,
                              chi = 2*abs(Lambda[j,h]))
      }
      phijh[j,] = Tjh/sum(Tjh)
    }
    
    # --- Update Sigma --- #
    Xtil = X - eta%*%Lambda.T
    ps = rgamma(p,as+0.5*n,1)
    ps = (1/(bs + 0.5*apply(Xtil^2,2,sum)))*ps
    Sigma = diag(1/ps)
    
    #store posterior samples
    if ((i %% thin == 0) & (i > burn)){
      
      # parameters of the factor model
      sigmasq_st[count] = sigmasq_y
      tau_st[count,] = tau
      
      #Bayesian estimators
      V_n = solve(Lambda.T%*%solve(Sigma)%*%Lambda+diag(rep(1,ncol(Lambda))))
      a_n = V_n%*%Lambda.T%*%solve(Sigma)
      dsVX = diag(sqrt(VX))
      dsVX_inv = diag(1/sqrt(VX))
      Omega_bayes[count,,] = dsVX_inv%*%t(a_n)%*%Psi%*%a_n%*%dsVX_inv
      beta_bayes[count,] = as.vector(t(phi)%*%a_n)
      alpha_bayes[count] = tr(Psi%*%V_n)
      Lambda_st[count,,] = dsVX %*% Lambda
      count = count + 1
    }
    
    if (i%%100==0){
      #    print(i)
      acp_mean = mean(acp)/100
      # print(acp_mean)
      if(acp_mean > 0.3){
        delta_rw = delta_rw*2
      }else if(acp_mean < 0.2){
        delta_rw = delta_rw*2/3
      }
      acp = numeric(n)
      #   print(delta_rw)
      #   #print(paste("time for last 100 iterations:",round(as.numeric(Sys.time()-t),0),
      #   #            "seconds",sep=" "))
      #   #t = Sys.time()
      #   #print(paste(i,"out of",nrun,sep=" "))
      #   #t_end = round(((as.numeric(difftime(Sys.time(),t0,units="secs")))/i)*(nrun-i),0)
      #   #print(paste("estimated time to end:",t_end,"seconds",sep=" "))
    }
  }
  
  beta_bayes = beta_bayes%*%dsVX_inv
  return(list(alpha_bayes = alpha_bayes,
              beta_bayes = beta_bayes,
              Omega_bayes = Omega_bayes,
              Lambda = Lambda_st,
              #acp = acp/(nrun-burn),
              tau = tau_st,
              sigmasq_st = sigmasq_st,
              a = a))
}

# rotate factors to enforce identifiability
# first method based on varimax rotation 
# second method from BADFM (Aßmann, Boysen- Hogrefe, and Pape 2014)

# ARGUMENTS: lambda: file path to factor matrix sample list (or just the list);
#            method: rotation method; one of c("varimax", "BADFM");
#            tolerance: rotation algorithm stopping tolerance;
#            maxiter: maximum number of algorithm iterations;
#            ncores: number of cores
#            normalize: logical. whether to normalize lambda samples
#            file: logical. whether lambda was passed directly or as an Rds

mcrotfact = function(lambda, method = "BADFM", 
                     tolerance = 1e-5, maxiter = 100, ncores = 1,
                     normalize = FALSE, file = TRUE){

  if(file) lambda = readRDS(lambdafile) # read in factor samples
  n = length(lambda)           # initialize shared attributes
  if(normalize) lambda = mclapply(lambda, scale,
                                  mc.cores = ncores)
  
  if(method == "varimax"){        # Varimax rotations
    Lr = mclapply(lambda,
                  function(lam) as(varimax(lam, eps = tolerance)[["loadings"]],
                                   "matrix"), 
                  mc.cores = ncores)
    LrMean = Reduce("+", Lr) / n
    class(LrMean) = "matrix"
    return(list(samples = Lr, mean = LrMean))
  }
  
  if(method == "BADFM"){          # BADFM rotations
    iter = 0
    err = tolerance + 1
    oldrot = lambda[[n]]
    
    while(err > tolerance & iter < maxiter){
      iter = iter + 1
      
      Lr = mclapply(lambda, function(lam){
        Sr = t(lam) %*% oldrot
        Ssvd = La.svd(Sr)
        D = Ssvd[["u"]] %*% Ssvd[["vt"]]
        return(lam %*% D)}, mc.cores = ncores)
      
      newrot = Reduce("+", Lr) / n
      err = norm(oldrot - newrot, type = "F")
      oldrot = newrot
    }
    return(list(samples = Lr, mean = newrot, iter = iter))
  }
}

# Solve Label Switching by Clustering
# Default method is multi-pivot
# Other method for probabalistic clustering

# ARGUMENTS: lambda: list of factor samples with constant dimension
#            method: one of c("multipivot","prob") [currently only multipivot supported]
#            intitialpivot: pivot sample (or matrix more generally)
#            stop: largest reasonable diff norm for an aligned cluster mean (see permfact)
#            itermax: max permutation index (see permfact)

clustalignplus = function(lambda, method = "multipivot", initialpivot = NULL, 
                          stop = NULL, itermax = 100000){
  
  if(is.null(initialpivot)){
    initialpivot = sample(lambda, 1)[[1]]
    pivot = initialpivot
  } else {
    pivot = sample(lambda, 1)[[1]]
  }                                                      # allow initial pivot value, or recursion
  
  if(length(lambda) < 5) return(permsignfact(lambda = lambda, 
                                             pivot = initialpivot, 
                                             stop = stop, 
                                             itermax = itermax))
  
  difflist = lapply(lambda, `-`, pivot)
  norms = sapply(difflist, norm, type = "2")             # norm differences between samples and a pivot (svd)
  norms[which.min(norms)] = min(norms[norms>1e-4])       # deal with pivot norm of 0
  
  clusters = kmeans(norms, 2, nstart = 1)$cluster        # cluster samples based on norm 
  
  baseCluster = which.min(c(var(norms[clusters==1]),     # find tightest cluster
                            var(norms[clusters==2])))
  if(is.null(stop)) {
    lower = which.min(c(mean(norms[clusters==1]),        # find base sampling noise for stopping criterion
                        mean(norms[clusters==2])))
    stop = mean(norms[clusters==lower]) + sd(norms[clusters==lower])
  }
  
  if(var(norms)/4 <= var(norms[clusters==baseCluster])){  # if only one cluster likely, align
    return(permsignfact(lambda = lambda, pivot = initialpivot, 
                        stop = stop, itermax = itermax))
  } else {                                               # if more than one cluster, reapply clustalign
    
    c1 = which(clusters == 1)                            # always align to the initial pivot
    lambda1 = lambda[c1]                                 # original sampling noise transmitted
    aligned1 = clustalignplus(lambda1, initialpivot = initialpivot, 
                              stop = stop, itermax = itermax)
    lambda[c1] = aligned1
    
    c2 = which(clusters == 2)
    lambda2 = lambda[c2]
    aligned2 = clustalignplus(lambda2, initialpivot = initialpivot, 
                              stop = stop, itermax = itermax)
    lambda[c2] = aligned2
    
    return(lambda)                                       # return collated, aligned samples
  }
}

# perform sign changes to minimise a norm

# ARGUMENTS: pivot: a pivot
#            permed: a permuted matrix
#            stop: the stopping criterion

signer = function(pivot, permed, stop, step = 1, 
                  start = rep(1, ncol(pivot)), 
                  minnorm = Inf, min = start){
  
  norm1 = norm(pivot - permed, "2")
  if(step == 1) minnom = norm1
  if(norm1 < stop) return(start)
  
  w = which.max(colSums((pivot - permed)^2))
  switch = rep(1, ncol(pivot))
  switch[w] = -1
  switched = t(t(permed) * switch)
  norm2 = norm(pivot - switched, "2")
  start = start * switch
  
  if(norm2 < minnorm) {
    minnorm = norm2
    min = start
  }
  
  if(step <= ncol(pivot)) return(signer(pivot, switched, stop, 
                                        step = step + 1, 
                                        start = start,
                                        minnorm = minnorm,
                                        min = min))
  if(step > ncol(pivot)) return(min)
}

# function to permute in minimal-switch order
# will output repeated orders for i > k!

# ARGUMENTS: vec: a vector to permute
#            i: index of permutation in minimal-switch order

permuter = function (vec, i) {
  k = length(vec)
  j = i %% (k^2) %/% k
  l = i %% k 
  if(! j) j = k
  if(! l) l = k
  temp = vec[l]
  vec[l] = vec[j]
  vec[j] = temp
  if(i %/% (k^2)) {
    return(permuter(vec, i %/% (k^2)))
  } else {return(vec)}
}

# ARGUMENTS: lambda: list of factor samples
#            pivot: matrix to align permutation to
#            stop: stopping criterion, largest reasonable norm for an aligned cluster mean
#            itermax: maximum number of permutations to search, can be larger than vector memory

permsignfact = function(lambda, pivot, stop, itermax = 100000){
  k = ncol(lambda[[1]])
  p = nrow(lambda[[1]])
  m = length(lambda)
  lamarray = unlist(lambda) %>% array(c(p, k, m))
  first = Reduce("+", lambda)/length(lambda)             # align median of cluster samples to the pivot
  k = ncol(first)
  
  i = 1
  mindiff = Inf
  minperm = 0L
  while(i < itermax){
    perm = permuter(1:k, i)                             # iteratively search 
    sign = signer(pivot, first[,perm], stop)
    signed = t(t(first[,perm]) * sign)
    diff = norm(pivot - signed, type = "2")
    if(diff < stop) break
    if(diff < mindiff){
      mindiff = diff
      minperm = perm
      minsign = sign
    }
    i = i + 1
  }
  
  if(i == itermax) {
    print(paste("permsignfact itermax of", i, "reached"))
    print(minperm * minsign)
    aligned = lapply(lambda, function(mat) t(t(mat[,minperm]) * minsign))
    return(aligned)
  }
  
  aligned = lapply(lambda, function(mat) t(t(mat[,perm]) * sign))
  print("cluster successfully permuted")
  print(perm * sign)
  return(aligned)
}







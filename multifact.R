# Multiple factor perturbation model
# See Isaac 2019

# ARGUMENTS: Y: list of data matrices, one per study, must have same, aligned columns
#            iter: number of total MCMC iterations
#            burn: number of burn in samples to discard
#            

multifact <- function(Y, iter = 1000, burn = 500, # Y is a list of data matrices from the different studies
                         prop = 1,                # Proportion of redundant elements before elimination of unnecessary column
                         epsilon = 1E-2,           # Cutoff for redundant element
                      phip = 2
){
                                # hyperparameters for:
  as = 1; bs = 0.3              # residual precision
  df = 3                        # local shrinkage factor
  dfp = 3                       # local shrinkage factor for the perturbation loadings matrix
  ad1 = 2; bd1 = 1              # delta_1
  ad2 = 6; bd2 = 1              # delta_h, h >= 2
  adp1 = 2; bdp1 = 1            # delta_1 for the perturbation loadings matrix
  adp2 = 6; bdp2 = 1            # delta_h, h >= 2 for the perturbation loadings matrix
  adf = 1; bdf = 1              # ad1 and ad2 or df
  b0 = 1; b1 = 0.0005           # Adaptation
  
  if(class(Y) != "list"){
    stop("Y should be a list of data matrices from multiple studies")
  }
  
  # --- data attributes --- #
  S = length(Y)
  n = sapply(Y, nrow)
  p = ncol(Y[[1]])
  if(all(p == sapply(Y, ncol))) stop("inconsistent predictor dimension")
  
  # ¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿¿Scale the data????????????????
  
  # --- init shared parameters --- #
  k = floor(log(p)*3)
  sigma = 1/rgamma(p, as,bs)     # residual noise after accounting for latent factors
  
  delta = c(rgamma(1, ad1, bd1), rgamma(k-1, ad2, bd2))
  tauh = cumprod(delta) 
  
  phi = matrix(rgamma(p*k, df/2, df/2), nrow = p, ncol = k)   
  
  Lambda = matrix(rnorm(p*k, 0, sqrt(1/t(t(phi)*tauh))),
                  nrow = p, ncol = k)
  
  # --- init study-specific params --- #
  eta = lapply(1:S, function(s) matrix(rnorm(n[s]*k),
                                       nrow = n[s], ncol = k))
  
  delta_perturb = lapply(1:S, function(s) c(rgamma(1, adp1, bdp1),
                                            rgamma(k-1, adp2, bdp2)))
  tauh_perturb = lapply(delta_perturb, cumprod)
  
  phi_perturb = lapply(1:S, function(s) matrix(rgamma(p*k, dfp/2, dfp/(2*phip)),
                                               nrow = p, ncol = k))
  
  Psi_perturb = mapply(function(phi, tauh) matrix(rnorm(p*k, 0, sqrt(1/t(t(phi)*tauh))),
                                                  nrow = p, ncol = k),
                       phi_perturb, tauh_perturb, SIMPLIFY = FALSE)
  ind = 0
  storage_variables = c("sigma", "eta", "lambda", "delta", "phi", "tauh", "omega", "k", "loadings", "psi")
  stores = lapply(storage_variables, function(v) vector("list", length=iter-burn))
  names(stores) = storage_variables
  
  for(i in 1:iter){
    
    # --- update eta --- #
    
    loadings_list = lapply(Psi_perturb, `+`, Lambda)
    Prec_list = lapply(loadings_list, function(loading) diag(k) + 
                         crossprod(loading * 1/sigma, loading))
    b_list = mapply(function(Y, loading) crossprod(loading * 1/mc$sigma, t(Y)),
                    Y, loadings_list, SIMPLIFY = FALSE)
    
    eta = mapply(function(n, b, prec){
      L = chol(prec)
      mean = solve(prec, b)
      t(mean + solve(L, matrix(rnorm(n*k), nrow=k, ncol=n))) # try backsolve
    }, n, b_list, Prec_list)
    
    # --- update lambda --- #
    
    b_cols = lapply(1:p, function(j) Reduce(`+`,
      mapply(function(e, w, y) crossprod(e, (y[,j] - e %*% w[j,])),
             eta, Psi_perturb, Y)))
    
    eta_t_eta = Reduce(`+`, lapply(eta, crossprod))
    
    phitau = t(t(phi) * tauh)
    phitau_rows = lapply(1:p, function(i) phitau[i,])
    
    Lambda = t(mapply(function(b, phitau, sigma){
      D = diag(phitau)
      prec = D + 1/sigma * eta_t_eta
      sds = sqrt(1/diag(prec))
      near0 = sds < 1E-3
      if(any(near0)){
        L = chol(prec[!near0, !near0])
        mean = solve(prec[!near0, !near0], 1/sigma * b[!near0])
        c(mean + solve(L, matrix(rnorm(sum(!near0)), 
                                 nrow=sum(!near0))), rep(0, sum(near0)))
      }else{
        L = chol(prec)
        mean = solve(prec, 1/sigma * b)
        mean + solve(L, matrix(rnorm(k), nrow=k)) # again try backsolve here
      }
      
    }, b_cols, phitau_rows, sigma))
    
    # --- update psi --- #
    
    b_cols_list = mapply(function(e, y) lapply(1:p, 
                           function(j) crossprod(e, (y[,j] - e %*% Lambda[j,]))),
                         eta, Y, SIMPLIFY = FALSE)
    
    eta_t_eta_list = lapply(eta, crossprod)
    
    phitau = mapply(function(phi, tauh) t(t(phi) * tauh), 
                    phi_perturb, tauh_perturb, SIMPLIFY = FALSE)
    
    phitau_rows_list = lapply(phitau, function(phitau_s) 
      lapply(1:p, function(i) phitau_s[i,]))
    
    Psi_perturb = mapply(function(b_cols, phitau_rows, eta_t_eta, eta) 
      t(mapply(function(b, phitau, sigma) {      
        D = diag(phitau)
        prec = D + 1/sigma * eta_t_eta
        sds = sqrt(1/diag(prec))
        near0 = sds < 1E-3
        if(any(near0)){
          L = chol(prec[!near0, !near0])
          mean = solve(prec[!near0, !near0], 1/sigma * b[!near0])
          c(mean + solve(L, matrix(rnorm(sum(!near0)), 
                                   nrow=sum(!near0))), rep(0, sum(near0)))
        }else{
          L = chol(prec)
          mean = solve(prec, 1/sigma * b)
          mean + solve(L, matrix(rnorm(k), nrow=k)) # again try backsolve here
        }}, b_cols, phitau_rows, sigma)),
      b_cols_list, phitau_rows_list, eta_t_eta_list, eta, SIMPLIFY = FALSE)
    
    # --- update phi --- #
    
    phi = t(apply(Lambda, 1, function(L){
      rgamma(k, (df + 1)/2, (df + L^2*tauh)/2)}))
    
    # --- update perturbed phi --- #
    
    phi_perturb = mapply(function(Psi, tauh) t(apply(Psi, 1,
                           function(w) rgamma(k, (dfp+1)/2, (dfp/phip + w^2*tauh)/2))),
                         Psi_perturb, tauh_perturb, SIMPLIFY = FALSE)
    
    # --- update delta and tau_h --- #
    
    Lambdasq_phi_sum = apply(Lambda^2 * phi, 2, sum)
    delta[1] = rgamma(1, ad1 + p*k/2, bd1 +
                        (1/2) * (1/delta[1]) * sum(tauh*Lambdasq_phi_sum))
    tauh = cumprod(delta)
    
    for(h in 2:k){
      delta[h] = rgamma(1, ad2 + p*(k - h + 1)/2, bd2 +
                          (1/2)*(1/delta[h]) * sum(tauh[h:k] * Lambdasq_phi_sum[h:k]))
      tauh = cumprod(delta)
    }
    
    # --- update perturbed delta and tau_h --- #
    Psisq_phi_sum_list = lapply(Psi_perturb, function(Psi) apply(Psi^2 * phi, 2, sum))
    delta_perturb = mapply(function(delta, tauh, Psisq_phi_sum){
      delta[1] = rgamma(1, adp1 + p*k/2, bdp1 + (1/2) * (1/delta[1]) * sum(tauh*Psisq_phi_sum))
      tauh = cumprod(delta)
      
      for(h in 2:k){
        delta[h] = rgamma(1, adp2 + p*(k - h + 1)/2, bdp2 * (1/delta[h-1]) + 
                            (1/2)*(1/delta[h]) * sum(tauh[h:k] * Psisq_phi_sum[h:k]))
        tauh = cumprod(delta)
      }
      return(delta)
    }, delta_perturb, tauh_perturb,
    Psisq_phi_sum_list, SIMPLIFY = FALSE)
    tauh_perturb = lapply(delta_perturb, cumprod)
    
    # --- update sigma --- #
    
    loadings_list = lapply(Psi_perturb, `+`, Lambda)
    err_list = mapply(function(Y, loading, eta) Y - tcrossprod(eta, loading),
                    Y, loadings_list, eta, SIMPLIFY = FALSE)
    err_sum = Reduce(`+`, lapply(err_list, function(err) apply(abs(err^2), 2, sum))) #F norm?
    
    sigma = 1/rgamma(p, as + sum(n)/2, bs + err_sum/2)
    
    # --- adapt k if necessary --- #
    
    prob = 1/exp(b0 + b1*i)
    u = runif(1)
    proportion_redundant = apply(abs(Lambda) < epsilon, 2, sum) / p
    redundant = proportion_redundant >= prop
    nredundant = sum(redundant)
    
    if( u < prob){
      if((i > 20) & (nredundant == 0) & all(proportion_redundant < 0.995)){
        k = k + 1
        Lambda = cbind(Lambda, 0)
        Psi_perturb = lapply(Psi_perturb, function(Psi) cbind(Psi, 0))
        eta = mapply(function(eta, n) cbind(eta, rnorm(n)), 
                     eta, n, SIMPLIFY = FALSE)
        phi = cbind(phi, rgamma(p, df/2, df/2))
        phi_perturb = lapply(phi_perturb, 
                             function(phi) cbind(phi, rgamma(p, dfp/2, dfp/2)))
        delta = c(delta, rgamma(1, ad2, bd2))
        delta_perturb = lapply(delta_perturb, 
                               function(delta) c(delta, rgamma(1, adp2, bdp2)))
        tauh = cumprod(delta)
        tauh_perturb = lapply(delta_perturb, cumprod)
      }
      else if(nredundant > 0){
        k = max(1, k - nredundant)
        keep = !redundant
        Lambda = Lambda[,keep]
        Psi_perturb = lapply(Psi_perturb, function(Psi) Psi[,keep])
        phi = phi[,keep]
        phi_perturb = lapply(phi_perturb, function(phi) phi[,keep])
        eta = lapply(eta, function(eta) eta[,keep])
        delta = delta[keep]
        delta_perturb = lapply(delta_perturb, function(delta) delta[keep])
        tauh = cumprod(delta)
        tauh_perturb = lapply(delta_perturb, cumprod)
      }
    }
    
    # Save samples
    if(i >= burn){
      ind = ind + 1
      stores$lambda[[ind]] = Lambda
      stores$sigma[[ind]] = sigma
      stores$eta[[ind]] = eta
      stores$phi[[ind]] = phi
      stores$delta[[ind]] = delta
      stores$tauh[[ind]] = tauh
      stores$psi[[ind]] =  Psi_perturb
      stores$loadings[[ind]] =  lapply(Psi_perturb, `+`, Lambda)
      stores$omega[[ind]] = Lambda %*% t(Lambda) + diag(sigma)
      stores$k[[ind]] = k
    }
    
    if((i %% 1000) == 0) {
      print(i)
    }
    
  }
  
  return(stores)
}

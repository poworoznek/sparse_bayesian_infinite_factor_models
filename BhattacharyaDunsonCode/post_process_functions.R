# Post-processing functions

# Implied covariance matrix
get_omega <- function(stores){
  tmp = abind(stores$omega_store, along=3)
  round(apply(tmp, c(1,2), mean), 2)
}

get_sigma <- function(stores){
  tmp = abind(map(stores$sigma_store, function(s) diag(s)), along=3)
  round(apply(tmp, c(1,2), mean), 2)
}


# Using the Orthogonal Procrustes approach to recover the Lambda factor loadings matrix
# Following all the steps from the De Vite and Parmigiani paper (at the moment)
# See the steps followed on github: https://github.com/rdevito/MSFA

get_lambda <- function(stores){
  get_lambda_bmsfa(stores$lambda_store)
}

get_lambda_bmsfa <- function(lambda_samps){
  # Get the estimated (median) covariance matrix from the lambda alone
  Lambda_cov = abind(map(lambda_samps, function(L) tcrossprod(L)), along=3)
  Lambda_cov = apply(Lambda_cov, c(1,2), median)
  
  # To choose a value of k, take the normalized eigenvalue decomposition, see how many eigenvalues are above a prespecified level
  k = sum(sp_eigen(Lambda_cov) > 0.05)
  
  # Now need to get samples of the matrix at this level of k...
  
  # In the vignette, they do something werid, keeping just the first k columns of all the samples
  # Not sure if that's correct, but I will do the same thing
  # If k is less, will just fill in 0's
  p = nrow(lambda_samps[[1]])
  
  add_k_cols <- function(x, p, k){
    old_k = ncol(x)
    if(k > old_k){
      for(i in (old_k+1):k){
        x = cbind(x, 0)
      }
    }
    x
  }
  
  tmp = lambda_samps %>%
    modify_if(function(x) ncol(x) < k, function(x) add_k_cols(x, p, k))  %>%
    modify_if(function(x) ncol(x) > k, function(x) matrix(x[,1:k], ncol=k))
  
  lambda_samps = abind(tmp, along=3)
  
  round(sp_OP(lambda_samps)$Phi, 2)
}
  
get_k <- function(stores){
  table(as.numeric(stores$k_store))
}

lambda_cor <- function(lambda, true_lambda){
  k = min(ncol(lambda), ncol(true_lambda))
  map_dbl(1:k, function(c) cor(lambda[,c], true_lambda[,c]))
}

# Similarity between 2 square matrices
RV <- function(m1, m2){
  tr <- function(m){ sum(diag(m))}
  tr(tcrossprod(m1, m2) %*% tcrossprod(m1, m2)) / (tr(tcrossprod(m1))^2 * tr(tcrossprod(m2))^2)
}

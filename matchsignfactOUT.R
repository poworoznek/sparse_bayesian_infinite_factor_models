# Find optimal column permutation and signs
# called from clustalign?
# quadratic loss columnwise

# ARGUMENTS: lambda: list of factor samples or single matrix sample
#            pivot: matrix with which to align permutation and signs

matchsignfactOUT = function(lambda, pivot){
  
  refr = cbind(lambda, -lambda)
  k = ncol(lambda)
  ind = c(1:k, -(1:k))
  permsign = c()
  
  for(i in 1:k){
    norms = apply((refr-pivot[,i])^2, 2, sum)
    w = which.min(norms)
    permsign[i] = ind[w]
    c = ncol(refr)/2
    if(w <= c) {r = c(w, w+c)} else {r = c(w, w-c)}
    ind = ind[-r]
    refr = refr[,-r]
  }
  
  return(permsign)
}


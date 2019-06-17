# Find optimal column permutation and signs
# called from clustalign?
# quadratic loss columnwise

# ARGUMENTS: lambda: list of factor samples or single matrix sample
#            pivot: matrix with which to align permutation and signs

matchsignfact = function(lambda, pivot){
  if(class(lambda)=="list"){
    repr = Reduce("+", lambda)/length(lambda)    # lambda representative
  } else {
    repr = lambda                                # or actual lambda if a single sample
  }
  refr = cbind(repr, -repr)
  k = ncol(repr)
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
  
  if(class(lambda)=="list"){
    aligned = lapply(lambda, function(mat) {t(t(mat[,abs(permsign)]) * sign(permsign))})
  } else {
    aligned = t(t(lambda[,abs(permsign)]) * sign(permsign))
  }
  return(aligned)
}


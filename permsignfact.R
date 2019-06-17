# Find optimal permutations
# called from clustalign
# svd loss

# ARGUMENTS: lambda: list of factor samples
#            pivot: matrix to align permutation to
#            stop: stopping criterion, largest reasonable norm for an aligned cluster mean
#            itermax: maximum number of permutations to search, can be larger than vector memory

permsignfact = function(lambda, pivot, stop, itermax = 100000){
  k = ncol(lambda[[1]])
  p = nrow(lambda[[1]])
  m = length(lambda)
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


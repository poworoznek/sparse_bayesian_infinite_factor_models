# Find optimal permeutations
# called from clustalign
# svd loss?

# ARGUMENTS: lambda: list of factor samples
#            pivot: matrix to align permutation to

permfact = function(lambda, pivot, stop, itermax = 1000000){
  first = sample(lambda, 1)[[1]]
  k = ncol(first)
  
  i = 1
  while(i < itermax){
    perm = permuter(1:k, i)
    diff = norm(pivot - first[,perm], type = "2")
    if(diff < stop) break
    i = i + 1
  }
  
  aligned = lapply(lambda, function(mat) mat[,perm])
  print("done")
  return(aligned)
}

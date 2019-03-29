# Find optimal permutations
# called from clustalign
# svd loss

# ARGUMENTS: lambda: list of factor samples
#            pivot: matrix to align permutation to
#            stop: stopping criterion, largest reasonable norm for an aligned cluster mean
#            itermax: maximum number of permutations to search, can be larger than vector memory

permfact = function(lambda, pivot, stop, itermax = 100000){
  k = ncol(lambda[[1]])
  p = nrow(lambda[[1]])
  m = length(lambda)
  lamarray = unlist(lambda) %>% array(c(p, k, m))
  first = apply(lamarray, c(1, 2), median)              # align median of cluster samples to the pivot
  k = ncol(first)
  
  i = 1
  mindiff = Inf
  minperm = 0L
  while(i < itermax){
    perm = permuter(1:k, i)                             # iteratively search 
    diff = norm(pivot - first[,perm], type = "2")
    if(diff < stop) break
    if(diff < mindiff){
      mindiff = diff
      minperm = perm
    }
    i = i + 1
  }
  
  if(i == itermax) {
    print(paste("permfact itermax of", i, "reached"))
    print(minperm)
    aligned = lapply(lambda, function(mat) mat[,minperm])
    return(aligned)
  }
  
  aligned = lapply(lambda, function(mat) mat[,perm])
  print("cluster successfully permuted")
  print(perm)
  return(aligned)
}


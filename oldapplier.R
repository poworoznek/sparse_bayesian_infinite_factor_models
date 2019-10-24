# Apply CAPout output to a list

# ARGUMENTS: l: list of matrices to apply to permutations to
#            p: permutation / sign output from  CAPout

applier = function(l, p){
  return(mapply(function(mat, perm){
    return(t(t(mat[,abs(perm)]) * sign(perm)))
  }, l, p, SIMPLIFY = FALSE))
}
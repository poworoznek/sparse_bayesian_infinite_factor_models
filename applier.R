# Apply CAPout output to a list

# ARGUMENTS: l: list of matrices to apply to permutations to
#            p: permutation / sign output from  CAPout

applier = function(l, p){
  return(mapply(aplr, l, p, SIMPLIFY = FALSE))
}

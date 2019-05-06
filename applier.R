# Apply CAPout output to a list

# ARGUMENTS: l: list of vectors to apply to permutations to
#            p: permutation / sign output from  CAPout

applier = function(l, p){
  return(mapply(function(vec, perm){
    return(vec[abs(perm)] * sign(perm))
  }, l, p))
}
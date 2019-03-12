# function to permute in minimal-switch order
# will output repeated orders for i > k!

# ARGUMENTS: vec: a vector to permute
#            i: index of permutation in minimal-switch order

permuter = function (vec, i) {
  k = length(vec)
  j = i %% (k^2) %/% k
  l = i %% k 
  if(! j) j = k
  if(! l) l = k
  temp = vec[l]
  vec[l] = vec[j]
  vec[j] = temp
  if(i %/% (k^2)) {
    return(permuter(vec, i %/% (k^2)))
  } else {return(vec)}
}

# perform sign changes to minimise a norm

# ARGUMENTS: pivot: a pivot
#            permed: a permuted matrix
#            stop: the stopping criterion

signer = function(pivot, permed, stop, step = 1, 
                  start = rep(1, ncol(pivot)), 
                  minnorm = Inf, min = start){
  
  norm1 = norm(pivot - permed, "2")
  if(step == 1) minnom = norm1
  if(norm1 < stop) return(start)

  w = which.max(colSums((pivot - permed)^2))
  switch = rep(1, ncol(pivot))
  switch[w] = -1
  switched = t(t(permed) * switch)
  norm2 = norm(pivot - switched, "2")
  start = start * switch
  
  if(norm2 < minnorm) {
    minnorm = norm2
    min = start
  }
  
  if(step <= ncol(pivot)) return(signer(pivot, switched, stop, 
                                        step = step + 1, 
                                        start = start,
                                        minnorm = minnorm,
                                        min = min))
  if(step > ncol(pivot)) return(min)
}


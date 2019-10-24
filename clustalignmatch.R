# Solve Label Switching by Clustering
# Default method is multi-pivot
# Other method for probabalistic clustering

# ARGUMENTS: lambda: list of factor samples with constant dimension
#            method: one of c("multipivot","prob") [currently only multipivot supported]
#            intitialpivot: pivot sample (or matrix more generally)
#            stop: largest reasonable diff norm for an aligned cluster mean (see permfact)
#            itermax: max permutation index (see permfact)

clustalignmatch = function(lambda, initialpivot = NULL){
  
  if(is.null(initialpivot)){
    initialpivot = sample(lambda, 1)[[1]]
    pivot = initialpivot
  } else {
    pivot = sample(lambda, 1)[[1]]
  }                                                      # allow initial pivot value, or recursion
  
  if(length(lambda) < 5) return(matchsignfact(lambda = lambda, 
                                             pivot = initialpivot))
  
  difflist = lapply(lambda, `-`, pivot)
  norms = sapply(difflist, norm, type = "2")             # norm differences between samples and a pivot (svd)
  norms[which.min(norms)] = min(norms[-(which.min(norms))])       # deal with pivot norm of 0
  
  clusters = kmeans(norms, 2, nstart = 1)$cluster        # cluster samples based on norm 
  
  baseCluster = which.min(c(var(norms[clusters==1]),     # find tightest cluster
                            var(norms[clusters==2])))
  
  if(var(norms)/3 <= var(norms[clusters==baseCluster])){  # if only one cluster likely, align
    return(matchsignfact(lambda = lambda, pivot = initialpivot))
  } else {                                               # if more than one cluster, reapply clustalign
    
    c1 = which(clusters == 1)                            # always align to the initial pivot
    lambda1 = lambda[c1]                                 # original sampling noise transmitted
    aligned1 = clustalignmatch(lambda1, initialpivot = initialpivot)
    lambda[c1] = aligned1
    
    c2 = which(clusters == 2)
    lambda2 = lambda[c2]
    aligned2 = clustalignmatch(lambda2, initialpivot = initialpivot)
    lambda[c2] = aligned2
    
    return(lambda)                                       # return collated, aligned samples
  }
}





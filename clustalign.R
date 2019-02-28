# Solve Label Switching by Clustering
# Default method is multi-pivot
# Other method for probabalistic clustering

# ARGUMENTS: lambda: list of factor samples with constant dimension
#            method: one of c("multipivot","prob")
#            intitialpivot: pivot sample (or matrix more generally)

clustalign = function(lambda, method = "multipivot", initialpivot = NULL){
  
  if(is.null(initialpivot)){
    initialpivot = sample(lambda, 1)[[1]]
    pivot = initialpivot
  } else {
    pivot = sample(lambda, 1)[[1]]
  }                                                      # allow initial pivot value, or recursion
  
  difflist = lapply(lambda, `-`, pivot)
  norms = sapply(difflist, norm, type = "2")             # norm differences between samples and a pivot (svd)
  norms[which.min(norms)] = min(norms[norms>1e-4])       # deal with pivot norm of 0
  
  clusters = kmeans(norms, 2, nstart = 1)$cluster        # cluster samples based on norm 
  
  baseCluster = which.min(c(var(norms[clusters==1]),     # find tightest clsuter
                            var(norms[clusters==2])))
  
  if(var(norms)/4 < var(norms[clusters==baseCluster])){  # if only one cluster likely, align
    return(permfact(lambda, initialpivot, 1))
  } else {                                               # if more than one cluster, reapply clustalign
    
    c1 = which(clusters == 1)
    lambda1 = lambda[c1]
    aligned1 = clustalign(lambda1, initialpivot = initialpivot)
    lambda[c1] = aligned1
    
    c2 = which(clusters == 2)
    lambda2 = lambda[c2]
    aligned2 = clustalign(lambda2, initialpivot = initialpivot)
    lambda[c2] = aligned2
    # return collated, aligned samples
    return(lambda)
  }
}





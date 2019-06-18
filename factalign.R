# Solve rotational ambiguity in factor loadings
# Variety of methods

# ARGUMENTS: lambda: List of samples
#            rotation: Method for rotation: one of c("varimax", "BADFM")
#            match: Logical. Use greedy matching?
#            cluster: Logical. Use greedy clustering?
#            verbose: Logical. Print updates?
#            ncores: How many threads to run in parallel
#            rerun: How many times to iterate clustering

factalign = function(lambda, rotation = "varimax", match = TRUE,
                     cluster = FALSE, verbose = TRUE, ncores = 1L,
                     rerun = 100L, itermax = 1000L){
  if(class(lambda) != "list"){if(class(lamdba)=="character"){lambda = readRDS(lambda)}
    else {stop("inappropriate lambda type")}}
  
  if(verbose) print("rotating")
  rotated = mcrotfact(lambda = lambda, method = rotation,
                      ncores = ncores, file = FALSE)
  rm(lambda)
  
  if(verbose) print("aligning")
  if(match){
    if(cluster){
      aligned = clustalignmatch(rotated$samples)
      if(rerun > 1) for(i in 2:rerun) aligned = clustalignmatch(aligned)
    } else {
      aligned = mclapply(rotated$samples, matchsignfact, 
                         rotated$samples[[round(length(rotated$samples)/2)]],
                         mc.cores = ncores, mc.preschedule = TRUE)
    }
  } else {
    if(cluster){
      prealign = lapply(rotated$samples[1:5], matchsignfact, 
                        rotated$samples[[1]])
      diff = lapply(prealign, "-", rotated$samples[[1]])
      stop = max(sapply(diff, norm, "2"))
      aligned = clustalignplus(lambda = rotated$samples, stop = stop, 
                               itermax = itermax)
      if(rerun > 1) for(i in 2:rerun) aligned = clustalignplus(lambda = aligned, 
                                                               stop = stop, 
                                                               itermax = itermax)
    } else {
      stop("must use clustering or matching")
    }
  }
  return(aligned)
}

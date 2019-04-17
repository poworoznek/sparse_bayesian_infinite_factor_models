# rotate factors to enforce identifiability
# first method based on varimax rotation 
# second method from BADFM (Aßmann, Boysen- Hogrefe, and Pape 2014)

# ARGUMENTS: lambda: file path to factor matrix sample list (or just the list);
#            method: rotation method; one of c("varimax", "BADFM");
#            tolerance: rotation algorithm stopping tolerance;
#            maxiter: maximum number of algorithm iterations;
#            ncores: number of cores
#            normalize: logical. whether to normalize lambda samples
#            file: logical. whether lambda was passed directly or as an Rds

mcrotfact = function(lambda, method = "BADFM", 
                     tolerance = 1e-5, maxiter = 100, ncores = 1,
                     normalize = FALSE, file = TRUE){
  library(parallel)
  if(file) lambda = readRDS(lambdafile) # read in factor samples
  n = length(lambda)           # initialize shared attributes
  if(normalize) lambda = mclapply(lambda, scale,
                                  mc.cores = ncores)
  
  if(method == "varimax"){        # Varimax rotations
    Lr = mclapply(lambda,
                  function(lam) as(varimax(lam, eps = tolerance)[["loadings"]],
                                   "matrix"), 
                  mc.cores = ncores)
    LrMean = Reduce("+", Lr) / n
    class(LrMean) = "matrix"
    return(list(samples = Lr, mean = LrMean))
  }
  
  if(method == "BADFM"){          # BADFM rotations
    iter = 0
    err = tolerance + 1
    oldrot = lambda[[n]]
    
    while(err > tolerance & iter < maxiter){
      iter = iter + 1
      
      Lr = mclapply(lambda, function(lam){
        Sr = t(lam) %*% oldrot
        Ssvd = La.svd(Sr)
        D = Ssvd[["u"]] %*% Ssvd[["vt"]]
        return(lam %*% D)}, mc.cores = ncores)
      
      newrot = Reduce("+", Lr) / n
      err = norm(oldrot - newrot, type = "F")
      oldrot = newrot
    }
    return(list(samples = Lr, mean = newrot, iter = iter))
  }
}


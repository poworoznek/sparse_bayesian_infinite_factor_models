# rotate factors to enforce identifiability
# first method based on varimax rotation 
# second method from BADFM (Aßmann, Boysen- Hogrefe, and Pape 2014)

# ARGUMENTS: lambdafile: file path to factor matrix sample list;
#            method: rotation method; one of c("varimax", "BADFM");
#            tolerance: rotation algorithm stopping tolerance;
#            maxiter: maximum number of algorithm iterations;

mcrotfact = function(lambdafile, method = "BADFM", 
                     tolerance = 1e-5, maxiter = 100, ncores = 1){
  library(parallel)
  lambda = readRDS(lambdafile) # read in factor samples
  n = length(lambda)           # initialize shared attributes
  
  if(method == "varimax"){        # Varimax rotations
    Lr = mclapply(lambda,
                  function(lam) varimax(lam, eps = tolerance)[["loadings"]], 
                  mc.cores = ncores)
    LrMean = Reduce("+", Lr) / n
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


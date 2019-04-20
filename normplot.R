# function to plot the norms of differences of
# factor loadings matrix samples

# ARGUMENTS: lambda: list of lambda samples

normplot = function(lambda){
  pivot = base::sample(lambda, 1)[[1]]
  difflist = lapply(lambda, "-", pivot)
  norms = sapply(difflist, norm, "2")
  norms = norms[-which.min(norms)]
  hist(norms, breaks = 100)
}

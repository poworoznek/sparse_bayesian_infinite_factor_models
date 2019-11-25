# only display  a matrix cell if posterior credible interval DNI 0
source("lmean.R")

summat = function(list){
  d1 = nrow(list[[1]])
  d2 = ncol(list[[2]])
  d3 = length(list)
  ar = array(unlist(list), dim = c(d1, d2, d3))
  low = apply(ar, c(1, 2), quantile, 0.025)
  high = apply(ar, c(1, 2), quantile, 0.975)
  cond = (low < 0) & (high > 0)
  m = lmean(list)
  m[cond] = 0
  return(m)
}
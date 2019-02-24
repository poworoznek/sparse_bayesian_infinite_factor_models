source("fixedfact.R")
library(mvtnorm)
library(tidyverse)
library(abind)

setwd("~/documents/factorR/")
source("./simulate_data_fxns.R")
source("./BhattacharyaDunsonCode/post_process_functions.R")
source("./BhattacharyaDunsonCode/mcmc.R")

sigmasq = 1
p = 100
n = 100
k = 10

output = simulate_x(n, p, k, sigmasq, simulate_lambda, return_lambda = TRUE)
Y = output$x

fixedfact(Y, nrun = 10000, burn = 2000, k = 10, output = "factSamples", factfilename = "factorsForLS.rds")

lambda = readRDS("factorsForLS.rds")
lamarray = unlist(lambda) %>% array(c(100, 10, 8000))
lammed = apply(lamarray, c(1, 2), median)

#### Test for label switches

difflist = lapply(lambda, function(L) L - lammed)
norms = sapply(difflist, norm, type = "f")
plot(density(norms))

w = which(norms > 6)
w2 = which(norms > 4.5 & norms < 6)

for(i in 1:10){
plot(lamarray[20, i, ], pch = ".")
points(w, lamarray[20, i, w], pch = ".", col = "red")
points(w2, lamarray[20, i, w2], pch = ".", col = "blue")
}

#### Force labels to switch on artificial data
lamvar = apply(lamarray, c(1, 2), var)

lambda.fake = lapply(1:8000, function(x) {lammed + 
    matrix(rnorm(1000, 0, mean(lamvar)), ncol = 10)})

difflist.fake = lapply(lambda.fake, function(L) L - lammed)
norms.fake = sapply(difflist.fake, norm, type = "2")
plot(density(norms.fake))

# perform switching

lambda.switch = lapply(1:8000, function(ind){
  if(ind < 2500 | ind > 4500){ 
    lambda.fake[[ind]]
  } else {
      lambda.fake[[ind]][, c(1:8, 10, 9)]}
})

difflist.switch = lapply(lambda.switch, function(L) L - lammed)
norms.switch = sapply(difflist.switch, norm, type = "2")
plot(density(norms.switch))

lamarray.switch = unlist(lambda.switch) %>% array(c(100, 10, 8000))

w =  which(norms.switch > 2)

for(i in 9:10){
  plot(lamarray.switch[20, i, ], pch = ".")
  points(w, lamarray.switch[20, i, w], pch = ".", col = "red")
}
















lambda2 = readRDS("Lambda2.rds")
lamarray2 = unlist(lambda2) %>% array(c(100, 10, 8000))
lammed2 = apply(lamarray2, c(1, 2), median)

difflist2 = lapply(lambda2, function(L) L - lammed2)
norms2 = sapply(difflist2, norm, type = "f")
plot(density(norms2))


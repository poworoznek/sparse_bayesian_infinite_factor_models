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

fixedfact(Y, nrun = 10000, burn = 200, k = 10, output = "factSamples", factfilename = "factorsForLS.rds")

lambda = readRDS("factorsForLS.rds")
lamarray = unlist(lambda) %>% array(c(100, 10, 8000))
lammed = apply(lamarray, c(1, 2), median)

#### Test for label switches

difflist = lapply(lambda, `-`, lammed)
norms = sapply(difflist, norm, type = "2")
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

par(mfrow = c(2,2))

w =  which(norms.switch > 2)

for(i in 9:10){
  plot(lamarray.switch[20, i, ], pch = ".",
       main = paste("column", i, "row 20"))
  points(w, lamarray.switch[20, i, w], pch = ".", col = "red")
}



lambda.aligned = clustalign2(lambda.switch)
difflist.aligned = lapply(lambda.aligned, function(L) L - lammed)
norms.aligned = sapply(difflist.aligned, norm, type = "2")
plot(density(norms.aligned))

lamarray.aligned = unlist(lambda.aligned) %>% array(c(100, 10, 8000))

w =  which(norms.aligned > 2)

for(i in 9:10){
  plot(lamarray.aligned[20, i, ], pch = ".",
       main = paste("column", i, "row 20 clustaligned"))
  points(w, lamarray.aligned[20, i, w], pch = ".", 
         col = "red", main = paste("column", i, "row 20"))
}


lambda.switch2 = lapply(1:8000, function(ind){
  if(ind < 4000 | ind > 7000){ 
    lambda.switch[[ind]]
  } else {
    lambda.switch[[ind]][, c(2, 1, 3:10)]}
})

difflist.switch2 = lapply(lambda.switch2, function(L) L - lammed)
norms.switch2 = sapply(difflist.switch2, norm, type = "2")
plot(density(norms.switch2))

lambda.aligned2 = clustalign(lambda.switch2)
difflist.aligned2 = lapply(lambda.aligned2, function(L) L - lammed)
norms.aligned2 = sapply(difflist.aligned2, norm, type = "2")
plot(density(norms.aligned2))


lamarray.aligned2 = unlist(lambda.aligned2) %>% array(c(100, 10, 8000))

w =  which(norms.aligned2 > 2)

for(i in 9:10){
  plot(lamarray.aligned2[20, i, ], pch = ".",
       main = paste("column", i, "row 20 clustaligned"))
  points(w, lamarray.aligned2[20, i, w], pch = ".", 
         col = "red", main = paste("column", i, "row 20"))
}


lambda = readRDS("factorsForLS.rds")
lamarray = unlist(lambda) %>% array(c(100, 10, 8000))
lammed = apply(lamarray, c(1, 2), median)

#### True test for label switches
source("clustalign.R")
source("permfact.R")
source("permuter.R")

difflist = lapply(lambda, function(L) L - lammed)
norms = sapply(difflist, norm, type = "2")
plot(density(norms))
lambda.test = clustalign(lambda)

difflist.test = lapply(lambda.test, function(L) L - lammed)
norms.test = sapply(difflist.test, norm, type = "2")
plot(density(norms.test))


# Burn in that looks like label switching
lambda2 = readRDS("Lambda2.rds")
lamarray2 = unlist(lambda2) %>% array(c(100, 10, 8000))
lammed2 = apply(lamarray2, c(1, 2), median)

difflist2 = lapply(lambda2, `-`, lammed2)
norms2 = sapply(difflist2, norm, type = "f")
plot(density(norms2))

lam2.test = clustalign(lambda2)

difflist2.test = lapply(lam2.test, `-`, lammed2)
norms2.test = sapply(difflist2.test, norm, type = "f")
plot(density(norms2.test))



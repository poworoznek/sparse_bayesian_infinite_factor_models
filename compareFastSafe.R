# compare fastfact and safefact
library(parallel)
library(mvtnorm)
library(abind)
source("fastfact.R")
source("safefact.R")
source("simulate_data_fxns.R")

samp = function(numfact){
  sigmasq = 1; p = 100; n = 100; ktrue = numfact
  output = simulate_x(n, p, ktrue, sigmasq, simulate_lambda, 
                      return_lambda = TRUE)
  Y = output$x
  
  safek = safefact(Y, nrun = 10000, burn = 0, output = "numFactors")
  fastk = fastfact(Y, nrun = 10000, burn = 0, output = "numFactors")
  
  return(list(safek = safek, fastk = fastk, k = numfact))
}

ks = round(runif(20, 5, 30))

test = mclapply(ks, samp, mc.cores = 20, mc.preschedule = FALSE)

saveRDS(test, "testFastSave.rds")

par(mfrow = c(2, 3))

plot(test[[c(2, 1, 1)]], type = "l")
lines(test[[c(2, 2, 1)]], col = "red")
abline(h = test[[c(2, 3)]], col = "blue", lwd = 2)

plot(test[[c(5, 1, 1)]], type = "l")
lines(test[[c(5, 2, 1)]], col = "red")
abline(h = test[[c(5, 3)]], col = "blue", lwd = 2)
legend(x = 6000, y = 10, col = c("black", "red"), 
       lty = c(1, 1), legend = c("safefact", "fastfact"))

plot(test[[c(6, 1, 1)]], type = "l")
lines(test[[c(6, 2, 1)]], col = "red")
abline(h = test[[c(6, 3)]], col = "blue", lwd = 2)

plot(test[[c(7, 1, 1)]], type = "l")
lines(test[[c(7, 2, 1)]], col = "red")
abline(h = test[[c(7, 3)]], col = "blue", lwd = 2)

plot(test[[c(11, 1, 1)]], type = "l")
lines(test[[c(11, 2, 1)]], col = "red")
abline(h = test[[c(11, 3)]], col = "blue", lwd = 2)

plot(test[[c(13, 1, 1)]], type = "l", ylim = c(4, 14))
lines(test[[c(13, 2, 1)]], col = "red")
abline(h = test[[c(13, 3)]], col = "blue", lwd = 2)

source("fastfact2.R")

ks2 = round(runif(200, 5, 30))

samp2 = function(numfact){
  sigmasq = 1; p = 100; n = 100; ktrue = numfact
  output = simulate_x(n, p, ktrue, sigmasq, simulate_lambda, 
                      return_lambda = TRUE)
  Y = output$x
  
  safek = safefact(Y, nrun = 10000, burn = 0, output = "numFactors")
  fastk = fastfact2(Y, nrun = 10000, burn = 0, output = "numFactors")
  
  return(list(safek = safek, fastk = fastk, k = numfact))
}

test2 = mclapply(ks2, samp2, mc.cores = 16, mc.preschedule = FALSE)

saveRDS(test2, "Fast2Safe200runs.rds")

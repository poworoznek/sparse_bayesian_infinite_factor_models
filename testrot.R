# test rotations
library(mvtnorm)
library(magick)
source("fixedfact.R")
source("mcrotfact.R")
source("spclone.R")

set.seed(13)

vals = rep(c(rep(1, 15), rep(0, 97)), 10)[1:1000]

lam = matrix(vals, nrow = 100)
lam[,10] = runif(100) * (-1)^(round(runif(100)))

image(t(lam))

cov = lam %*% t(lam) + diag(100)

dat = rmvnorm(1000, mean = rep(0, 100), sigma = cov)

fixedfact(Y = dat, nrun = 10000, burn = 2000, k = 10, 
          output = "factSamples", factfilename = "Lambda2.rds")

mcrotated = mcrotfact(lambdafile = "Lambda2.rds", method = "BADFM", maxiter = 100, ncores = 4, tol = 1e-5)

image(t(mcrotated$mean))

library(label.switching)
lambda = array(unlist(mcrotated$samples), dim = c(100,10,8000))
lambda2 = aperm(lambda, c(3, 2, 1))
image(lambda2[7500,,])
# perm = pra(mcmc.pars = lambda2, pivot = lambda2[7500,,])
difflist = mapply(`-`, mcrotated$samples, mcrotated$samples[[7500]], SIMPLIFY = FALSE)
norm = sapply(difflist, norm, type = "f")
plot(density(norm))
w = which(norm>10)
plot(c(1:8000)[-w], lambda[1, 1, -w], pch = ".")
points(w, lambda[1, 1, w], pch =".", col = "red")

for(i in 1:10){
  plot(c(1:8000)[-w], lambda[(i-1)*10+1, i, -w], pch = ".")
  points(w, lambda[(i-1)*10+1, i, w], pch =".", col = "red")
}


fact = readRDS("Lambda2.rds")
factl = array(unlist(fact), dim = c(100,10,8000))
lambda.image =  image_scale(image_read(abs(factl[,,1,  drop = F] / max(factl))), "600x600!")
for(i in seq(2, 8000,  by = 100)){
lambda.add =  image_scale(image_read(abs(factl[,,i,  drop = F] / max(factl))), "600x600!")
lambda.image =  c(lambda.image, lambda.add)
}
lambda.animated = image_animate(lambda.image, fps = 20, dispose = "background")
image_write(lambda.animated, "LambdaAnimation.gif")


rot.image = image_scale(image_read(abs(lambda[,,1,  drop = F] / max(lambda))), "600x600!")
for(i in seq(2, 8000,  by = 100)){
  rot.add =  image_scale(image_read(abs(lambda[,,i,  drop = F] / max(lambda))), "600x600!")
  rot.image =  c(rot.image, rot.add)
}
rot.animated = image_animate(rot.image, fps = 20, dispose = "background")
image_write(rot.animated, "rot.gif")


fmean = Reduce("+", fact)

image(t(fmean))

sprotated = spclone(lambdafile = "Lambda2.rds", method = "BADFM", maxiter = 100, ncores = 4, tol = 1e-5)

image(t(sprotated$mean))

image(t(mcrotated$mean))
image(t(sprotated$mean))

norm(mcrotated$mean - sprotated$mean, type = "F")

vmrotated = mcrotfact(lambdafile = "Lambda2.rds", method = "varimax", maxiter = 100, ncores = 4, tol = 1e-5)

image(t(vmrotated$mean))
image(t(fmean))

lambda = array(unlist(vmrotated$samples), dim = c(100,10,8000))
plot(lambda[10,10,], pch = ".")

frames = image_graph(width = 600, height = 600)

walk(1000:1100, ~{
  print(plotmat(lambda_bayes[[.x]], "green"))
})

dev.off()

ttry = image_animate(frames, fps = 5)

image_write(ttry, "raw.gif")

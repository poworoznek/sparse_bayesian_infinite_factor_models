# swiss roll ex
library(plot3D)

phi = runif(1000, 1.5 * pi, 4.5 * pi)
psi = runif(1000, 1, 10)

x = phi * cos(phi)
y = phi * sin(phi)
z = psi

points3D(x, y, z, colvar = phi)

kern = function(X, Y, c){
  return(exp(-c*sum((X-Y)^2)))
}

dat = cbind(x, y, z)

k = apply(dat, 1, function(X){
  apply(dat, 1, kern, X, 0.1)
})

range(k)

eig = eigen(k)

lambda = lm(phi ~ eig$vectors[,1:3])

summary(lambda)

wrong = lm(phi ~ x + y + z)

summary(wrong)




library(tidyverse)
library(bayesSurv)
library(GIGrvg)
library(statmod)
library(MCMCpack)
library(reshape2)
library(ggplot2)
library(latex2exp)
library(psych)
library(parallel)
source("https://raw.githubusercontent.com/poworoznek/sparse_bayesian_infinite_factor_models/master/demo_source.R")
source("plotmat.R")
set.seed(1)
k0 = 9
p = 200

lambda = matrix(rnorm(p*k0, 0, 0.01), ncol = k0)
lambda[sample.int(p, 90, replace = TRUE) +
         p*(sample.int(k0, 90, replace = TRUE)-1)] = rnorm(90, 0, 1)
lambda[1:40, 1] = rnorm(40, 2, 0.5)
lambda[41:60, 2] = rnorm(20, -2, 0.5)
lambda[61:100, 3] = rnorm(40, 2, 0.5)
lambda[101:170, 4] = rnorm(70, 2, 0.5)
lambda[171:200, 5] = rnorm(30, -2, 0.5)
lambda[,6] = rnorm(p, 0, 0.5)
lambda[,7] = rnorm(p, 0, 0.5)
lambda[,8] = rnorm(p, 0, 0.5)
lambda[,9] = rnorm(p, 0, 0.5)

lvar = varimax(lambda)$loadings

class(lvar)="matrix"

plotmat(lvar, color = "green")

n = 500;
X = matrix(rnorm(n*k0),n,k0)%*%t(lambda) + bayesSurv::rMVNorm(n,Sigma = 3*diag(p))

##### SAMPLER

##### Factor model for X ####
library(bayesSurv)
library(GIGrvg)
library(statmod)
library(MCMCpack)
library(psych)

#### Generate X ### 

k = 9
nrun = 20000; burn = 10000; thin = 1
sp = floor((nrun - burn)/thin)

a = 1/2
b0 = 1            
b1 = 0.0005
epsilon = 1e-3
prop = 1.00
as = 1
bs = 0.3
df = 150
ad1 = 3
bd1 = 1
ad2 = 4
bd2 = 1
adf = 1
bdf = 1

Omega_bayes = array(0,c(sp,p,p))  # sample storage memory allocation
lambda_bayes=list()
Sigma_st = matrix(0,sp,p)
tau_st = array(0,c(sp,p))
sigmasq_st = numeric(sp)

ps = rgamma(p,as,bs)
Sigma = diag(1/ps)
Lambda = matrix(0,p,k)
eta = matrix(rnorm(n*k),n,k)
meta = matrix(0,n,k)
veta = diag(rep(1,k))

tau = rgamma(p,p*a,1/2)
psijh = matrix(rexp(p*k,1/2),p,k)
phijh = matrix(0,p,k)
for(j in 1:p){
  phijh[j,] = rdirichlet(1,rep(a,k))
}
Plam = psijh*(phijh^2)*matrix(rep(tau^2,k),p,k,byrow=F)

t = t0 = Sys.time()               # begin sample timing
count = 1 

for(i in 1:nrun){
  
  # --- Update eta --- #
  Lmsg = Lambda * ps
  Veta1 = diag(k) + t(Lmsg) %*% Lambda
  Tmat = chol(Veta1)
  S = solve(Tmat)
  Veta = tcrossprod(S)                                               # Veta = inv(Veta1)
  Meta = X %*% Lmsg %*% Veta                                      # n x k 
  eta = Meta + tcrossprod(matrix(rnorm(n*k), nrow = n, ncol = k), S)  
  eta.T = t(eta)
  
  # --- Update Lambda --- #
  Plam = psijh*(phijh^2)*matrix(rep(tau^2,k),p,k,byrow=F)
  eta2 = eta.T%*%eta
  zlams = rnorm(k*p)       # generate normal draws all at once 
  
  for(j in 1:p) {
    Llamt = chol(diag(Plam[j,]) + ps[j]*eta2)
    Lambda[j,] = t(backsolve(Llamt,
                             zlams[1:k + (j-1)*k]) + 
                     backsolve(Llamt,
                               forwardsolve(t(Llamt),
                                            ps[j] * eta.T %*% X[,j])))
  }
  Lambda.T = t(Lambda)
  
  # --- Update psijh --- #
  mujh = phijh*matrix(rep(tau,k),p,k,byrow = F)/abs(Lambda)
  psijh = matrix(rinvgauss(k*p, mujh),p,k,byrow = T)
  
  # --- Update tau --- #
  for(j in 1:p){
    tau[j] = GIGrvg::rgig(n=1,lambda = 1-k, psi = 1, 
                          chi = 2*sum(abs(Lambda[j,])/phijh[j,]))
  }
  
  # --- Update phijh --- #
  Tjh = numeric(k)
  for(j in 1:p){
    for(h in 1:k){
      Tjh[h] = GIGrvg::rgig(n=1,lambda = a-1,psi = 1,
                            chi = 2*abs(Lambda[j,h]))
    }
    phijh[j,] = Tjh/sum(Tjh)
  }
  
  # --- Update Sigma --- #
  Xtil = X - eta%*%Lambda.T
  ps = rgamma(p,as+0.5*n,1)
  ps = (1/(bs + 0.5*colSums(Xtil^2)))*ps
  Sigma = diag(1/ps)
  
  #store posterior samples
  if ((i %% thin == 0) & (i > burn)){
    
    # parameters of the factor model
    tau_st[count,] = tau
    #Omega_bayes[count,,] = Lambda%*%t(Lambda) + diag(1/ps)
    lambda_bayes[[count]] = Lambda
    Sigma_st[count,] = 1/ps
    count = count + 1
  }
  
  if (i%%100==0){
    print(i)
  }
}

lsampvar = varimax(lambda_bayes[[length(lambda_bayes)]])$loadings
class(lsampvar) = "matrix"
plotmat(lsampvar)

lsampvar2 = varimax(lambda_bayes[[length(lambda_bayes)-1000]])$loadings
class(lsampvar2) = "matrix"
plotmat(lsampvar2)

plotmat(lambda_bayes[[length(lambda_bayes)]])
plotmat(lambda_bayes[[length(lambda_bayes)-1000]])

sample_mean = Reduce(`+`, lambda_bayes)/length(lambda_bayes)
plotmat(sample_mean, color = "green")
rotated = mcrotfact(lambda_bayes, method = "varimax", file = FALSE, ncores = 4)
rotmean = Reduce(`+`, rotated$samples)/length(rotated$samples)
plotmat(rotmean, color = "green")
aligned = clustalignplus(rotated$samples, itermax = 10000, stop = 0)

for(i in 1:100) aligned = clustalignplus(aligned, itermax = 10000, stop = 0, initialpivot = lambda)
plotmat(Reduce(`+`, aligned)/length(aligned), color = "green")

library(gridExtra)


grid.arrange(plotmat(sample_mean, color = "green"), plotmat(rotmean, color = "green"),
             plotmat(Reduce(`+`, aligned)/length(aligned), color = "green"), ncol = 3)

grid.arrange(plotmat(lvar, color = "green"), 
             plotmat(rotated$samples[[1200]], color = "green"),
             plotmat(rotated$samples[[9000]], color = "green"),
             ncol = 3)

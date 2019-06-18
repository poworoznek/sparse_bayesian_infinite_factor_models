##### Factor model for X ####
library(bayesSurv)
library(GIGrvg)
library(statmod)
library(MCMCpack)
library(psych)

#### Generate X ### 
p = 10; 
W = riwish(p+1,diag(p))
n = 50000; 
svd_W = svd(W)
d = svd_W$d; u = svd_W$u; v = svd_W$v;
X = bayesSurv::rMVNorm(n=n,Sigma=W)
eig = eigen(W)

k = 9
d_hat = c(d[1:k],rep(0,p-k))
W_hat = u%*%diag(d_hat)%*%t(v)

nrun = 1000; burn = 500; thin = 1
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
    Omega_bayes[count,,] = Lambda%*%t(Lambda) + diag(1/ps)
    Sigma_st[count,] = 1/ps
    count = count + 1
  }
  
  if (i%%100==0){
    print(i)
  }
}

Omega_hat = apply(Omega_bayes, c(2,3), mean)
err = mean(Omega_hat-W)
err_sq = mean((Omega_hat-W)^2)
err_svd = mean(Omega_hat-W_hat)
err_svd_sq = mean((Omega_hat-W_hat)^2)
#print(c(err,err_sq,err_svd,err_svd_sq))
mean(W-W_hat);mean((W-W_hat)^2)

# Frob
Fr_true = Omega_hat - W; sqrt(tr(t(Fr_true)%*%Fr_true))
#Fr_svd = Omega_hat - W_hat; sqrt(tr(t(Fr_svd)%*%Fr_svd))
Fr_svd_true = W - W_hat; sqrt(tr(t(Fr_svd_true)%*%Fr_svd_true))

# svd from the Cov of (X)
C = cov(X); svd_C = svd(C); d_c = svd_C$d; u_c = svd_C$u; v_c = svd_C$v;
d_c_hat = c(d_c[1:k],rep(0,p-k))
W_c_hat = u_c%*%diag(d_c_hat)%*%t(v_c)
#Fr_c = C - W; sqrt(tr(t(Fr_c)%*%Fr_c))
Fr_c_hat = W_c_hat - W; sqrt(tr(t(Fr_c_hat)%*%Fr_c_hat))



#apply(X,2,var);apply(Sigma_st,2,mean)

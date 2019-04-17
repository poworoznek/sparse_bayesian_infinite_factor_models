# function to perform factorized regression with interactions
# See Federico's paper (2019?)
# Dirichlet Laplace prior with random walk met-hastings

# ARGUMENTS: y: vector of responses;
#            X: Matrix of predictors
#            delta_rw: MH tuning parameter
#            epsilon_rw:  MH tuning parameter
#            a: Dirichlet Laplace prior hyperparameter
#            k: number of factors

library(bayesSurv)
library(psych)
library(GIGrvg)
library(statmod)
library(MCMCpack)

gibbs_DL = function(y, X ,nrun, burn, thin = 1, 
                    delta_rw = 0.0526749, epsilon_rw = 0.5,
                    a = 1/2, k = NULL){
   n = nrow(X)
   p = ncol(X)                    # collect data attributes
   
   if(is.null(k)) k = floor(log(p)*3)
   if(length(y) != n) stop("Mismatching input lengths")
   
   VX = apply(X, 2, var)          # prepare data matrix
   X = as.matrix(scale(X))
   
   sp = floor((nrun - burn)/thin)
   
   b0 = 1            # factorization hyperparameters
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
   tau_st = array(0,c(sp,p))
   sigmasq_st = numeric(sp)
   alpha_bayes = numeric(sp)
   beta_bayes = matrix(0,sp,p)
   Lambda_st = array(0,c(sp,p,k))
   acp = numeric(n)
   
   sigmasq_y = 1                     # initial values
   phi = numeric(k)
   ps = rgamma(p,as,bs)
   Sigma = diag(1/ps)
   Lambda = matrix(0,p,k)
   Lambda.T = t(Lambda)
   eta = matrix(rnorm(n*k),n,k)
   meta = matrix(0,n,k)
   veta = diag(rep(1,k))
   Psi = matrix(0,k,k)
   
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
      aMH = phi%*%t(phi)/sigmasq_y + Lambda.T%*%diag(ps)%*%Lambda + diag(k)
      
      for (h in 1:n){                # Metropolis hastings step 
         
         eta_star = bayesSurv::rMVNorm(1,eta[h,],diag(k)*delta_rw)
         eta_star.T = t(eta_star)     # avoid repeated transpose calls
         eta.T = t(eta[h,])
         
         logr = eta_star.T%*%(aMH - 2*Psi*y[h]/sigmasq_y)%*%eta_star -
            2*eta_star.T%*%(Lambda.T%*%diag(ps)%*%X[h,] + phi*y[h]/sigmasq_y) +
            2*eta_star.T%*%phi*(eta_star%*%Psi%*%eta_star)/sigmasq_y + 
            (1/sigmasq_y)*(eta_star%*%Psi%*%eta_star)^2 - 
            (eta.T%*%(aMH - 2*Psi*y[h]/sigmasq_y)%*%eta[h,] -
                2*eta.T%*%(Lambda.T%*%diag(ps)%*%X[h,] + phi*y[h]/sigmasq_y) +
                2*eta.T%*%phi*(eta[h,]%*%Psi%*%eta[h,])/sigmasq_y +
                (1/sigmasq_y)*(eta[h,]%*%Psi%*%eta[h,])^2)
         logr = logr*(-0.5)
         
         logu = log(runif(1))
         
         if (logr > logu){
            eta[h,] = eta_star
            acp[h] = acp[h] + 1
         }
      }
      
      # --- Update Psi --- #
      MM = model.matrix(y~.^2 - 1,as.data.frame(eta))   # perform factorized regression
      X_reg = cbind(eta^2,MM[,(k+1):ncol(MM)])
      X_reg.T = t(X_reg)
      Lambda_n = X_reg.T%*%X_reg/sigmasq_y + diag(rep(1,ncol(X_reg)))
      Vcsi = solve(Lambda_n)
      Mcsi = Vcsi%*%X_reg.T%*%(y-eta%*%phi)/sigmasq_y
      csi = bayesSurv::rMVNorm(n=1,mean=Mcsi,Sigma=Vcsi)
      lambda_diag = csi[1:k]
      lambda = csi[(k+1):length(csi)]
      Psi[lower.tri(Psi)] = lambda/2
      Psi[upper.tri(Psi)] = 0
      Psi = Psi + t(Psi)
      diag(Psi) = lambda_diag
      
      # --- Update phi --- #
      eta.T = t(eta)
      Lambda_n = eta.T%*%eta/sigmasq_y + diag(rep(1,ncol(eta)))
      Vcsi = solve(Lambda_n)
      Mcsi = Vcsi%*%eta.T%*%(y-diag(eta%*%Psi%*%eta.T))/sigmasq_y     # using updated psi
      phi = bayesSurv::rMVNorm(n = 1, mean = Mcsi, Sigma = Vcsi)
      
      # --- Update sigmasq_y --- #
      an = 0.5 + n/2
      bn = 0.5 + 0.5*t(y-eta%*%phi-diag(eta%*%Psi%*%eta.T))%*%
         (y-eta%*%phi-diag(eta%*%Psi%*%eta.T))
      sigmasq_y = 1/rgamma(1,an,bn)
      
      # --- Update Lambda --- #
      Plam = psijh*(phijh^2)*matrix(rep(tau^2,k),p,k,byrow=F)
      eta2 = eta.T%*%eta
      zlams = rnorm(k*p, sd = 3)       # generate normal draws all at once 
      
      for(j in 1:p) {
         Llamt = chol(diag(Plam[j,]) + ps[j]*eta2)
         Lambda[j,] = t(solve(Llamt,
                              zlams[1:k + (j-1)*k]) + 
                           solve(Llamt,
                                 solve(t(Llamt),
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
      ps = (1/(bs + 0.5*apply(Xtil^2,2,sum)))*ps
      Sigma = diag(1/ps)
      
      #store posterior samples
      if ((i %% thin == 0) & (i > burn)){
         
         # parameters of the factor model
         sigmasq_st[count] = sigmasq_y
         tau_st[count,] = tau
         
         #Bayesian estimators
         V_n = solve(Lambda.T%*%solve(Sigma)%*%Lambda+diag(rep(1,ncol(Lambda))))
         a_n = V_n%*%Lambda.T%*%solve(Sigma)
         dsVX = diag(sqrt(VX))
         Omega_bayes[count,,] = dsVX%*%t(a_n)%*%Psi%*%a_n%*%dsVX
         beta_bayes[count,] = as.vector(t(phi)%*%a_n)
         alpha_bayes[count] = tr(Psi%*%V_n)
         Lambda_st[count,,] = dsVX %*% Lambda
         count = count + 1
      }
      
      #if (i%%100==0){
      #   print(i)
      #   acp_mean = mean(acp)/100
      #   print(acp_mean)
      #   if(acp_mean > 0.3){
      #      delta_rw = delta_rw*2
      #   }else if(acp_mean < 0.2){
      #      delta_rw = delta_rw*2/3
      #   }
      #   acp = numeric(n)
      #   print(delta_rw)
      #   #print(paste("time for last 100 iterations:",round(as.numeric(Sys.time()-t),0),
      #   #            "seconds",sep=" "))
      #   #t = Sys.time()
      #   #print(paste(i,"out of",nrun,sep=" "))
      #   #t_end = round(((as.numeric(difftime(Sys.time(),t0,units="secs")))/i)*(nrun-i),0)
      #   #print(paste("estimated time to end:",t_end,"seconds",sep=" "))
      #}
   }
   
   beta_bayes = beta_bayes%*%diag(sqrt(VX))
   return(list(alpha_bayes = alpha_bayes,
               beta_bayes = beta_bayes,
               Omega_bayes = Omega_bayes,
               Lambda = Lambda_st,
               acp = acp/(nrun-burn),
               tau = tau_st,
               sigmasq_st = sigmasq_st,
               a = a))
}

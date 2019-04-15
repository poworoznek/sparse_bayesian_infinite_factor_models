generate_DL_reg = function(S = 10000, k = 15, p = 20,
                           a_DL = 0.5, var_omega = 5,
                           as = 1, bs = 0.3,
                           trace = FALSE, every = 500){

   library(MCMCpack)
   
   # Storage
   # third order interactions in the model
   beta1_X_DL =  matrix(0,p,S)
   beta2_X_DL = array(0,c(p,p,S))
   beta3_X_DL = array(0,c(1,2,3,S))

   
   for (i in 1:S){
      
      Lambda_DL = matrix(0,p,k)
      ps = rgamma(p,as,bs)
      for(j in 1:p){
         phi_DL = rdirichlet(1, rep(a_DL,p))
         for(h in 1:k){
            
            #Dirichlet-Laplace
            psi = rexp(1,0.5)
            tau_DL = rgamma(1,k*a_DL,0.5)
            Lambda_DL[j,h] = rnorm(1,0,tau_DL*phi_DL[j]*sqrt(psi))
         }
      }
      
      omega1 = rnorm(k,0,var_omega)
      omega2 = rnorm(k,0,var_omega)
      omega3 = rnorm(k,0,var_omega)
      A_DL = solve(t(Lambda_DL)%*%diag(ps)%*%Lambda_DL + 
                      diag(k))%*%t(Lambda_DL)%*%diag(ps)

      
      
      # up to 3rd order interactions
      beta1_X_DL[,i] = t(A_DL)%*%omega1 + 3*diag(ps)%*%t(A_DL)%*%omega3

      beta2_X_DL[,,i] = 2*(2*t(A_DL)%*%diag(omega2)%*%A_DL + 
                              6*diag(1/ps)%*%t(A_DL)%*%diag(omega2)%*%A_DL)

      beta3_X_DL[1,2,3,i] = 6*sum(omega3*A_DL[,1]*A_DL[,2]*A_DL[,3])

      if(trace){
         if (i%%every==0){
            print(i)
         }
      }
      
   }
   
   
   xlim = c(-3,3)
   par(mfrow = c(1,3))
   hist(beta1_X_DL[1,],freq = F,breaks = seq(-1000,1000,by =0.05),xlim = xlim,
        xlab = "main effect",main = "")
   abline(v = quantile(beta1_X_DL[1,],probs = c(0.05,0.95)),col="red",lty="dotted")
   abline(v = quantile(beta1_X_DL[1,],probs = c(0.25,0.75)),col="green",lty="dotted")
   hist(beta2_X_DL[1,2,],freq = F,breaks = seq(-1000,1000,by =0.05),xlim = xlim,
        xlab = "2nd order interaction",main =  paste("k =",k,",","p =",p))
   abline(v = quantile(beta2_X_DL[1,2,],probs = c(0.05,0.95)),col="red",lty="dotted")
   abline(v = quantile(beta2_X_DL[1,2,],probs = c(0.25,0.75)),col="green",lty="dotted")
   hist(beta2_X_DL[1,2,],freq = F,breaks = seq(-1000,1000,by =0.05),xlim = xlim,
        xlab = "3rd order interaction",main = "")
   abline(v = quantile(beta3_X_DL[1,2,3,],probs = c(0.05,0.95)),col="red",lty="dotted")
   abline(v = quantile(beta3_X_DL[1,2,3,],probs = c(0.25,0.75)),col="green",lty="dotted")
   
   
}


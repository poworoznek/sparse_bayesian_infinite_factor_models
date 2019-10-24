####### CLEAN CODE for simulations ######
library(mvtnorm)
library(MASS)
library(beepr)
library(psych)
library(bayesSurv)
library('PIE')
library('glmnet')
library(RAMP)
library(hierNet)
library(FAMILY)
library(RCurl)
library(stargazer)
library(R.utils)
sourceDirectory("~/Desktop/Code_factor_simulations/Functions")

#p = 10, k = 8, delta = 0.015

S = 50
err = err_pred = FR = err_beta = matrix(0,ncol = 6, nrow = S)
alpha = 0.05
TP_main = TN_main = TP_int = TN_int = matrix(0,ncol = 6, nrow = S)
acp_min = acp_max = acp_mean = matrix(0,ncol = 2,nrow = S)

for(s in 1:S){
   
   # generate the data
   #data = generate_factor_model(p = 10,n = 100,k_true = 8)
   data = generate_indep_model(p = 10,n = 100)
   y = data$y; X = data$X; beta_true = data$beta_true; 
   Omega_true = data$Omega_true;
   y_test = data$y_test; X_test = data$X_test
   p = ncol(X)
   
   #Factor models
   nrun = 2500; burn = 1500;
   gibbs_DLk = gibbs_factor_dirichlet_laplace(y, X ,nrun = nrun, burn = burn, thin = 1, 
                                                delta_rw = 0.08, epsilon_rw = 0.5,
                                                a = k, k = floor(log(ncol(X))*3))
   gibbs_DLp = gibbs_factor_dirichlet_laplace(y, X ,nrun = nrun, burn = burn, thin = 1, 
                                             delta_rw = 0.065, epsilon_rw = 0.5,
                                             a = p, k = floor(log(ncol(X))*3))

   
   
   #plot_acp(nrun,burn,gibbs_DLk)
   #plot_acp(nrun,burn,gibbs_DLp)
   
   
   
   #acp 
   acp_1 = gibbs_DLk$acp;acp_2 = gibbs_DLp$acp
   acp_min[s,1] = min(acp_1/(nrun-burn));acp_max[s,1] = max(acp_1/(nrun-burn));acp_mean[s,1] = mean(acp_1/(nrun-burn))
   acp_min[s,2] = min(acp_2/(nrun-burn));acp_max[s,2] = max(acp_2/(nrun-burn));acp_mean[s,2] = mean(acp_2/(nrun-burn))
   
   #Competitors
   hiernet = Hiernet_fct(y, X, X_test, y_test)
   Family = FAMILY_fct(y, X, X_test, y_test)
   PIE = PIE_fct(y, X, X_test, y_test)
   RAMP = RAMP_fct(y, X, X_test, y_test)
   
   #Errors
   errors = compute_errors(hiernet,Family,PIE,RAMP,
                           y,y_test,Omega_true,beta_true,
                           gibbs_DLk,gibbs_DLp)
   
   #Error 
   err[s,] = errors$err_insample[-1]
   err_pred[s,] = errors$err_pred[-1]
   FR[s,] = errors$FR_norm[-1]
   err_beta[s,] = errors$beta_MSE[-1]
   
   #TP and TNs
   rate_main = rate_recovery_maineff(gibbs_DLk,gibbs_DLp,alpha = alpha,beta_true = beta_true,
                                      hiernet$beta,Family$beta,PIE$beta,RAMP$beta)
   rate_int = rate_recovery_interactions(gibbs_DLk,gibbs_DLp,alpha = alpha,Omega_true=Omega_true,
                                          hiernet$Omega,Family$Omega,PIE$Omega,RAMP$Omega)
   
   TP_main[s,] = rate_main$TP
   TN_main[s,] = rate_main$TN
   TP_int[s,] = rate_int$TP
   TN_int[s,] = rate_int$TN
   
   print(s)
}

colnames(err_beta) = colnames(err_pred) = colnames(FR) = colnames(err) =
   c("Hiernet","Family","Pie","RAMP","BFM_DLk","BFM_DLp")
colnames(TP_main) = colnames(TN_main) = colnames(TP_int) = colnames(TN_int) = 
   c("BFM_DLk","BFM_DLp","Hiernet","Family","Pie","RAMP")

apply(err,2,mean)
apply(err_pred,2,mean)
apply(err_beta,2,mean)
apply(FR,2,mean)

apply(TP_main,2,mean)
apply(TN_main,2,mean)
apply(TP_int,2,mean)
apply(TN_int,2,mean)

apply(err,2,sd)
min_err = apply(err,1, function(x) (x == min(x)) )
min_err = apply(min_err,1,sum)

apply(err_pred,2,mean)
apply(err_pred,2,sd)
min_err_pred = apply(err_pred,1, function(x) (x == min(x)) )
min_err_pred = apply(min_err_pred,1,sum)

apply(FR,2,mean)
apply(FR,2,sd)
min_FR = apply(FR,1, function(x) (x == min(x)) )
min_FR = apply(min_FR,1,sum)

apply(err_beta,2,mean)
apply(err_beta,2,sd)
min_beta = apply(err_beta,1, function(x) (x == min(x)) )
min_beta = apply(min_beta,1,sum)

G = rbind(apply(err,2,mean),min_err/100,
          apply(err_pred,2,mean),min_err_pred/100,
          apply(FR,2,mean),min_FR/100,
          apply(err_beta,2,mean),min_beta/100)
rownames(G) = c("MSE","","MSE prediction","",
                "Frobenious norm","","MSE beta","")

stargazer(G)




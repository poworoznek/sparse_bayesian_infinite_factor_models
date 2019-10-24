library(tidyverse)
library(mvtnorm)
library(abind)
source("multifact.R")
source("mcmc_perturb_model.R")
source("simulate_data_fxns 2.R")

n=100
p=50
k=10
D=10
sigmasq = 1

out = simulate_x_perturbed(n=n, p=p, k=k, D=D, 
                           sigmasq=sigmasq, 
                           lambda_generator = simulate_lambda())


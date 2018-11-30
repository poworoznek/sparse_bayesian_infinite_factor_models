# Simulating data from a factor analysis
# These latent factor matrices are the same as used in the simulated example of Expandable Factor Analysis

library(mvtnorm)
library(tidyverse)

sigmasq = 1
p = 20
n = 10
k = 5

x = simulate_x(n, p, k, sigmasq, simulate_lambda_EFA_example)

corrplot(cov2cor(cov(x)))

x = simulate_x(n, p, k = 3, sigmasq, simulate_lambda)

corrplot(cov2cor(cov(x)))

# Simulate data with 2 extra, study specific, latent variables and add it to the phthalate data
p = 9
n = 382
k = 2
sigmasq = .5


# Testing out the simulated part of the data
output = simulate_x(n, p, k, sigmasq, simulate_lambda, return_lambda = TRUE)
x = output$x
lambda = output$lambda

corrplot(cov2cor(cov(x)))

# Standardize the Mt. Sinai Phthalate Data
standardize <- function(x){ (x - mean(x)) / sd(x)}
data = apply(phth_data, 2, standardize)

# Create datasets by adding in simulated data to the Mt. Sinai data
add_simulated_latent_factors <- function(data, k, sigmasq, lambda_generator){
  n = nrow(data)
  p = ncol(data)
  output = simulate_x(n, p, k, sigmasq, lambda_generator, return_lambda = TRUE)
  output[[1]] = output[[1]] + data
  output
}

S = 3
data_sets = lapply(1:S, function(s) add_simulated_latent_factors(data, k, sigmasq, simulate_lambda))
lambda_j = map(data_sets, function(x) pluck(x, 2))
data_sets = map(data_sets, function(x) pluck(x, 1))

map(data_sets, function(x) print(corrplot(cor(x))))

map(data_sets, function(x) apply(x, 2, var))

map(data_sets, head)

setwd("~/Homework/Fall 2018 Courses/Latent Factor/Data/")
j_s = k # For the multi_study package, this is the number of study-specific latent factors
save(data_sets, j_s, S, file="simulated_multi_study_datasets3.RData")

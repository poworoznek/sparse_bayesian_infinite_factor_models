# mcmc

mcmc <- function(Y, iter = 1000, burn = 500,
                 as = 1, bs = 0.3,              # gamma hyperparameters for residual precision
                 df = 3,                        # gamma hyperparameters for local shrinkage factor
                 ad1 = 2, bd1 = 1,              # gamma hyperparameters for delta_1
                 ad2 = 6, bd2 = 1,              # gamma hyperparameters delta_h, h >= 2
                 adf = 1, bdf = 1,              # gamma hyperparameters for ad1 and ad2 or df
                 b0 = 1, b1 = 0.0005,           # Adaptation hyperparameters
                 proportion = 1,                # Proportion of redundant elements before elimination of unnecessary column
                 epsilon = 1E-3                 # Cutoff for redundant element
){
  
  n = nrow(Y)
  p = ncol(Y)
  
  # Scale the data
  # var_y = apply(Y, 2, var)
  # scaleMat = sqrt((var_y) %*% t(var_y))
  # Y = scale(Y)
  # scaleMat = 1
  
  # consts = list("n" = n, "p" = p, "df" = df, "as" = as, "bs" = bs, "ad1" = ad1)
  # attach(consts)
  
  # local = true means that the fxn are evaluated from within this 'mcmc' function, so they can see the constant variables
  # Such as n, p, df, as, bs, ad1, bd1, etc... without passing them into the sampling fxns as arguments
  source("./BhattacharyaDunsonCode/sampling_functions.R", local = TRUE)

  # Set up the blank matrices to hold samples
  sigma_store = vector("list", length = iter)
  eta_store = vector("list", length = iter)
  lambda_store = vector("list", length = iter)
  delta_store = vector("list", length = iter)
  phi_store = vector("list", length = iter)
  tauh_store = vector("list", length = iter)
  omega_store = vector("list", length = iter)
  k_store = vector("list", length = iter)
  
  stores = list("lambda_store" = lambda_store, "sigma_store" = sigma_store, "eta_store" = eta_store,
                "phi_store" = phi_store, "delta_store" = delta_store, "tauh_store" = tauh_store,
                "omega_store" = omega_store, "k_store" = k_store)
  
  # Initialize the parameters
  k = round(log(p)*3, 0)
  sigma = 1/rgamma(p, as,bs)                    # residual noise after accounting for latent factors
  eta = matrix(rnorm(n*k), nrow = n, ncol = k)
  delta = c(rgamma(1, ad1, bd1), rgamma(k-1, ad2, bd2))
  tauh = cumprod(delta) # Global shrinkage param
  phi = matrix(rgamma(p*k, df/2, df/2), nrow = p, ncol = k) # Local shrinkage params
  Lambda = matrix(0, nrow = p, ncol = k)
  # If I simulated this data, then initialize the gibbs sampler at the true latent factors, give it help
  if(exists("output$lambda")){
    Lambda = cbind(output$lambda, matrix(0, nrow = p, ncol = k - ncol(output$lambda)))  
  }
  
  # Set up list of mcmc variables
  mc = list("sigma" = sigma, "eta" = eta, "delta" = delta, "tauh" = tauh, "phi" = phi,
            "Lambda" = Lambda, "k" = k)
  
  rm(list = c("sigma", "eta", "delta", "tauh", "phi", "Lambda", "k"))
  
  # Run the Gibbs Sampler
  for(i in -burn:iter){
    # for(i in 30:40){
    # Update Latent Factors (eta)
    mc = update_latent_factors(mc, Y)
    
    # Update Factor Loadings Matrix (Lambda)
    mc = update_factor_loadings(mc, Y)
    
    # Update phi
    mc = update_phi(mc)
    
    # Update delta & tau_h
    mc = update_delta(mc)
    
    # Update Sigma
    mc = update_sigma(mc, Y)
    
    # Make adaptations to 'k*', then number of latent factors
    mc = update_k(mc)
    
    # Save samples
    if(i >= 1){
      stores = save_samples(stores, mc, i)
    }
    
    if((i %% 1000) == 0) {
      print(i)
    }
    
  }
  
  # detach(consts)
  
  return(stores)
}

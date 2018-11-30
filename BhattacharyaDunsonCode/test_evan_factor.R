# Checking the latent factoring of evan

# First, read in his sampling function

# Then, run my code to generated simulated data, at the start of 'run_sim_example'

as = 1
bs = 0.3
df = 3
ad1 = 2
bd1 = 1
ad2 = 6
bd2 = 1
adf = 1
bdf = 1
b0 = 1
b1 = 0.0005
proportion = 1
epsilon = 1E-2

# NOTE THERE WAS ONE LINE IN FACT THAT I CHANGED - FIND THE COMMENT IN ALL CAPS
# IT'S ON UPDATING DELTA i.e. THETA
store_evan = fact(Y, b0 = b0, b1 = b1, as = as, bs = bs, df = df, ad1 = ad1, bd1 = bd1, ad2 = ad2, bd2 = bd2, adf = adf, bdf = bdf,
                  prop = proportion, epsilon = epsilon, nrun = iter, burn = burn,
                  thin = 1, kinit = NULL, output = c("covMean", "covSamples", "factSamples", "numFactors"))

round(store_evan$covMean, 1)
round(Omega, 1)

round(store_evan$covMean - Omega, 1) #Our omega estimates are very close

round(Omega - output$sigma, 1)

round(store_evan$covMean - output$sigma, 1)

table(store_evan$numFactors)
table(as.numeric(stores$k_store)) # But he consistently estimates a smaller k than I do... hmmm

mean(store_evan$numFactors)
mean(as.numeric(stores$k_store)) # Even after changing one small bug in his code, I think I am slightly consistently higher

get_lambda(store_evan)

# Testing the time - his code is about 20% faster, probably because mine is broken down into functions? I haven't looked at optimizing my code
# b/c so far sampler has been quite fast
iter = 1000
library(microbenchmark)
microbenchmark(mcmc(Y, iter, burn),
               fact(Y, b0, b1, as, bs, df, ad1, bd1, ad2, bd2, adf, bdf, prop = proportion, epsilon = epsilon, nrun = iter, burn = burn,
                    thin = 1, kinit = NULL, output = c("covMean", "covSamples", "factSamples", "numFactors")),
               times = 5)

# Settings --------------------------------------------
print("Settings...")

library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE)
library(gtools)
options(mc.cores = parallel::detectCores())

# Read data ---------------------------------------------
load("resBin1.RData")
ext = rstan::extract(res$stan_output, pars = c("beta", "gamma", "theta"), permuted = TRUE, inc_warmup = FALSE, include = TRUE)

# Compute probability of success -----------------------

niter = nrow(ext$beta)
prob_success = matrix(nrow = niter, ncol = res$datalist$N)
rr_c = vector(length = niter)
rr_t = vector(length = niter)
lor = vector(length = niter)


print("Computing probability of success...")
for (i in 1:niter) {
if(i %% 1000 == 0) {print(paste("Iteration ", i, sep = ""))}
   for (n in 1:res$datalist$N){
      prob_success[i,n] = inv.logit(ext$beta[i] + ext$gamma[i]*res$datalist$x[n] + ext$theta[i]*res$datalist$t[n])
   }
   rr_c[i] = mean(prob_success[i, res$datalist$t == 0])
   rr_t[i] = mean(prob_success[i, res$datalist$t == 1])
   lor[i] = log(rr_t[i]) - log(1 - rr_t[i]) - log(rr_c[i]) + log(1 - rr_c[i])
}

ps <- list(lor = lor)
save(ps, file  = "psBin1.RData")

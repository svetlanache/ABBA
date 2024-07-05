# Settings -------------------------------------
print("Settings...")
library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE)
library(gtools)
options(mc.cores = parallel::detectCores())

# Read data --------------------------------
load("resBin.RData")
ext = rstan::extract(res$stan_output, pars = c("beta", "gamma", "theta"), permuted = TRUE, inc_warmup = FALSE, include = TRUE)

# Compute probability of success -----------------------
print ("Computing probability of success...")
niter = nrow(ext$beta)
prob_success = matrix(nrow = niter, ncol = res$datalist$N)
rr_c = matrix(nrow = niter, ncol = res$datalist$nsub)
rr_t = matrix(nrow = niter, ncol = res$datalist$nsub)
lor = matrix(nrow = niter, ncol = res$datalist$nsub)

for (i in 1:niter) {
if(i %% 1000 == 0) {print(paste("Iteration ", i, sep = ""))}
   for (n in 1:res$datalist$N){
      prob_success[i,n] = inv.logit(ext$beta[i,res$datalist$sub[n]] + ext$gamma[i,res$datalist$sub[n]]*res$datalist$x[n] + ext$theta[i,res$datalist$sub[n]]*res$datalist$t[n])
   }
   for ( n in 1:res$datalist$nsub) {
     rr_c[i,n] = mean(prob_success[i,res$datalist$sub == n & res$datalist$t == 0])
     rr_t[i,n] = mean(prob_success[i,res$datalist$sub == n & res$datalist$t == 1])
     lor[i,n] = log(rr_t[i,n]) - log(1 - rr_t[i,n]) - log(rr_c[i,n]) + log(1 - rr_c[i,n])
   }
}

# save the results ------------------------------
ps <- list(lor = lor)
save(ps, file  = "psBin.RData")

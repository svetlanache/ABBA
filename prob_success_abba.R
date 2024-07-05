
# Settings --------------------------------------------
print("Settings...")
library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE)
library(mvtnorm)
options(mc.cores = parallel::detectCores())

# Read data ---------------------------------------------
load("resAbba.RData")
ext = rstan::extract(res$stan_output, pars = c("beta", "gamma", "theta", "L_Corr", "sd_y"), permuted = TRUE, inc_warmup = FALSE, include = TRUE)

# Compute probability of success -------------------------------

y1_lower = log(20)
y1_upper = Inf
y2_lower = 0
y2_upper = Inf
niter = nrow(ext$beta)
prob_success = matrix(nrow = niter, ncol = res$datalist$N)
rr_c = matrix(nrow = niter, ncol = res$datalist$nsub)
rr_t = matrix(nrow = niter, ncol = res$datalist$nsub)
lor = matrix(nrow = niter, ncol = res$datalist$nsub)

covmtx = array(dim = c(niter,2,2))
print ("Computing probability of success...")
for (i in 1:niter) {
  if(i %% 1000 == 0) {print(paste("Iteration ", i, sep = ""))}
  Sigma = c(ext$sd_y[i], 1)
  covmtx[i,,] = diag(Sigma) %*% (ext$L_Corr[i,,] %*% t(ext$L_Corr[i,,]))%*% diag(Sigma)
  for (n in 1:res$datalist$N){
    mu1 <- ext$beta[i,res$datalist$sub[n],1] + ext$gamma[i,res$datalist$sub[n],1]*res$datalist$x[n] + ext$theta[i,res$datalist$sub[n],1]*res$datalist$t[n]
    mu2 <- ext$beta[i,res$datalist$sub[n],2] + ext$gamma[i,res$datalist$sub[n],2]*res$datalist$x[n] + ext$theta[i,res$datalist$sub[n],2]*res$datalist$t[n]
    prob_success[i,n] <- pmvnorm(c(y1_lower, y2_lower), c(y1_upper, y2_upper), mean = c(mu1, mu2), sigma = covmtx[i,,])[1]
  }
  for ( n in 1:res$datalist$nsub) {
     rr_c[i,n] = mean(prob_success[i,res$datalist$sub == n & res$datalist$t == 0])
     rr_t[i,n] = mean(prob_success[i,res$datalist$sub == n & res$datalist$t == 1])
     lor[i,n] = log(rr_t[i,n]) - log(1 - rr_t[i,n]) - log(rr_c[i,n]) + log(1 - rr_c[i,n])
  }
}

# Save the results --------------------------------------------------------------

ps <- list(lor = lor)
save(ps, file  = "psAbba.RData")

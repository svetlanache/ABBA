

# Settings --------------------------------------------
print("Settings...")
library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE)
library(mvtnorm)
options(mc.cores = parallel::detectCores())

# Read data ---------------------------------------------
load("resAbba1.RData")
ext = rstan::extract(res$stan_output, pars = c("beta", "gamma", "theta", "L_Corr", "sd_y"), permuted = TRUE, inc_warmup = FALSE, include = TRUE)

# Compute probability of success -----------------------

y1_lower = log(20)
y1_upper = Inf
y2_lower = 0
y2_upper = Inf
niter = nrow(ext$beta)
prob_success = matrix(nrow = niter, ncol = res$datalist$N)
rr_c = vector(length = niter)
rr_t = vector(length = niter)
lor = vector(length = niter)

covmtx = array(dim = c(niter,2,2))
print ("Computing probability of success...")
for (i in 1:niter) {
  if(i %% 1000 == 0) {print(paste("Iteration ", i, sep = ""))}
  Sigma = c(ext$sd_y[i], 1)
  covmtx[i,,] = diag(Sigma) %*% (ext$L_Corr[i,,] %*% t(ext$L_Corr[i,,]))%*% diag(Sigma)
  for (n in 1:res$datalist$N){
    mu1 <- ext$beta[i,1] + ext$gamma[i,1]*res$datalist$x[n] + ext$theta[i,1]*res$datalist$t[n]
    mu2 <- ext$beta[i,2] + ext$gamma[i,2]*res$datalist$x[n] + ext$theta[i,2]*res$datalist$t[n]
    prob_success[i,n] <- pmvnorm(c(y1_lower, y2_lower), c(y1_upper, y2_upper), mean = c(mu1, mu2), sigma = covmtx[i,,])[1]
  }
  rr_c[i] = mean(prob_success[i, res$datalist$t == 0])
  rr_t[i] = mean(prob_success[i, res$datalist$t == 1])
  lor[i] = log(rr_t[i]) - log(1 - rr_t[i]) - log(rr_c[i]) + log(1 - rr_c[i])
}

## Save the results

ps <- list(lor = lor)
save(ps, file  = "psAbba1.RData")

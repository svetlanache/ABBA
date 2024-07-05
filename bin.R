
# Implementation of the BIN model -----------------------
# --------------------------------------------------------

# Settings -----------------------------------------------

print("Settings...")

library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE)
library(gtools)
options(mc.cores = parallel::detectCores())

# Load data --------------------------------------------
load("datalist.RData")

# Model ------------------------------------------------

model_code <- "
functions {
  real to_real(int x) { return x; }
}


data {
  int<lower=1> N;
  int<lower=0, upper=1>  b[N];
  vector[N]  y;
  int nsub;
  int  sub[N];
  int<lower=0,upper=1> t[N];
  vector[N] x;
}

transformed data {
  int<lower=0, upper=1> outcome[N];

  for (i in 1:N){
    outcome[i] = 0;
    if ((y[i] > log(20)) && (b[i] == 1)) {outcome[i] = 1;}
  }
}

parameters {
  real beta[nsub];
  real gamma[nsub];
  real theta[nsub];
  real beta_mean;
  real gamma_mean;
  real theta_mean;
  real<lower=0.3> sd_beta;
  real<lower=0.3> sd_gamma;
  real<lower=0.3> sd_theta;
}

transformed parameters {
  real outcome_mean[N];

  for (n in 1:N) {

    outcome_mean[n] = beta[sub[n]] + x[n]*gamma[sub[n]] + t[n]*theta[sub[n]];
  }
}

model {

  beta_mean ~ normal(0,5);
  gamma_mean ~ normal(0,5);
  theta_mean ~ normal(0,5);
  sd_beta ~ exponential(2);
  sd_gamma ~ exponential(2);
  sd_theta ~ exponential(2);

  to_vector(beta) ~ normal(beta_mean,sd_beta);
  to_vector(gamma) ~ normal(gamma_mean,sd_gamma);
  to_vector(theta) ~ normal(theta_mean,sd_theta);

  for (n in 1:N) 
  {
     outcome[n] ~ bernoulli_logit(outcome_mean[n]);
  }  
}
"

# Run the model ------------------------------------------

print("Run model...")

iter = 15000
warmup = 5000
chains = 2
stan_output <- stan(model_code = model_code, data = datalist, iter = iter, warmup = warmup, chains = chains, cores = 2, 
       par = c("beta", "gamma", "theta", "beta_mean", "gamma_mean", "theta_mean", "sd_beta", "sd_gamma", "sd_theta"), 
       include = TRUE)


# Save the results --------------------------------------------------------------------

res <- list(stan_output = stan_output, datalist = datalist)
save(res, file = "resBin.RData")


# Implementation of the ABBA model ------------
# ----------------------------------------------

# Settings-------------------------------------

print("Settings...")

library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE)
library(bayestestR)
library(mvtnorm)
options(mc.cores = parallel::detectCores())

# Data load -----------------------------------

load("datalist.RData")

# Model ----------------------------------------

model_code <- "
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
  int<lower=0> N_pos;
  int<lower=0> N_neg;
  int<lower=1, upper=(N+1)> ipos;
  int<lower=1, upper=(N+1)> ineg;
  int<lower=1, upper=N>  n_pos[sum(b)];
  int<lower=1, upper=N>  n_neg[N - size(n_pos)];
  int<lower=0> threshold;

  N_pos = size(n_pos);
  N_neg = size(n_neg);

  threshold = 20;
  ipos = 1;
  ineg = 1;
  for (n in 1:N) {
    if (b[n] == 1) {
       n_pos[ipos] = n;
       ipos += 1;
    } else {
       n_neg[ineg] = n;
       ineg += 1;
    }
  }
}

parameters {
  vector<lower=0>[N_pos] z_pos;
  vector<upper=0>[N_neg] z_neg;
  matrix[nsub, 2] beta;
  matrix[nsub, 2] gamma;
  matrix[nsub, 2] theta;
  vector[2] mean_beta;
  vector[2] mean_gamma;
  vector[2] mean_theta;
  cholesky_factor_corr[2] L_Corr;
  real<lower=0.00001> sd_y;
  real<lower=0.1> sd_beta1;
  real<lower=0.1> sd_beta2;
  real<lower=0.1> sd_gamma1;
  real<lower=0.1> sd_gamma2;
  real<lower=0.1> sd_theta1;
  real<lower=0.1> sd_theta2;
}

transformed parameters {

  matrix[N, 2] outcome_mean;
  matrix[N, 2] outcome;
  vector<lower=0>[2] outcome_sd;
  matrix[2, 2] L_Sigma;

  outcome_sd[1] = sd_y;
  outcome_sd[2] = 1;

  L_Sigma = diag_pre_multiply(outcome_sd, L_Corr);

  outcome[,1] = y;
  for (n in 1:N_pos) {
    outcome[n_pos[n], 2] = z_pos[n];
  }

  for (n in 1:N_neg) {
    outcome[n_neg[n], 2] = z_neg[n];
  }

  for (n in 1:N) {
    outcome_mean[n,1] = beta[sub[n],1] + x[n]*gamma[sub[n],1] + t[n]*theta[sub[n],1];
    outcome_mean[n,2] = beta[sub[n],2] + x[n]*gamma[sub[n],2] + t[n]*theta[sub[n],2];
  }
}

model {

  sd_y ~ inv_gamma(0.5, 0.005);
  L_Corr ~ lkj_corr_cholesky(5);

  mean_beta ~ normal(0,5);
  mean_gamma ~ normal(0,5);
  mean_theta ~ normal(0,5);

  sd_beta1 ~ exponential(2);
  sd_beta2 ~ exponential(2);
  sd_gamma1 ~ exponential(2);
  sd_gamma2 ~ exponential(2);
  sd_theta1 ~ exponential(2);
  sd_theta2 ~ exponential(2);

  to_vector(beta[,1]) ~ normal(mean_beta[1],sd_beta1);
  to_vector(beta[,2]) ~ normal(mean_beta[2],sd_beta2);

  to_vector(gamma[,1]) ~ normal(mean_gamma[1],sd_gamma1);
  to_vector(gamma[,2]) ~ normal(mean_gamma[2],sd_gamma2);

  to_vector(theta[,1]) ~ normal(mean_theta[1],sd_theta1);
  to_vector(theta[,2]) ~ normal(mean_theta[2],sd_theta2);

  for (n in 1:N) {
     outcome[n,] ~ multi_normal_cholesky(outcome_mean[n,], L_Sigma);
  }

}
"
# Initialise ------------------------------------

z_pos_init = as.list(rep(1, sum(datalist$b == 1)))
z_neg_init = as.list(rep(-1, sum(datalist$b == 0)))
init_fun <- function() {
       list(z_pos = z_pos_init, z_neg = z_neg_init,
            beta = matrix(nrow = 3, ncol = 2, data = rep(1,6)),
            gamma = matrix(nrow = 3, ncol = 2, data = rep(1,6)),
            theta = matrix(nrow = 3, ncol = 2, data = rep(1,6)),
            mean_beta = c(0,1), mean_gamma = c(0,1), mean_theta = c(0,1),
            L_Corr = t(chol(matrix(nrow = 2, ncol = 2, data = c(1, 0.3, 0.3, 1), byrow = T))),
            sd_y = 0.5, sd_beta1 = 1.1, sd_beta2 = 1.1, sd_gamma1 = 1.1, sd_gamma2 = 1.1, sd_theta1 = 1.1, sd_theta1 = 1.1)}

# Run the model -----------------------------------

print("Run model...")
iter = 15000
warmup = 5000
chains = 2
stan_output <- stan(model_code = model_code, data = datalist, iter = iter, warmup = warmup, chains = chains, cores = 2, init = init_fun, 
       par = c("mean_beta", "mean_gamma", "mean_theta","sd_beta1", "sd_beta2", "sd_gamma1", "sd_gamma2", 
               "sd_theta1", "sd_theta2", "beta", "gamma", "theta", "L_Corr", "sd_y"),  include = TRUE) 
  

# Save the results ------------------------------------------------

res <- list(stan_output = stan_output, datalist = datalist)
save(res, file  = "resAbba.RData")

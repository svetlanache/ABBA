
# Implementation of the ABBA model for stratified analysis -----------------------
# --------------------------------------------------------------------------------

# Settings ----------------------------------------------------------------------

print("Settings...")
library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load data ---------------------------------------------------------------------

load("datalist1.RData")

# Model -----------------------------------------------------------------------

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
  vector[2] beta;
  vector[2] gamma;
  vector[2] theta;
  cholesky_factor_corr[2] L_Corr;
  real<lower=0.00001> sd_y;
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
    outcome_mean[n,1] = beta[1] + x[n]*gamma[1] + t[n]*theta[1];
    outcome_mean[n,2] = beta[2] + x[n]*gamma[2] + t[n]*theta[2];
  }
}

model {

  sd_y ~ inv_gamma(0.5,0.005);
  L_Corr ~ lkj_corr_cholesky(5);

  beta[1] ~ normal(0,5);
  beta[2] ~ normal(0,5);

  gamma[1] ~ normal(0,5);
  gamma[2] ~ normal(0,5);

  theta[1] ~ normal(0,10);
  theta[2] ~ normal(0,10);

  for (n in 1:N) {
     outcome[n,] ~ multi_normal_cholesky(outcome_mean[n,], L_Sigma);
  }

}

generated quantities {

     matrix[N, 2] outcome_post;
     vector<lower=0, upper=1>[N] responder;
     vector<lower=0, upper=N>[nsub] resp_control;
     vector<lower=0, upper=N>[nsub] resp_treat;
     vector<lower=1, upper=N>[nsub] count_control;
     vector<lower=1, upper=N>[nsub] count_treat;
     vector<lower=0, upper=1>[nsub] rr_control;
     vector<lower=0, upper=1>[nsub] rr_treat;
     vector[nsub] logoddsratio;
     real nudge;
     nudge=0.01;
     
     for (n in 1:nsub) {
       resp_control[n] = 0;
       resp_treat[n] = 0;
       count_control[n] = 0;
       count_treat[n] = 0;
     }

     for (n in 1:N)
     {
       responder[n] = 0;
       outcome_post[n,1:2] = to_row_vector(multi_normal_cholesky_rng(outcome_mean[n,], L_Sigma));
       if ((outcome_post[n,1] > log(threshold)) && (outcome_post[n,2] > 0)) {responder[n] = 1;}
       if (t[n] == 0) {
          resp_control = resp_control + responder[n];
          count_control = count_control + 1;
        } else if (t[n] == 1) {
          resp_treat = resp_treat + responder[n];
          count_treat  = count_treat + 1;
        }
     }
      
     rr_control = resp_control ./ count_control; 
     rr_treat = resp_treat ./ count_treat;
     logoddsratio = log(rr_treat + nudge) + log(1 - rr_control + nudge) - log(1 - rr_treat + nudge) - log(rr_control + nudge);
}
"

# Initialise ---------------------------------------------------------------------------------

z_pos_init = as.list(rep(1, sum(datalist$b == 1)))
z_neg_init = as.list(rep(-1, sum(datalist$b == 0)))
init_fun <- function() {
       list(z_pos = z_pos_init, z_neg = z_neg_init,
            beta = c(1,1),
            gamma = c(1,1),
            theta = c(1,1),
            L_Corr = t(chol(matrix(nrow = 2, ncol = 2, data = c(1, 0.3, 0.3, 1), byrow = T))),
            sd_y = 0.5)
}


# Run model -----------------------------------------------------------------------------------

print ("Running model...")

stan_output <- stan(model_code = model_code, data = datalist, iter = 15000, warmup = 5000, chains = 2, cores = 2, init = init_fun, 
       par = c("beta", "gamma", "theta", "L_Corr", "sd_y", "rr_control", "rr_treat", "logoddsratio"), include = TRUE) 
       
res <- list(stan_output = stan_output, datalist = datalist)

save(res, file  = "resAbba1.RData")

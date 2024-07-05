
# Implementation of the BIN model for stratified analysis --------------
# ----------------------------------------------------------------------

# Settings -------------------------------------------------------------

print("Settings...")
library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load ata--------------------------------------------------------------

load("datalist1.RData")

# Model ----------------------------------------------------------------

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
  real beta;
  real gamma;
  real theta;
}

transformed parameters {
  real outcome_mean[N];

  for (n in 1:N) {

    outcome_mean[n] = beta + x[n]*gamma + t[n]*theta;
  }
}

model {

  beta ~ normal(0,5);
  gamma ~ normal(0,5);
  theta ~ normal(0,10);

  for (n in 1:N) 
  {
     outcome[n] ~ bernoulli_logit(outcome_mean[n]);
  }  
}

generated quantities {

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

     for (n in 1:(nsub)) {
       resp_control[n] = 0;
       resp_treat[n] = 0;
       count_control[n] = 0;
       count_treat[n] = 0;
     }

     for (n in 1:N)
     {
       responder[n] = bernoulli_logit_rng(outcome_mean[n]);
       if (t[n] == 0) {
          resp_control = resp_control + responder[n];
          count_control = count_control + 1;
        } else if (t[n] == 1) {
          resp_treat = resp_treat + responder[n];
          count_treat = count_treat + 1;
        }
     }

     rr_control = resp_control ./ count_control;
     rr_treat = resp_treat ./ count_treat;
     logoddsratio = log(rr_treat + nudge) + log(1 - rr_control + nudge) - log(1 - rr_treat + nudge) - log(rr_control + nudge);

}
"
stan_output <- stan(model_code = model_code, data = datalist, iter = 15000, warmup = 5000, chains = 2, cores = 2, 
       par = c("beta", "gamma", "theta", "rr_control", "rr_treat", "logoddsratio"), 
       include = TRUE)

res <- list(stan_output = stan_output, datalist = datalist)
save(res, file  = "resBin1.RData")

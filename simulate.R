
# Description ----------------------------
# Simulated data set for a basket trial with 3 subtrials, 50 patients per subtrial, equally randomised to treatment and control.
# The data set includes:
# n - number of patients in each subtrial
# N - overall number of patients
# sub - an indicator of a subtrial (1, 2, or 3)
# x - covariate randomly distributed 0-100, on the log scale
# y, b - countinous and binary components of the outcome.
# (b is set as 0 or 1 depending on the latent variable. Latent variable simulated from a multivariate normal distribution.
# The dataset is saved as datalist.RData R object.


#Setting --------------------------------
library(MASS)
set.seed(123)

# Load parameters ------------------------
load("param_scenario1.RData")

nsub <- 3
n <- 50
N <- n*nsub
sub <- rep(1:nsub, each = n)
t <- rep (0:1, each = n/2, nsub)
x <- log(runif(N, 1, 100)) # log baseline score
thr = 20 # threshold for the log score

# Create data ----------------------------

print ("Simulating data...")

beta1 <- param$beta1
gamma1 <- param$gamma1
theta1 <- param$theta1

beta2 <- param$beta2
gamma2 <- param$gamma2
theta2 <- param$theta2

y_mtx = matrix(nrow = N, ncol = 2)
y = vector(length = N)
b = vector(length = N)
mean_y1 <- vector(length = N)
mean_y2 <- vector(length = N)
resp <- vector(length = N)
sd1 = 0.5
sd2 = 1
rho = 0.3
Sigma_y = matrix(data = c(sd1^2, sd1*sd2*rho, sd1*sd2*rho, sd2^2), nrow = 2, ncol = 2, byrow = T)

for (i in 1:N) {
  mean_y1[i] <- beta1[sub[i]] + gamma1[sub[i]]*x[i] + theta1[sub[i]]*t[i]
  mean_y2[i] <- beta2[sub[i]] + gamma2[sub[i]]*x[i] + theta2[sub[i]]*t[i]
  y_mtx[i,] = mvtnorm::rmvnorm(n = 1, mean = c(mean_y1[i], mean_y2[i]), sigma = Sigma_y)
  if ((y_mtx[i,1] > log(thr)) & (y_mtx[i,2] > 0)) {resp[i] = 1}
}

y = y_mtx[,1]
b = ifelse(y_mtx[,2] > 0, 1, 0)

rr_c = vector(length = nsub)
rr_t = vector(length = nsub) 
for (i in 1:nsub) {
  rr_c[i] <- mean(resp[sub == i & t == 0])
  rr_t[i] <- mean(resp[sub == i & t == 1])
}
  
datalist = list(N = N, x = x, y = y, b = b, t = t, sub = sub, nsub = nsub,
                resp = resp, rr_c = rr_c, rr_t = rr_t, param = param)
save(datalist, file = "datalist.RData")

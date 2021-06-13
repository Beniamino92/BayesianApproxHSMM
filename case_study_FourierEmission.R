# (your current directory must contain both /stan and /include)
setwd("/Users/beniamino/Desktop/reviewed_BayesianApproxHSMM/")

library("rstan")
library("bayesplot")
library("bridgesampling")
library("matrixStats")
library("lubridate")


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


source("include/util_hmm.R")
source("include/util_negativeBinomial.R")
source("include/util_poisson.R")
#source("include/util_oscillatory.R")


RA.data <- read.csv(paste0(getwd(), "/data/S16.csv"))

# --- NAs (mean adjecent point)
NAs <- which(is.na(RA.data$activity)) # only 3 NA's
for (i in 1:length(NAs)) {
  ind_lower <- NAs[i] - 1
  ind_upper <- NAs[i] + 1
  RA.data$activity[NAs[i]] <- mean(RA.data$activity[c(ind_lower, ind_upper)])
}

# -- transformed data (sqrt)
obs <- sqrt(RA.data$activity) 
N <- length(obs)

## Comparable Priors - Set Up

K <- 3 # n of states
alpha_0 <- matrix(NA, K, K) # hyperparms dirichlet
a_0 <- numeric(K) # hyperparms gamma
b_0 <- numeric(K) #       "


state <- 1
mean.target <- 7.5*12
var.target <- (0.5*12)^2

dir.parms <- optim(par = rep(5, 2),
                   fn = function(par) {
                     geom.dwell.error(exp(par[1]), exp(par[2]), 
                                      mean.target, var.target)
                   })

# we assume equal probability of going to state 2 (MA) or 3 (HA)
est <- exp(dir.parms$par)
alpha_0[state, ] <- c(est[1], est[2]/2, est[2]/2)

gamma.parms <- optim(par = rep(5, 2), 
                     fn = function(par) {
                       nb.dwell.error(exp(par[1]), exp(par[2]), 
                                      mean.target, var.target)
                     })
est <- exp(gamma.parms$par)
a_0[state] <- est[1]; 
b_0[state] <- est[2];


geom.dwell.info(alpha_0[state, state], sum(alpha_0[state, -state]))
nb.dwell.info(a_0[state], b_0[state])

# mean: 1x12 = 12; var: (1.5x12)^2 = 324. 
# prior information:  the other 16-17 hours are split between the other two states, i.e.
# moderately active and highly active states. In general, we have much less information 
# on how long peopel spend in  these states. We fix mean 2 hours and a large variance.
state <- 2
mean.target <- 2*12
var.target <- (1.5*12)^2

dir.parms <- optim(par = c(5, 2),
                   fn = function(par) {
                     geom.dwell.error(exp(par[1]), exp(par[2]), 
                                      mean.target, var.target)
                   })
# MA and HA alternates between each other more frequently than with IA, 
# Hence we may assume transitions from MA->HA is 4 times more likely 
# than going from MA->IA
est <- exp(dir.parms$par)
alpha_0[state, ] <- c(est[2]/5, est[1], est[2]*(4/5))

gamma.parms <- optim(par = c(5, 5), 
                     fn = function(par) {
                       nb.dwell.error(exp(par[1]), exp(par[2]), 
                                      mean.target, var.target)
                     })
est <- exp(gamma.parms$par)
a_0[state] <- est[1]; 
b_0[state] <- est[2];


geom.dwell.info(alpha_0[state, state], sum(alpha_0[state, -state]))
nb.dwell.info(a_0[state], b_0[state])


# mean: 1x12 = 12; var: (1.5x12)^2 = 324. 
# prior information:  the other 16-17 hours are split between the other two states, i.e.
# MA and HA states. In general, we have much less information 
# on how long peopel spend in these states. We fix mean 2 hours and a large variance.

state <- 3
mean.target <- 2*12
var.target <- (1.5*12)^2

dir.parms <- optim(par = rep(5, 2),
                   fn = function(par) {
                     geom.dwell.error(exp(par[1]), exp(par[2]), 
                                      mean.target, var.target)
                   })

# MA and HA alternates between each other more frequently than with IA, 
# Hence we may assume transitions from HA->MA is 4 times more likely 
# than going from HA->IA
est <- exp(dir.parms$par)
alpha_0[state, ] <- c(est[2]/5, est[2]*(4/5),est[1])

gamma.parms <- optim(par = rep(5, 2), 
                     fn = function(par) {
                       nb.dwell.error(exp(par[1]), exp(par[2]), 
                                      mean.target, var.target)
                     })
est <- exp(gamma.parms$par)
a_0[state] <- est[1]; 
b_0[state] <- est[2];

geom.dwell.info(alpha_0[state, state], sum(alpha_0[state, -state]))
nb.dwell.info(a_0[state], b_0[state])



## Bayesian HSMM Approx - Negative Binomial Dwell

K <- 3 # n of states
m <- c(150, 10, 10) # dwell threshold

# hyperparms
sigma_0 <- 2
# mu_0 <- rep(mean(obs), K)
mu_0 <- c(mean(obs), 0, 0)
a_0; b_0
v_0 <- matrix(c(alpha_0[1, 2:3], alpha_0[2, c(1, 3)], alpha_0[3, 2:3]),
              nrow = K, ncol = K -1, byrow = TRUE)

# test for sparsity 
path.stan <- "stan/bayesHSMMapprox_HarmonicEmis_NegBinomDur_sparse.stan"

HSMM.data <- list(N = length(obs), K = K, y = obs, omega = 0.003453863, m = m, mu_0 = mu_0,
                  sigma_0 = sigma_0, a_0_lambda = a_0, b_0_lambda = b_0,
                  a_0_phi = 2, b_0_phi = 2, alpha_0 = v_0)
HSMM.stan <- stan(file = path.stan, data = HSMM.data,
                  warmup = 1000, iter=(1+5)*1000, chains=1, cores=1)
print(HSMM.stan, probs = c(0.05, 0.95))
HSMM.bridge <- bridge_sampler(HSMM.stan)
HSMM.lml <- HSMM.bridge$logml
save(HSMM.stan, HSMM.lml, file = paste0(getwd(), "harmonic_HSMM_negbinom.RData"))



## Bayesian HMM - Geometric Dwell

K <- 3 # n of states

# hyperparams
# mu_0 <- rep(mean(obs), K)
mu_0 <- c(mean(obs), 0, 0)
sigma_0 <- 2
alpha_0



HMM.data <- list(N = length(obs), K = K, y = obs, omega = 0.003453863, 
                 mu_0 = mu_0, sigma_0 = sigma_0, alpha_0 = alpha_0)
HMM.stan <- stan(file = "stan/bayesHMM_HarmonicEmis.stan", data = HMM.data,
                 warmup = 1000, iter=(1+5)*1000, chains=1, cores=1)
HMM.bridge <- bridge_sampler(HMM.stan)
HMM.bridge$logml
print(HMM.stan, probs = c(0.05, 0.95))
HMM.lml <- HMM.bridge$logml
save(HMM.stan, HMM.lml, file = paste0(getwd(), "harmonic_HMM.RData"))


## HSMM -- Poisson Dwell: 


# hyperparms
K <- 3 # n of states
m <- c(150, 10, 10) # dwell threshold

sigma_0 <- 2
mu_0 <- c(mean(obs), 0, 0)
a_0; b_0
v_0 <- matrix(c(alpha_0[1, 2:3], alpha_0[2, c(1, 3)], alpha_0[3, 2:3]),
              nrow = K, ncol = K -1, byrow = TRUE)


path.stan <- "stan/bayesHSMMapprox_HarmonicEmis_PoissDur_sparse.stan"
HSMM.data <- list(N = length(obs), K = K, y = obs, omega = 0.003453863, m = m, mu_0 = mu_0,
                  sigma_0 = sigma_0, a_0 = rep(0.01, K), b_0 = rep(0.01, K),
                  alpha_0 = matrix(1, nrow = K, ncol = K-1))
HSMM.stan <- stan(file = path.stan, data = HSMM.data,
                  warmup = 1000, iter=(1+5)*1000, chains=1, cores=1)
HSMM.bridge <- bridge_sampler(HSMM.stan)
HSMM.lml <- HSMM.bridge$logml
save(HSMM.stan, HSMM.lml, file = paste0(getwd(), "harmonic_HSMM_pois.RData"))


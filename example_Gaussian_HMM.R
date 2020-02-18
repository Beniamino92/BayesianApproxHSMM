library("rstan")
library("bayesplot")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
setwd("/Users/beniamino/Desktop/My_PostDoc/HSMM_Approx/")



# Generate a realization of a Gaussian–HMM (fixed state-specific sd)
gauss.HMM.generate_sample <- 
  function(n, m, mu, gamma, sd = 1, delta = NULL)
  {
    if(is.null(delta))  delta <- solve(t(diag(m)-gamma + 1), rep(1,m))
    mvect <- 1:m
    state <- numeric(n)
    state[1] <- sample(mvect, 1, prob = delta)
    for (i in 2:n)
      state[i] <- sample(mvect, 1, prob = gamma[state[i-1], ])
    x <- rnorm(n, mean = mu[state], sd = sd)
    list(x = x, state = state)
  }



set.seed(1)

# True generating parameters
m = 2
n <- 200 
gamma <- matrix(c(0.9, 0.1, 0.1, 0.9), m, m, byrow = T) 
#lambda <-  c(1, 5)
mu <- c(1, 6)
delta <- c(1.0, 0.0)

# Simulated data and sequence
simul <- gauss.HMM.generate_sample(n, m, mu, gamma)
data <- simul$x
state <- simul$state

# Plot data and sequence
par(mfrow = c(2, 1))
plot(data, type= "p", pch = 20, cex = 0.6, col = state)
plot(1:n, state, lwd = 2, type = "l")

# Stan model
stan_data <- list(N = n, 
                  K = m, 
                  y = data)

hmm_fit <- stan("Gaussian_HMM.stan", data = stan_data, iter = 1e3, chains = 4)
print(hmm_fit, pars = "z_star", include = FALSE, probs = c(0.05, 0.95))

# Posterior sample
sims <- extract(hmm_fit)
z_est <- sims$z_star[1, ]

# Plot: traces
color_scheme_set("viridisA")
mcmc_trace(as.array(hmm_fit), regex_pars = "^gamma\\[|^mu\\[", facet_args = list(nrow = 2))

# Plot: estimated state sequence
par(mfrow = c(2, 1))
plot(1:n, z_est, lwd = 2, type = "l")
plot(data, type = "p", cex = 0.7, col = state, pch = 20)


# Diagnostic: posterior predictive check.

samples <- as.matrix(hmm_fit)
mu <- samples[,grep("^mu",colnames(samples))]
z_star <- samples[,grep("^z_star",colnames(samples))]

y_hat <- list()
for (i in 1:nrow(samples)) {
  psi_seq <- sapply(z_star[i, ], function(x) {mu[i, x]})
  y_hat[[i]] <- rnorm(n, psi_seq, 1)
}

indxs <- sample(length(y_hat), 100, replace = FALSE)

# Plot: observed vs predicted output
par(mfrow = c(1, 1))
plot(data, type = "n", ylim = c(-4, 12),
     main = "Observed vs Predicted Output", 
     ylab = "Observation Value", 
     xlab = "Time")
for (i in indxs) {
  lines(y_hat[[i]], col = "#ff668890")
}
lines(data, lwd = 2)
legend("bottomright", c("Observed","Predicted"), 
       col = c("#000000","#ff668890"), 
       lty = c(1,1), lwd = c(2,1), cex = 0.8)


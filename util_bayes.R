#  simplex to tpm 
a_to_mat <- function(a) 
{
  if (length(a) == 2) {
    K <- 2
  } else {
    K <- nrow(a)
  }
  gamma_out <- matrix(0, K, K)
  gamma_out[!diag(K)] <- as.vector(t(a))
  t(gamma_out)
}

# get gamma sims from a_sims
get.gamma_sims <- function(a) 
{
  n_iter = nrow(a)
  get.gamma <- function(t, a) {a_to_mat(a[t, , ])}
  simplify2array(lapply(1:n_iter, get.gamma, a = sims$gamma))
}


# trace predictive likelihoood (log-likelihood) 
HSMM.trace.mllk <- function(sims, plt = FALSE) 
{
  # posterior samples
  lambda_sims <- sims$lambda
  mu_sims <- sims$mu
  sigma_sims <- sims$sigm
  gamma_sims <- get.gamma_sims(sims$a)
  
  n_iter <- dim(lambda_sims)[1]
  mllk_trace <- numeric(n_iter)
  
  for (t in 1:n_iter) {
    if((t %% (n_iter/5)) == 0) {
      cat(" ...", as.integer((t/n_iter)*100), "% \n")
    }
    lambda.hat <- lambda_sims[t, ]
    mu.hat <- mu_sims[t, ]
    sigma.hat <- sigma_sims[t, ]
    gamma.hat <- gamma_sims[, , t]
      
    mllk_trace[t] <- HSMM.mllk(obs, m, lambda.hat, mu.hat, 
                               sigma.hat, gamma.hat)
  }
  
  if (plt == TRUE) {
    plot(mllk_trace, type = "l", xlab = "MCMC Iteration", ylab = "mllk")
  }
  mllk_trace
}

# predictive accuracy measures: AIC, BIC, DIC, WAIC (see BDA Ch.7)
HSMM.performance <- function(sims, obs, m) 
{
  lambda_sims <- sims$lambda
  mu_sims <- sims$mu
  sigma_sims <- sims$sigma
  gamma_sims <- get.gamma_sims(sims$gamma)
  
  lambda.hat <- colMeans(lambda_sims)
  mu.hat <- colMeans(mu_sims)
  sigma.hat <- colMeans(sigma_sims)
  gamma.hat <- apply(gamma_sims, c(1, 2), mean)
  
  n_iter <- nrow(lambda_sims)
  K <- ncol(lambda_sims)
  N <- length(obs)
  np <- if(K == 2) 3*K else K^2 - 2*K + 3*K
  
  temp <- matrix(NA, N, n_iter)
  llk.trace <- numeric(n_iter)
  cat(" ... computing predictive performance ... \n")
  for (t in 1:n_iter) {
    if((t %% (n_iter/10)) == 0) {
      cat(" ...", as.integer((t/n_iter)*100), "% \n")
    }
    lambda <- lambda_sims[t, ]
    mu <- mu_sims[t, ]
    sigma <- sigma_sims[t, ]
    gamma <- gamma_sims[, , t]
    z <- HSMM.viterbi(obs, m, lambda, mu, sigma, gamma, draw = TRUE)
    llk.trace[t] <- -HSMM.mllk(obs, m, lambda, mu, sigma, gamma)
    temp[, t] <- dnorm(obs, mu[z], sigma[z])
  }
  
  lppd <- sum(log(rowMeans(temp)))
  llk <- -HSMM.mllk(obs, m, lambda.hat, mu.hat, sigma.hat, gamma.hat)
  p.DIC <- 2*(llk - mean(llk.trace))
  p.WAIC <- 2*sum(log(rowMeans(temp)) - rowMeans(log(temp)))
  
  AIC <- -2*(llk - np)
  BIC <- -2*llk+np*log(N)
  DIC <- -2*llk + 2*p.DIC 
  WAIC <- lppd - p.WAIC
  
  list(mmlk = -llk, AIC = AIC, BIC = BIC, DIC = DIC, WAIC = -WAIC)
}

# AIC: uses dplug-in Bayes estimate
# DIC: uses plug-in Bayes estimate
# WAIC: uses all posterior sample.
# BIC: (not recommended by Gelman: it is not inteded to predict out-of-sample model
#      performance but rather for other purposes, BDA p.175).


# posterior predictive check: draw from posterior predictive + plot
HSMM.predictive.plot <- function(sims, obs, m, ndraw = 100)
{
  y_hat <- list()
  niter <- nrow(sims$mu)
  
  lambda_sims <- sims$lambda
  mu_sims <- sims$mu
  sigma_sims <- sims$sigma
  gamma_sims <- get.gamma_sims(sims$gamma)
  
  cat(" ... sampling from posterior predictive... \n")
  for (t in 1:nrow(sims$mu))  {
    if((t %% (niter/10)) == 0) {
      cat(" ...", as.integer((t/niter)*100), "% \n")
    }
    lambda_temp <- lambda_sims[t, ]
    mu_temp <- mu_sims[t, ]
    sigma_temp <- sigma_sims[t, ]
    gamma_temp <- gamma_sims[, , t]
    
    z_star <- HSMM.viterbi(obs, m, lambda_temp, mu_temp, 
                           sigma_temp, gamma_temp, draw = TRUE)
    y_hat[[t]] <- rnorm(N, mu_temp[z_star], sigma_temp[z_star])
  }
  
  # plot: posterior predictive check
  par(mfrow = c(1, 1))
  lambda.hat <- colMeans(lambda_sims)
  mu.hat <- colMeans(mu_sims)
  sigma.hat <- colMeans(sigma_sims)
  gamma.hat <- apply(gamma_sims, c(1, 2), mean)
  
  z_star <- HSMM.viterbi(obs, m, lambda.hat, mu.hat, 
                         sigma.hat, gamma.hat, draw = TRUE)
  indxs <- sample(length(y_hat), ndraw, replace = FALSE)
  plot(obs, type = "n", ylab = "Observations", 
       xlab = "Time", ylim = c(min(obs) - 3, max(obs) + 3))
  for (i in indxs) {
    lines(y_hat[[i]], col = "gray85")
  }
  points(1:N, obs, pch = 20, cex = 0.6, col = z_star)
}


# hist: transtion matrix posterior histograms
HSMM.transitions.hist <- function(gamma_sims)
{
  K <- nrow(gamma_sims)
  if (K == 2) stop("tpm is deterministic")
  par(mfrow = c(K, K))
  for (i in 1:K) {
    for (j in 1:K) {
      if (i == j) {
        plot(0,type='n',axes=FALSE,ann=FALSE)
        box(col = 'black')
      } else {
        hist(gamma_sims[i, j, ], col = "pink", breaks = 20,
             main = bquote(gamma[.(i)*.(j)]), xlim = c(0, 1), freq = TRUE, 
             ylab = "", xlab = "")
        box(col = 'black')
      }
    }
  }
  par(mfrow = c(1, 1))
}



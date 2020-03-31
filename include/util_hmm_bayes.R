# posterior predictive check: draw from posterior predictive + plot
HMM.predictive.plot <- function(sims, obs, ndraw = 100)
{
  y_hat <- list()
  niter <- nrow(sims$mu)
  
  mu_sims <- sims$mu
  sigma_sims <- sims$sigma
  gamma_sims <- sims$gamma
  
  cat(" ... sampling from posterior predictive... \n")
  for (t in 1:nrow(sims$mu))  {
    if((t %% (niter/10)) == 0) {
      cat(" ...", as.integer((t/niter)*100), "% \n")
    }
    
    mu_temp <- mu_sims[t, ]
    sigma_temp <- sigma_sims[t, ]
    gamma_temp <- gamma_sims[t, , ]
    
    z_star <- HMM.viterbi(obs, mu_temp, sigma_temp, gamma_temp, draw = TRUE)
    y_hat[[t]] <- rnorm(N, mu_temp[z_star], sigma_temp[z_star])
  }
  
  # plot: posterior predictive check
  par(mfrow = c(1, 1))
  mu.hat <- colMeans(mu_sims)
  sigma.hat <- colMeans(sigma_sims)
  gamma.hat <- apply(gamma_sims, c(2, 3), mean)
  
  z_star <- HMM.viterbi(obs, mu.hat, sigma.hat, gamma.hat)
  indxs <- sample(length(y_hat), ndraw, replace = FALSE)
  plot(obs, type = "n", ylab = "Observations", 
       xlab = "Time", ylim = c(min(obs) - 3, max(obs) + 3))
  for (i in indxs) {
    lines(y_hat[[i]], col = "gray85")
  }
  points(1:N, obs, pch = 20, cex = 0.6, col = z_star)
  
}

#  predictive accuracy measures: AIC, BIC, DIC, WAIC (see BDA Ch.7)
HMM.stan.performance <- function(sims, obs) 
{
  mu_sims <- sims$mu
  sigma_sims <- sims$sigma
  gamma_sims <- sims$gamma
  
  mu.hat <- colMeans(mu_sims)
  sigma.hat <- colMeans(sigma_sims)
  gamma.hat <- apply(gamma_sims, c(2, 3), mean)
  
  n_iter <- nrow(mu_sims)
  K <- length(mu.hat)
  N <- length(obs)
  np <- K*(K-1) + 2*K
  
  temp <- matrix(NA, N, n_iter)
  llk.trace <- numeric(n_iter)
  cat(" ... computing predictive performance ... \n")
  for (t in 1:n_iter) {
    if((t %% (n_iter/10)) == 0) {
      cat(" ...", as.integer((t/n_iter)*100), "% \n")
    }
    mu <- mu_sims[t, ]
    sigma <- sigma_sims[t, ]
    gamma <- gamma_sims[t, , ]
    z <- HMM.viterbi(obs,  mu, sigma, gamma, draw = TRUE)
    llk.trace[t] <- -HMM.performance(obs, mu, sigma, gamma)$mllk
    temp[, t] <- dnorm(obs, mu[z], sigma[z])
  }
  
  lppd <- sum(log(rowMeans(temp)))
  llk <- -HMM.performance(obs, mu.hat, sigma.hat, gamma.hat)$mllk
  p.DIC <- 2*(llk - mean(llk.trace))
  p.WAIC <- 2*sum(log(rowMeans(temp)) - rowMeans(log(temp)))
  
  AIC <- -2*(llk - np)
  BIC <- -2*llk+np*log(N)
  DIC <- -2*llk + 2*p.DIC 
  WAIC <- lppd - p.WAIC
  
  list(mmlk = -llk, AIC = AIC, BIC = BIC, DIC = DIC, WAIC = -WAIC)
}



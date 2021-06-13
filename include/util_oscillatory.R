# - design matrix of fourier basis (no intercept)
fourierBasis <- function(omega, N) 
{
  M = length(omega)
  X = c()
  for (j in 1:M) {
    out = cbind(cos(2*pi*(1:N)*omega[j]), 
                sin(2*pi*(1:N)*omega[j]))
    X = cbind(X, out)
  }
  return(X)
}

# - simulate oscillatory time series 
generateOscillatory <- function(N, pars)
{
  omega <- pars$omega # freq
  beta <- pars$beta # linear basis
  sigma <- pars$sigma # res var 
  
  signal = fourierBasis(omega, N) %*% beta
  noise = rnorm(N, 0, sigma)
  data = signal + noise
  return(list(data = data, signal = signal))
}

# - ?
generateHierarchicalOscillatory <- function(parms)
{
  omega <- parms$omega
  omega.sd <- parms$omega.sd
  beta <- parms$beta 
  N.sbj <- parms$N.sbj
  N.all <- parms$N.all
  sigma.all <- parms$sigma.all
  
  M <- length(omega)
  simul <- list()
  
  for (sbj in 1:N.sbj) {
    
    N <- N.all[sbj]
    omega.sbj <- rnorm(1, omega, omega.sd)
    X.sbj <- fourierBasis(omega.sbj, N)
    beta.sbj <- rmvnorm(1, beta, parms$beta.sd*diag(2*M))
    signal <- X.sbj %*% t(beta.sbj)
    noise <- rnorm(N, 0, sigma.all[sbj])
    data.sbj <- signal + noise
    
    simul[[sbj]] <- list(data = data.sbj, 
                         omega = omega.sbj,
                         beta = beta.sbj,
                         sigma = sigma.all[sbj])
  }
  
  return(simul)
}


# (stan) - oscillatory model - posterior predictive plot 
pred.plot.oscillatory <- function(sims, obs, ndraw=50) 
{
  N <- length(obs)
  y_hat <- list()
  niter <- nrow(sims$omega)
  
  omega.sims <- sims$omega 
  beta.sims <- sims$beta
  sigma.sims <- sqrt(sims$sigma2)
  
  idxs <- sample(1:niter, ndraw, replace = FALSE)
  i <- 1
  cat(" ... sampling from posterior predictive... \n")
  for (t in idxs) {
    if((i %% (ndraw/10)) == 0) {
      cat(" ...", as.integer((i/ndraw)*100), "% \n")
    }
    omega.temp <- omega.sims[t, ]
    beta.temp <- beta.sims[t, ]
    sigma.temp <- sigma.sims[t]
    
    X <- fourierBasis(omega.temp, N)
    mu <- as.vector(X %*% beta.temp)
    y_hat[[i]] <- rnorm(N, mu, sigma.temp) 
    i <- i + 1
  }
  X.hat <- fourierBasis(omega.hat, N)
  mu.hat <- as.vector(X.hat %*% beta.hat)
  
  par(mfrow = c(1, 1))
  plot(obs, type = "n", ylab = "Obs", 
       xlab = "Time", ylim = c(min(obs, na.rm = T) - 3,
                               max(obs, na.rm = T) + 3))
  for (i in 1:ndraw) {
    lines((y_hat[[i]]), col = "gray85")
  }
  
  lines(1:N, mu.hat, col = "red", lwd = 1)
  points(1:N, obs, pch = 20, cex = 0.35)
}

# (stan) - oscillatory model - posterior predictive plot 
pred.plot.oscillatory.new <- function(sims, obs, ndraw=50) 
{
  N <- length(obs)
  y_hat <- list()
  niter <- nrow(sims$omega)
  
  omega.sims <- sims$omega 
  beta.sims <- sims$beta
  sigma.sims <- sims$sigma
  
  # - bayes estimates
  omega.hat <- mean(omega.sims)
  beta.hat <- apply(beta.sims, 2, mean)
  sigma.hat <- mean(sigma.sims)
  
  idxs <- sample(1:niter, ndraw, replace = FALSE)
  i <- 1
  cat(" ... sampling from posterior predictive... \n")
  for (t in idxs) {
    if((i %% (ndraw/10)) == 0) {
      cat(" ...", as.integer((i/ndraw)*100), "% \n")
    }
    omega.temp <- omega.sims[t]
    beta.temp <- beta.sims[t, ]
    sigma.temp <- sigma.sims[t]
    
    X <- fourierBasis(omega.temp, N)
    mu <- as.vector(X %*% beta.temp)
    y_hat[[i]] <- rnorm(N, mu, sigma.temp) 
    i <- i + 1
  }
  X.hat <- fourierBasis(omega.hat, N)
  mu.hat <- as.vector(X.hat %*% beta.hat)
  
  par(mfrow = c(1, 1))
  plot(obs, type = "n", ylab = "Obs", 
       xlab = "Time", ylim = c(min(obs, na.rm = T) - 3,
                               max(obs, na.rm = T) + 3))
  for (i in 1:ndraw) {
    lines((y_hat[[i]]), col = "gray85")
  }
  
  lines(1:N, mu.hat, col = "red", lwd = 1)
  points(1:N, obs, pch = 20, cex = 0.35)
}

# - ? 
pred_plot_multisubj_oscillatory <- function(sbj, sims, obs, N_all, ndraw = 50)
{
  
  y <- obs[sbj, 1:N_all[sbj]]
  N <- length(y)
  y_hat <- list()
  niter <- nrow(sims$omega)
  
  omega.sims <- sims$omega[, sbj]
  beta.sims <- sims$beta[, sbj, ]
  sigma.sims <- sims$sigma[, sbj]
  
  # - bayes estimates
  omega.hat <- mean(omega.sims)
  beta.hat <- apply(beta.sims, 2, mean)
  sigma.hat <- mean(sigma.sims)
  
  idxs <- sample(1:niter, ndraw, replace = FALSE)
  i <- 1
  cat(" ... sampling from posterior predictive... \n")
  for (t in idxs) {
    if((i %% (ndraw/10)) == 0) {
      cat(" ...", as.integer((i/ndraw)*100), "% \n")
    }
    omega.temp <- omega.sims[t]
    beta.temp <- beta.sims[t, ]
    sigma.temp <- sigma.sims[t]
    
    X <- fourierBasis(omega.temp, N)
    mu <- as.vector(X %*% beta.temp)
    y_hat[[i]] <- rnorm(N, mu, sigma.temp) 
    i <- i + 1
  }
  
  X.hat <- fourierBasis(omega.hat, N)
  mu.hat <- as.vector(X.hat %*% beta.hat)
  
  par(mfrow = c(1, 1))
  plot(y, type = "n", ylab = "Obs", 
       xlab = "Time", ylim = c(min(y, na.rm = T) - 1.5,
                               max(y, na.rm = T) + 1.5))
  for (i in 1:ndraw) {
    lines((y_hat[[i]]), col = "gray85")
  }
  
  lines(1:N, mu.hat, col = "red", lwd = 1)
  points(1:N, y, pch = 20, cex = 0.35)
}


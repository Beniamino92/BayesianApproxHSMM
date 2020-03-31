# generate gaussian HMM 
gauss.HMM.generate_sample <- function(N, parms)
{
  
  K <- parms$K
  mu <- parms$mu
  sigma <- parms$sigma
  gamma <- parms$gamma
  delta <- parms$delta
  if(is.null(delta))  delta <- solve(t(diag(K)-gamma + 1), rep(1,K))
  
  state <- numeric(N)
  state[1] <- sample(1:K, 1, prob = delta)
  for (i in 2:N) {
    state[i] <- sample(1:K, 1, prob = gamma[state[i-1], ])
  }
  x <- rnorm(N, mean = mu[state], sd = sigma[state])
  list(x = x, state = state)
}



# extract tpm and gauss parms from hmm depmix output
HMM.summary.pars <- function(fit, K)
{
  pars <- getpars(fit)
  gamma.hat <- matrix(pars[(K+1):(K+K*K)], K, K, byrow = TRUE)
  temp <- matrix(pars[(K+K*K+1):length(pars)], K, 2, byrow = TRUE)
  mu.hat <- temp[, 1]
  sigma.hat <- temp[, 2]
  list(gamma = gamma.hat, mu = mu.hat, sigma = sigma.hat)
}

# gaussian HMM likelihood
HMM.performance <- function(obs, mu, sigma, gamma)
{
  N <- length(obs)
  K <- length(mu)
  
  allprobs <- get.gauss.emission(obs, mu, sigma)
  lscale <- 0
  #foo <- solve(t(diag(M)-B + 1), rep(1,M))
  foo <- rep(1, K)
  for (i in 1:N) {
    foo <- foo%*%gamma*allprobs[i,]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
  }
  llk <- lscale
  np <- K*(K-1) + 2*K
  AIC <- -2*(llk -np)
  BIC <- -2*llk+np*log(N)
  list(mllk = -llk, AIC = AIC, BIC = BIC)
}


# forward backward messages (aggregates)
forward_backwards <- function(obs, allprobs, gamma) 
{
  
  K <- nrow(gamma)
  N <- length(obs)
  
  lalpha <- lbeta <- matrix(NA, K, N)
  delta <- rep(1, K)
  
  foo <- delta*allprobs[1, ]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[, 1] <- log(foo) + lscale
  
  for (i in 2:N) 
  {
    foo <- foo%*%gamma*allprobs[i,]
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
    lalpha[, i] <- log(foo) +lscale
  }
  llk <- lscale
  
  
  lbeta[, N] <- rep(0, K)
  foo <- rep(1/K, K)
  lscale <- log(K)
  
  for (i in (N-1):1) {
    foo <- gamma%*%(allprobs[i+1,]*foo)
    lbeta[,i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale + log(sumfoo)
  }
  
  list(la = lalpha, lb = lbeta, llk = llk)
}

# viterbi algorithm (see Zucchini et al., p.82) + plt
HMM.viterbi <- function(obs, mu, sigma, gamma, draw = FALSE, plt = FALSE)
{
  N <- length(obs)
  K <- length(mu)
  
  allprobs <- get.gauss.emission(obs, mu, sigma)
  xi <- matrix(0, N, K)
  delta <- solve(t(diag(K)-gamma + 1), rep(1,K))
  foo <- delta * allprobs[1, ]
  xi[1, ] <- foo/sum(foo)
  
  for (t in 2:N) {
    foo <- apply(xi[t-1, ] * gamma, 2, max) * allprobs[t, ]
    xi[t, ] <- foo/sum(foo)
  }
  z_star <- numeric(N)
  
  # - Resampling
  if (draw == TRUE) {
    z_star[N] <- sample(1:K, size = 1, prob = xi[N, ])
    for (t in (N-1):1) {
      z_star[t] <- sample(1:K, size = 1, prob = gamma[, z_star[t+1]] * xi[t, ])
    }
  } 
  # Maximizing
  else {
    z_star[N] <- which.max(xi[N, ])
    for (t in (N-1):1) {
      z_star[t] <- which.max(gamma[, z_star[t+1]] * xi[t, ])
    }
  }

  
  if (plt == TRUE) {
    par(mfrow = c(2, 1))
    plot(z_star, type = "l", lwd = 2, main = "", 
         xlab = "Time", ylab = "State Sequence")
    plot(obs, col = z_star, cex = 0.5, pch = 20, type = "p", 
         main = "", xlab = "Time", ylab = "Observations")
    par(mfrow = c(1, 1))
  }
  z_star
}

# full conditionals observations (see Zucchini et al., p.76)
HMM.conditionals <- function(obs, mu, sigma, gamma, xrange = NULL)
{
  
  N <- length(obs)
  K <- length(mu)
  
  delta <- solve(t(diag(K)-gamma+1),rep(1,K))
  if(is.null(xrange))
    xrange <- seq(from = qnorm(0.001, min(mu), max(sigma)),
                  to = qnorm(0.999, max(mu), max(sigma)), len = 1e2)
  
  allprobs <- get.gauss.emission(obs, mu, sigma)
  fb <- forward_backwards(obs, allprobs, gamma)
  la <- fb$la
  lb <- fb$lb
  la <- cbind(log(delta), la)
  lafact <- apply(la, 2, max)
  lbfact <- apply(lb, 2, max)
  
  w <- matrix(NA, ncol = N, nrow = K)
  for (i in 1:N) {
    foo <- (exp(la[, i] - lafact[i]) %*% gamma) * 
      exp(lb[, i] - lbfact[i])
    w[, i] <- foo/sum(foo)
  }
  allprobs <- get.gauss.emission(xrange, mu, sigma)
  cdists <- allprobs %*% w
  cdists <- cdists/apply(cdists, 2, sum)
  
  return(list(xrange = xrange, cdists = cdists))
}

# ordinary pseudo residuals + plot, qqplot and hist
#  (see Zucchini et al., p.96)
HMM.pseudo_residuals <- function(obs, mu, sigma, gamma, xrange = NULL, plt = TRUE)
{
  
  N <- length(obs)
  K <- length(mu)
  
  delta <- solve(t(diag(K)-gamma+1),rep(1,K))
  
  if(is.null(xrange))
    xrange <- seq(from = qnorm(0.001, min(mu), max(sigma)),
                  to = qnorm(0.999, max(mu), max(sigma)), len = 1e2)
  
  cdists <- HMM.conditionals(obs, mu, sigma, gamma, xrange = xrange)$cdists
  cumdists <- rbind(rep(0, N), apply(cdists, 2, cumsum))
  u <- rep(NA, N)
  for (t in 1:N) {
    u[t] <- cumdists[tail(which(xrange <= obs[t]), 1), t]
  }
  npsr <- qnorm(u)
  
  if (plt == TRUE) {
    par(mfrow = c(1, 3))
    # pseudo residuals
    plot(npsr, type = "o", pch = 20, cex = 0.5,
         xlab = "Time", ylab = "Pseudo-residuals")
    # hist pseudo residual
    grid <- seq(-4, 4, len = 1e2)
    hist(npsr, breaks = 20, xlim = c(-4, 4), main ="",
         xlab = "Pseudo-residuals",col = "grey", freq = F)
    lines(grid, dnorm(grid), col = "red", lwd = 2)
    # qq plot
    qqnorm(npsr)
    par(mfrow = c(1, 1))
  }
  
  return(npsr)
}













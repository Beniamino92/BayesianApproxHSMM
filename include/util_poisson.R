# generate HSSM with gaussian emission and poisson duration
gauss.HSMM.generate_sample <- function(N, parms)
{
  K <- parms$K
  lambda <- parms$lambda
  mu <- parms$mu
  sigma <- parms$sigma
  gamma <- parms$gamma
  delta <- parms$delta
  if(is.null(delta))  delta <- solve(t(diag(K)-gamma + 1), rep(1,K))
  
  state_seq = c()
  observations = c()
  state = 1
  while (length(state_seq) < N) {
    duration = rpois(1, lambda[state]) + 1
    state_seq = c(state_seq, rep(state, duration))
    obs  = rnorm(duration, mu[state], sigma[state])
    observations = c(observations, obs)
    probs = gamma[tail(state_seq, n = 1), ]
    state = sample(1:K, 1, prob = probs)
  }
  list(obs = observations[1:N], state = state_seq[1:N])
}

# generate HSSM with student-t emission and gamma duration
student.HSMM.generate_sample <- function(N, parms)
{
  K <- parms$K
  a <- parms$a
  b <- parms$b
  mu <- parms$mu
  df <- parms$df
  gamma <- parms$gamma
  delta <- parms$delta
  if(is.null(delta))  delta <- solve(t(diag(K)-gamma + 1), rep(1,K))
  
  state_seq = c()
  observations = c()
  state = 1
  while (length(state_seq) < N) {
    duration = ceiling(rgamma(1, a[state], b[state]))
    state_seq = c(state_seq, rep(state, duration))
    obs  = mu[state] + rt(duration, df[state])
    observations = c(observations, obs)
    probs = gamma[tail(state_seq, n = 1), ]
    state = sample(1:K, 1, prob = probs)
  }
  list(obs = observations[1:N], state = state_seq[1:N])
}



# ECM (Expected Conditional Maximisation) algorithm for approx HSMM
HSMM.ECM <- function(K, m, obs,  parms_init, niter = 1e3)
{
  
  mu <- parms_init$mu
  sigma <- parms_init$sigma
  lambda <- parms_init$lambda
  gamma <- parms_init$gamma
  
  N <- length(obs)
  M <- sum(m)
  
  mu.next = mu
  sigma.next = sigma
  lambda.next = lambda
  gamma.next = gamma
  
  for (iter in 1:niter)  {
    if((iter %% (niter/5)) == 0) {
      cat(" ...", as.integer((iter/niter)*100), "% \n")
    }
    allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
    # - transitions
    gamma.next <- update.transitions(K, m, gamma, lambda, allprobs_aggr)
    # - emissions
    fb <- forward_backwards.aggr(obs, m, allprobs_aggr, lambda,  gamma.next)
    la <- fb$la
    lb <- fb$lb
    temp <- update.gauss.emission(m, obs, fb)
    mu.next <- temp$mu.next
    sigma.next <- temp$sigma.next
    # - durations
    temp <- nlm(HSMM.mllk.lambda, log(lambda), 
                m = m, obs = obs, mu = mu.next,
                sigma = sigma.next, gamma = gamma.next)
    lambda.next <- exp(temp$estimate)
    
    mu <- mu.next
    sigma <- sigma.next
    lambda <- lambda.next
    gamma <- gamma.next
  }
  
  # z <- HSMM.viterbi(obs, m, lambda, mu, sigma, gamma)
  # lppd <- sum(dnorm(obs, mu[z], sigma[z], log = TRUE)) # log pointwise predictive density, BDA p. (168-169)
  llk <- -HSMM.mllk(obs, m, lambda,  mu, sigma, gamma) 
  
  np <- if(K < 3) 3*K else K^2 - 2*K + 3*K
  AIC <- -2*(llk -np)
  BIC <- -2*llk+np*log(N)
  
  return(list(mu = mu, sigma = sigma, lambda = lambda, m = m, 
              gamma = gamma, mllk = -llk, AIC = AIC, BIC = BIC))
}


# update gauss emission (see Zucchini et al., p.67)
update.gauss.emission <- function(m, obs, fb)
{
  K <- length(m)
  N <- length(obs)
  mu.next <- numeric(K)
  sigma.next <- numeric(K)
  
  la <- fb$la
  lb <- fb$lb
  llk <- fb$llk
  
  # evaluate without la
  cm = c(0, cumsum(m))
  
  for (j in 2:(K+1)) {
    idx = (cm[j - 1]+1):cm[j]
    if (m[j-1] == 1) {
      u = exp(la[idx, ] + lb[idx, ] - llk)
    } else {
      u = colSums(exp(la[idx, ] + lb[idx, ] - llk))
    }
    mu.next[j-1] <- sum(u*obs)/sum(u)
    sigma.next[j-1] <- sqrt(sum(((obs - mu.next[j-1])^2)*u)/sum(u))
  }
  list(mu.next = mu.next, sigma.next  = sigma.next)
}


# update transition probs
update.transitions <- function(K, m, gamma, lambda, allprobs_aggr)
{
  if (K < 3) {
    gamma.next <- matrix(0, K, K)
    gamma.next[!diag(K)] = 1
  } else {
    parvect <- gamma_mat.to.log_vec(K, gamma)
    fit.mle <-nlm(HSMM.mllk.gamma, parvect, m = m, 
                  lambda = lambda, allprobs_aggr = allprobs_aggr)
    gamma.next <- log_vec.to.gamma_mat(K, fit.mle$estimate)
  }
  return(gamma.next)
}

# forward backward messages (aggregates)
forward_backwards.aggr <- function(obs, m, allprobs_aggr, lambda, gamma) 
{
  
  M <- sum(m)
  N <- length(obs)
  lalpha <- lbeta <- matrix(NA, M, N)
  
  B <- B.matrix_poisson(K, m, gamma, lambda)
  delta <- rep(1, M)
  
  foo <- delta*allprobs_aggr[1, ]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[, 1] <- log(foo) + lscale
  
  for (i in 2:N) 
  {
    foo <- foo%*%B*allprobs_aggr[i,]
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
    lalpha[, i] <- log(foo) +lscale
  }
  llk <- lscale
  
  lbeta[, N] <- rep(0, M)
  foo <- rep(1/M, M)
  lscale <- log(M)
  
  for (i in (N-1):1) {
    foo <- B%*%(allprobs_aggr[i+1,]*foo)
    lbeta[,i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale + log(sumfoo)
  }
  
  list(la = lalpha, lb = lbeta, llk = llk)
}

#  generate gaussian emissions
get.gauss.emission <- function(obs, mu, sigma) 
{
  N <- length(obs)
  K <- length(mu)
  allprobs <- matrix(NA, N, K)
  for (j in 1:K) {
    allprobs[, j] = dnorm(obs, mu[j], sigma[j])
  }
  allprobs <- ifelse(!is.na(allprobs), allprobs, 1)
  allprobs
}

#  generate gaussian emission for state aggregates
get.gauss.emission_aggr <- function(obs, m, mu, sigma) 
{
  N <- length(obs)
  K <- length(m)
  # - emission probs
  allprobs <- get.gauss.emission(obs = obs, mu = mu, sigma = sigma)  
  # - emission aggregates
  allprobs_aggr = c()
  for (i in 1:K) {
    temp = matrix(rep(allprobs[, i], times = m[i]), nrow = N, ncol = m[i])
    allprobs_aggr = cbind(allprobs_aggr, temp)
  }
  allprobs_aggr
}


#  (approx) HSSM likelihood, only tpm as parameter
HSMM.mllk.gamma <- function(parvect, m, lambda, allprobs_aggr) 
{
  N <- nrow(allprobs_aggr)
  M <- sum(m)
  K <- length(m)
  gamma <- log_vec.to.gamma_mat(K, parvect)
  B <- B.matrix_poisson(K, m, gamma, lambda)
  
  lscale <- 0
  #foo <- solve(t(diag(M)-B + 1), rep(1,M))
  foo <- rep(1, M)
  for (i in 1:N) {
    foo <- foo%*%B*allprobs_aggr[i,]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}


# (approx) HSSM likelihood, only lambda as parameter
HSMM.mllk.lambda <- function(parvect, m, obs, mu, sigma, gamma) {
  
  N <- length(obs)
  K <- length(m)
  M <- sum(m)
  
  lambda <- exp(parvect)
  B <- B.matrix_poisson(K, m, gamma, lambda)
  allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
  
  lscale <- 0
  # foo <- solve(t(diag(M)-B + 1), rep(1,M))
  foo <- rep(1, M) # singularity problem if above
  
  for (i in 1:N) {
    foo <- foo%*%B*allprobs_aggr[i,]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}

# (approx) HSSM likelihood
HSMM.mllk <- function(obs, m, lambda, mu, sigma, gamma) 
{
  N <- length(obs)
  M <- sum(m)
  K <- length(m)
  B <- B.matrix_poisson(K, m, gamma, lambda)
  allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
  
  lscale <- 0
  #foo <- solve(t(diag(M)-B + 1), rep(1,M))
  foo <- rep(1, M)
  for (i in 1:N) {
    foo <- foo%*%B*allprobs_aggr[i,]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  
  return(mllk)
}

# ordinary pseudo residuals + plot, qqplot and hist
#  (see Zucchini et al., p.96)
HSMM.pseudo_residuals <- function(obs, m, lambda, mu, sigma, gamma, xrange = NULL, plt = TRUE)
{
  N <- length(obs)
  K <- length(m)
  
  if(is.null(xrange))
    xrange <- seq(from = qnorm(0.001, min(mu), max(sigma)),
                  to = qnorm(0.999, max(mu), max(sigma)), len = 5e2)
  
  cdists <- HSMM.conditionals(obs, m, lambda, mu, sigma, gamma,
                              xrange = xrange)$cdists
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
    hist(npsr, breaks = 15, xlim = c(-4, 4), main ="", 
         xlab = "Pseudo-residuals",col = "grey", freq = F)
    lines(grid, dnorm(grid), col = "red", lwd = 2)
    # qq plot
    qqnorm(npsr)
    par(mfrow = c(1, 1))
  }
  
  return(npsr)
}


# full conditionals observations (see Zucchini et al., p.76)
HSMM.conditionals <- function(obs, m, lambda, mu, sigma, gamma, xrange = NULL)
{
  
  N <- length(obs)
  M <- sum(m)
  K <- length(m)
  if(is.null(xrange))
    xrange <- seq(from = qnorm(0.001, min(mu), max(sigma)),
                  to = qnorm(0.999, max(mu), max(sigma)), len = 5e2)
  
  B <- B.matrix_poisson(K, m, gamma, lambda)
  delta <- solve(t(diag(M)-B+1),rep(1,M))
  allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
  fb <- forward_backwards.aggr(obs, m, allprobs_aggr, lambda,  gamma)
  la <- fb$la
  lb <- fb$lb
  la <- cbind(log(delta), la)
  lafact <- apply(la, 2, max)
  lbfact <- apply(lb, 2, max)
  
  w <- matrix(NA, ncol = N, nrow = M)
  for (i in 1:N) {
    foo <- (exp(la[, i] - lafact[i]) %*% B) * 
      exp(lb[, i] - lbfact[i])
    w[, i] <- foo/sum(foo)
  }
  allprobs <- get.gauss.emission_aggr(xrange, m, mu, sigma)
  cdists <- allprobs %*% w
  cdists <- cdists/apply(cdists, 2, sum)
  
  return(list(xrange = xrange, cdists = cdists))
}

# Viterbi algorithm (see Zucchini et al., p.82) + plt
HSMM.viterbi <- function(obs, m, lambda, mu, sigma, gamma, draw = FALSE, plt = FALSE) {
  
  N <- length(obs)
  K <- length(m)
  M <- sum(m)
  
  B <- B.matrix_poisson(K, m, gamma, lambda)
  allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
  xi <- matrix(0, N, M)
  delta <- solve(t(diag(M)-B + 1), rep(1,M))
  foo <- delta * allprobs_aggr[1, ]
  xi[1, ] <- foo/sum(foo)
  
  for (t in 2:N) {
    foo <- apply(xi[t-1, ] * B, 2, max) * allprobs_aggr[t, ]
    xi[t, ] <- foo/sum(foo)
  }
  z_star <- numeric(N)
  
  # - Resampling
  if (draw == TRUE) {
    z_star[N] <- sample(1:M, size = 1, prob = xi[N, ])
    for (t in (N-1):1) {
      z_star[t] <- sample(1:M, size = 1, prob = B[, z_star[t+1]] * xi[t, ])
    }
  }
  
  # - Maximizing
  else {
    z_star[N] <- which.max(xi[N, ])
    for (t in (N-1):1) {
      z_star[t] <- which.max(B[, z_star[t+1]] * xi[t, ])
    }
  }

  # standard state sequence (not aggregate)
  cm = c(0, cumsum(m))
  for (j in 2:(K+1)) {
    idx = (cm[j - 1]+1):cm[j]
    z_star[z_star %in% idx] = j - 1
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


# avoid over/under flow.
logsumexp <- function(vec) 
{
  maxvec = max(vec)
  return(maxvec + log(sum(exp(vec-maxvec))))
}


# c_hazard_dwell_poisson
c.hazard.dwell_poisson <- function(r, lambda)
{
  dpois(r - 1, lambda, )/(1-ppois(r-2, lambda))
}


# B_ii_poisson;
B.ii_poisson <- function(m_i, lambda_i) 
{
  B_ii = matrix(0, m_i, m_i)
  if (m_i != 1) {
    for (i in 1:(m_i - 1)) {
      B_ii[i, i+1] = 1 - c.hazard.dwell_poisson(i, lambda_i)
    }
  }
  B_ii[m_i, m_i] = 1-c.hazard.dwell_poisson(m_i, lambda_i)
  
  return(B_ii)
}


# B_ij_poisson
B.ij_poisson <- function(m_i, m_j, a_ij, lambda_i) 
{
  B_ij = matrix(0, m_i, m_j)
  for( i in 1:m_i){
    B_ij[i, 1] =  a_ij * c.hazard.dwell_poisson(i, lambda_i);
  }
  return(B_ij)
}


# B_matrix_poisson
B.matrix_poisson <- function(K, m, gamma, lambda) 
{
  M = sum(m)
  m_temp = numeric(K+1)
  B = matrix(0, M, M)
  
  m_temp[1] = 0
  for (i in 1:K) {
    m_temp[i+1] = m[i]
  }
  for(i in 1:K) {
    for(j in 1:K) {
      if(i ==j) {
        B[(sum(m_temp[1:i])+1):sum(m_temp[1:(i+1)]), 
          (sum(m_temp[1:j])+1):sum(m_temp[1:(j+1)])] =
          B.ii_poisson(m_temp[i+1], lambda[i]);
      } 
      else {
        B[(sum(m_temp[1:i])+1):sum(m_temp[1:(i+1)]), 
          (sum(m_temp[1:j])+1):sum(m_temp[1:(j+1)])] =
          B.ij_poisson(m_temp[i+1], m_temp[j+1], gamma[i,j], lambda[i]);
      }
    }
  }
  
  return(B)
}


# transform tpm gamma to vector of 
# unconstrained parameter (see Zucchini et al. p.49)
gamma_mat.to.log_vec <- function(K, gamma) 
{
  lpar = c()
  for (i in 1:K) {
    temp = log(gamma[i, -i]/(tail(gamma[i, -i], 1) +1e-30))[1:(K-2)]
    lpar = c(lpar, temp)
  }
  return(lpar)
}


# transform vec of unconstrained parameter to tpm gamma
log_vec.to.gamma_mat <- function(K, lpar) 
{
  gamma_out <- matrix(0, K, K)
  a = 1
  for (i in 1:K) {
    b = (a + (K-2) - 1)
    tau <- lpar[a:b]
    temp <- c(exp(tau)/(1 + sum(exp(tau))), 1/(1 + sum(exp(tau))))
    gamma_out[i, -i] <- temp
    a = b + 1
  }
  return(gamma_out)
}

# init ECM/MLE; gaussian emission
emis.init <- function(obs, K)
{
  clasif <- kmeans(obs, K)
  init.mu <- by(obs, clasif$cluster, mean)
  init.sigma <- by(obs, clasif$cluster, sd) 
  init.order <- order(init.mu)
  list(mu = init.mu[init.order], 
       sigma = init.sigma[init.order]) 
}

# init ECM/MLE; tpm
gamma.init <- function(K) 
{
  gamma <- matrix(0, K, K)
  gamma[!diag(K)] <- rep(1/(K-1), K)
  gamma 
}

# init all pars ECM/MLE 
# (might need to write something 'automatic' for lambda)
HSMM.init <- function(obs, K, lambda0)
{
  emis_init <- emis.init(obs, K)
  pars <- c()
  pars$lambda <- lambda0
  pars$mu <- emis_init$mu
  pars$sigma <- emis_init$sigma
  pars$gamma <- gamma.init(K)
  pars
}


#  MLE for (approx) HSSM
HSMM.mle <- function(K, obs, parms_init) 
{
  parvect.0 <- HSMM.pn2pw(K, parms_init$lambda, parms_init$mu, 
                          parms_init$sigma, parms_init$gamma)
  mle.fit <-nlm(HSMM.mllk.aux, parvect.0, obs = obs, m = m)
  pn <- HSMM.pw2pn(K, mle.fit$estimate)  
  mllk <- HSMM.mllk.aux(mle.fit$estimate, obs, m)
  np <- length(parvect.0)
  AIC <- 2*(mllk + np)
  n <- sum(!is.na(obs))
  BIC <- 2*mllk + np*log(n)
  
  list(lambda=pn$lambda, mu = pn$mu, sigma = pn$sigma, 
       gamma = pn$gamma, delta = pn$delta, 
       mllk = mllk, AIC=AIC, BIC=BIC)
}

# aux llk fun for MLE (approx) HSMM
HSMM.mllk.aux <- function(parvect, obs, m) 
{
  
  N <- length(obs)
  M <- sum(m)
  K <- length(m)
  pn <- HSMM.pw2pn(K, parvect)
  if (K == 2) {
    gamma.temp <- matrix(c(0, 1, 1, 0), K, K, byrow = T)
    B <- B.matrix_poisson(K, m, gamma.temp, pn$lambda)
  } else {
    B <- B.matrix_poisson(K, m, pn$gamma, pn$lambda)
  }
  
  allprobs_aggr <- get.gauss.emission_aggr(obs, m, pn$mu, pn$sigma)
  
  lscale <- 0
  # foo <- solve(t(diag(M)-B + 1), rep(1,M))
  foo <- rep(1, M)
  for (i in 1:N) {
    foo <- foo%*%B*allprobs_aggr[i,]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}

#  natural parameter to  working parameter (for nlm) in MLE
HSMM.pn2pw <- function(K, lambda, mu, sigma, gamma)
{
  if (K == 2)  {
    parvect <- c(log(lambda), mu, log(sigma))
    return(parvect) 
  }
  else {
    lgamma <- gamma_mat.to.log_vec(K, gamma)
    parvect <- c(log(lambda), mu, log(sigma), lgamma)
    return(parvect)
  }
}

#  working parameter to natural (for nlm) in MLE
HSMM.pw2pn <- function(K, parvect) 
{
  
  lambda <- exp(parvect[1:K])
  mu <- parvect[(K+1):(2*K)]
  sigma <- exp(parvect[(2*K+1):(3*K)])
  lgamma <- tail(parvect, K^2 - 2*K)
  if (K == 2) {
    gamma <- matrix(c(0, 1, 1, 0), K, K, byrow = T)
  } else {
    gamma <- log_vec.to.gamma_mat(K, lgamma)
  }
  #delta <- solve(t(diag(K) - gamma + 1), rep(1, K))
  list(lambda = lambda, mu = mu, sigma = sigma,
       gamma = gamma)
}



# ------ STAN utilities 


# starting values STAN. 
HSMM.init.stan <- function(K, obs, lambda0) 
{
  out <- emis.init(obs, K)
  out$lambda <- lambda0
  out
}


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
  sigma_sims <- sqrt(sims$sigma2)
  gamma_sims <- get.gamma_sims(sims$gamma)
  
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
HSMM.stan.performance <- function(sims, obs, m) 
{
  lambda_sims <- sims$lambda
  mu_sims <- sims$mu
  sigma_sims <- sqrt(sims$sigma2)
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
  sigma_sims <- sqrt(sims$sigma2)
  gamma_sims <- get.gamma_sims(sims$gamma)
  
  indxs <- sample(1:niter, ndraw, replace = FALSE)
  i <- 1
  cat(" ... sampling from posterior predictive... \n")
  for (t in indxs)  {
    if((i %% (ndraw/10)) == 0) {
      cat(" ...", as.integer((i/ndraw)*100), "% \n")
    }
    lambda_temp <- lambda_sims[t, ]
    mu_temp <- mu_sims[t, ]
    sigma_temp <- sigma_sims[t, ]
    gamma_temp <- gamma_sims[, , t]
    
    z_star <- HSMM.viterbi(obs, m, lambda_temp, mu_temp, 
                           sigma_temp, gamma_temp, draw = TRUE)
    y_hat[[i]] <- rnorm(N, mu_temp[z_star], sigma_temp[z_star])
    i <- i + 1
  }
  
  # plot: posterior predictive check
  par(mfrow = c(1, 1))
  lambda.hat <- colMeans(lambda_sims)
  mu.hat <- colMeans(mu_sims)
  sigma.hat <- colMeans(sigma_sims)
  gamma.hat <- apply(gamma_sims, c(1, 2), mean)
  
  z_star <- HSMM.viterbi(obs, m, lambda.hat, mu.hat, 
                         sigma.hat, gamma.hat, draw = FALSE)
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







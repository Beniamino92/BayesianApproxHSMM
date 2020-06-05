
## 1) We model the negative binomial in terms of phi = 1/size as this means that when phi = 0 we get a Poisson - Gamma prior - NO
## 2) We model the negative binomial in terms of phi = size as this means that when phi = 0 we get a Poisson - Inverse-Gamma prior - YEs

## The hazard function for the Negative-Binomail
c_hazard_NegBinom <- function(r, lambda, phi){
  ## Shifted so that the minimum value is 1
  return(dnbinom(r-1, mu = lambda, size = phi)/(1-pnbinom(r-2, mu = lambda, size = phi)))
}


## Density under the HSMM-HMM approximation to the Negative-Binomial dwell
NegBinom_dwell_approx <- function(r, lambda, phi, m){
  return(
    as.numeric((r <= m))*dnbinom(r-1, mu = lambda, size = phi) + 
      as.numeric(r > m)*dnbinom(m-1, mu = lambda, size = phi)*(1-c_hazard_NegBinom(m, lambda, phi))^(r-m))
}

## Sampling from the HSMM-HMM approximation to the Negative-Binomial dwell
rbinom_approx <- function(n, lambda, phi, m, trun = 100){
  return(sample(1:trun, n, replace = TRUE, prob = (NegBinom_dwell_approx(r = 1:trun, lambda, phi, m))))
}

## Generating the HSMM with Gaussian dwell and the approximate neative-binomial dwell distribution. 
HSMM_generate_NegBinom_approx <- function(N, K, m_vect ,lambda_vec, phi_vec, a_mat, mu_vec, sigma_vec){
  ## N is the number of univariate observations
  ## K is the number of states
  ## m_vect is the vector fo state aggregates
  ## lambda_vec and phi_vec are K x 1 matricies containing the parameters of the negative-binomail dwell distributions for each state 
  ## a_mat is the state change transition matrix with 0's on the diagnonal.
  ## mu_vec is the parameter of the emmision distribution for each state

  
  ## We start in state 1
  z_curr <- 1
  y <- c()
  z <- c()
  while(length(y) <= N){
    #r <- 1 + rnbinom(1, mu = lambda_vec[z_curr], size = phi_vec[z_curr] )
    r <- rbinom_approx(1, lambda = lambda_vec[z_curr], phi = phi_vec[z_curr], m_vect, trun = 100)
    y <- c(y, rnorm(r, mu_vec[z_curr], sigma_vec[z_curr]))
    z <- c(z, rep(z_curr, r))
    z_curr <- sample(1:K, 1, replace = TRUE, prob = a_mat[z_curr,])
  }
  return(list("y" = y[1:N], "z" = z[1:N]))
}

# generate gaussian HSSM with negative Binomial duration
gauss.HSMM.NegBinom.generate_sample <- function(N, parms)
  
{
  K <- parms$K
  lambda <- parms$lambda
  phi <- parms$phi
  mu <- parms$mu
  sigma <- parms$sigma
  gamma <- parms$gamma
  delta <- parms$delta
  
  if(is.null(delta))  delta <- solve(t(diag(K)-gamma + 1), rep(1,K))
  state_seq = c()
  observations = c()
  state = 1
  while (length(state_seq) < N) {
    duration = 1 + rnbinom(1, mu = lambda[state], size = phi[state])
    state_seq = c(state_seq, rep(state, duration))
    obs  = rnorm(duration, mu[state], sigma[state])
    observations = c(observations, obs)
    probs = gamma[tail(state_seq, n = 1), ]
    state = sample(1:K, 1, prob = probs)
  }
  list(obs = observations[1:N], state = state_seq[1:N])
}


# EM algorithm for approx HSMM
#HSMM.EM <- function(K, obs, m, lambda, mu, sigma, gamma, niter = 1e3)
negBinom.HSMM.ECM <- function(K, m, obs, parms_init, niter = 1e3)
{
  
  mu <- parms_init$mu
  sigma <- parms_init$sigma
  lambda <- parms_init$lambda
  phi <- parms_init$phi
  gamma <- parms_init$gamma
  
  N <- length(obs)
  M <- sum(m)
  
  mu.next <- mu
  sigma.next <- sigma
  lambda.next <- lambda
  phi.next <- phi
  gamma.next <- gamma
  
  for (iter in 1:niter)  {
    
    if((iter %% (niter/10)) == 0) {
      cat(" ...", as.integer((iter/niter)*100), "% \n")
    }
    
    allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
    
    # - transitions
    # gamma.next <- update.transitions(K, m, gamma, lambda, allprobs_aggr)
    gamma.next <- negBinom.update.transitions(K, m, gamma, lambda, phi, allprobs_aggr)
    # - emissions
    # fb <- forward_backwards.aggr(obs, m, allprobs_aggr, lambda,  gamma.next)
    fb <- negBinom.forward_backwards.aggr(obs, m, allprobs_aggr, lambda, phi,  gamma.next)
    la <- fb$la
    lb <- fb$lb
    
    temp <- negBinom.update.gauss.emission(m, obs, fb)
    mu.next <- temp$mu.next
    sigma.next <- temp$sigma.next
    
    # - durations
    #temp <- nlm(HSMM.mllk.lambda, log(lambda), 
    #            m = m, obs = obs, mu = mu.next,
    #            sigma = sigma.next, gamma = gamma.next)
    #lambda.next <- exp(temp$estimate)
    temp <- nlm(negBinom.HSMM.mllk.lambda_phi, c(log(lambda), log(phi)), 
                m = m, obs = obs, mu = mu.next,
                sigma = sigma.next, gamma = gamma.next)
    lambda.next <- exp(temp$estimate[1:K])
    phi.next <- exp(temp$estimate[(K+1):(2*K)])
    
    mu <- mu.next
    sigma <- sigma.next
    lambda <- lambda.next
    phi <- phi.next
    gamma <- gamma.next
  }
  
  # z <- HSMM.viterbi(obs, m, lambda, mu, sigma, gamma)
  # lppd <- sum(dnorm(obs, mu[z], sigma[z], log = TRUE)) # log pointwise predictive density, BDA p. (168-169)
  #llk <- -HSMM.mllk(obs, m, lambda,  mu, sigma, gamma) 
  llk <- -negBinom.HSMM.mllk(obs, m, lambda, phi,  mu, sigma, gamma) 
  
  #np <- if(K < 3) 3*K else K^2 - 2*K + 3*K
  np <- if(K < 3) 4*K else K^2 - 2*K + 4*K
  AIC <- -2*(llk - np)
  BIC <- -2*llk + np*log(N)
  
  #return(list(mu = mu, sigma = sigma, lambda = lambda,
  return(list(mu = mu, sigma = sigma, lambda = lambda, phi = phi,
              gamma = gamma, mllk = -llk, AIC = AIC, BIC = BIC))
}


# update gauss emission (see Zucchini et al., p.67)
negBinom.update.gauss.emission <- function(m, obs, fb)
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
    #u = colSums(exp(la[idx, ] + lb[idx, ] - llk))
    ## I edited this to deal with m == 1
    if(m[j-1] == 1){
      u <- (exp(la[idx, ] + lb[idx, ] - llk))
    }
    else{
      u <- colSums(exp(la[idx, ] + lb[idx, ] - llk))
    }
    
    mu.next[j-1] <- sum(u*obs)/sum(u)
    sigma.next[j-1] <- sqrt(sum(((obs - mu.next[j-1])^2)*u)/sum(u))
  }
  list(mu.next = mu.next, sigma.next  = sigma.next)
}


# update transition probs
#update.transitions <- function(K, m, gamma, lambda, allprobs_aggr)
negBinom.update.transitions <- function(K, m, gamma, lambda, phi, allprobs_aggr)
{
  if (K < 3) {
    gamma.next <- matrix(0, K, K)
    gamma.next[!diag(K)] = 1
  } else {
    parvect <- gamma_mat.to.log_vec(K, gamma)
    #fit.mle <-nlm(HSMM.mllk.gamma, parvect, m = m, 
    #              lambda = lambda, allprobs_aggr = allprobs_aggr)
    fit.mle <-nlm(negBinom.HSMM.mllk.gamma, parvect, m = m, 
                  lambda = lambda, phi = phi, allprobs_aggr = allprobs_aggr)
    gamma.next <- negBinom.log_vec.to.gamma_mat(K, fit.mle$estimate)
  }
  return(gamma.next)
}

# forward backward messages (aggregates)
# forward_backwards.aggr <- function(obs, m, allprobs_aggr, lambda, gamma) 
negBinom.forward_backwards.aggr <- function(obs, m, allprobs_aggr, lambda, phi, gamma) 
{
  
  M <- sum(m)
  N <- length(obs)
  lalpha <- lbeta <- matrix(NA, M, N)
  
  #B <- B.matrix_poisson(K, m, gamma, lambda)
  B <- B_matrix_negativeBinomial(K, m_vect = m, a_mat = gamma, lambda, phi)
  #delta <- solve(t(diag(M)-B+1), rep(1, M))
  delta <- rep(1, M)
  ## mu and sigma aren't inputs allprobs_aggr are! This was an errro and cuases the code not to converge anymore!
  # allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
  
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



#  (approx) HSSM likelihood, only tpm (gamma) as parameter
# HSMM.mllk.gamma <- function(parvect, m, lambda, allprobs_aggr) 
negBinom.HSMM.mllk.gamma <- function(parvect, m, lambda, phi, allprobs_aggr)
{
  N <- nrow(allprobs_aggr)
  M <- sum(m)
  K <- length(m)
  gamma <- log_vec.to.gamma_mat(K, parvect)
  #B <- B.matrix_poisson(K, m, gamma, lambda)
  B <- B_matrix_negativeBinomial(K, m_vect = m, a_mat = gamma, lambda, phi)
  
  lscale <- 0
  foo <- solve(t(diag(M)-B + 1), rep(1,M))
  #foo <- rep(1e-10, M)
  for (i in 1:N) {
    foo <- foo%*%B*allprobs_aggr[i,]
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}


# (approx) HSSM likelihood, only theta_raw (lambda) as parameter
#HSMM.mllk.lambda <- function(parvect, m, obs, mu, sigma, gamma) {
negBinom.HSMM.mllk.lambda_phi <- function(parvect, m, obs, mu, sigma, gamma) {  
  
  N <- length(obs)
  K <- length(m)
  M <- sum(m)
  
  #lambda <- exp(parvect)
  lambda <- exp(parvect[1:K])
  phi <- exp(parvect[(K+1):(2*K)])
  #B <- B.matrix_poisson(K, m, gamma, lambda)
  B <- B_matrix_negativeBinomial(K, m_vect = m, a_mat = gamma, lambda, phi)
  allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
  
  lscale <- 0
  # foo <- solve(t(diag(M)-B + 1), rep(1,M))
  foo <- rep(1/M, M) # singularity problem if above
  
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
# HSMM.mllk <- function(obs, m, lambda, mu, sigma, gamma) 
negBinom.HSMM.mllk <- function(obs, m, lambda, phi, mu, sigma, gamma) 
{
  N <- length(obs)
  M <- sum(m)
  K <- length(m)
  #B <- B.matrix_poisson(K, m, gamma, lambda)
  B <- B_matrix_negativeBinomial(K, m_vect = m, a_mat = gamma, lambda, phi)
  allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
  
  lscale <- 0
  foo <- solve(t(diag(M)-B + 1), rep(1,M))
  #foo <- rep(1, M)
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
#HSMM.pseudo_residuals <- function(obs, m, lambda, mu, sigma, gamma,
negBinom.HSMM.pseudo_residuals <- function(obs, m, lambda, phi, mu, sigma, gamma,
                                  delta = NULL, xrange = NULL, plt = TRUE)
{
  N <- length(obs)
  K <- length(m)
  
  if(is.null(delta))
    delta <- solve(t(diag(K)-gamma+1),rep(1,K))
  if(is.null(xrange))
    xrange <- seq(from = qnorm(0.001, min(mu), max(sigma)),
                  to = qnorm(0.999, max(mu), max(sigma)), len = 5e2)
  
  #cdists <- HSMM.conditionals(obs, m, lambda, mu, sigma, gamma,
  cdists <- negBinom.HSMM.conditionals(obs, m, lambda, phi, mu, sigma, gamma,
                                    delta = NULL, xrange = xrange)$cdists
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
}


# full conditionals observations (see Zucchini et al., p.76)
#HSMM.conditionals <- function(obs, m, lambda, mu, sigma, gamma, delta = NULL, xrange = NULL)
negBinom.HSMM.conditionals <- function(obs, m, lambda, phi, mu, sigma, gamma, delta = NULL, xrange = NULL)
{
  
  N <- length(obs)
  M <- sum(m)
  K <- length(m)
  if(is.null(delta))
    delta <- solve(t(diag(K)-gamma+1),rep(1,K))
  if(is.null(xrange))
    xrange <- qnorm(0.001, min(mu), 
                    max(sigma)):qnorm(0.999, max(mu), max(sigma))
  
  #B <- B.matrix_poisson(K, m, gamma, lambda)
  B <- B_matrix_negativeBinomial(K, m_vect = m, a_mat = gamma, lambda, phi)
  allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
  #fb <- forward_backwards.aggr(obs, m, allprobs_aggr, lambda,  gamma)
  fb <- negBinom.forward_backwards.aggr(obs, m, allprobs_aggr, lambda, phi,  gamma)
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
#HSMM.viterbi <- function(obs, m, lambda, mu, sigma, gamma, draw = FALSE, plt = FALSE) {
negBinom.HSMM.viterbi <- function(obs, m, lambda, phi, mu, sigma, gamma, draw = FALSE, plt = FALSE) {
  
  N <- length(obs)
  K <- length(m)
  M <- sum(m)
  
  ## B <- B.matrix_poisson(K, m, gamma, lambda)
  B <- B_matrix_negativeBinomial(K, m_vect = m, a_mat = gamma, lambda, phi)
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
  
  # - Maximizing
  if (draw == TRUE) {
    z_star[N] <- sample(1:M, size = 1, prob = xi[N, ])
    for (t in (N-1):1) {
      z_star[t] <- sample(1:M, size = 1, prob = B[, z_star[t+1]] * xi[t, ])
    }
  }
  # - Resampling
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

# posterior predictive check: draw from posterior predictive + plot
negBinom.HSMM.predictive.plot <- function(sims, obs, m, ndraw = 100)
{
  y_hat <- list()
  niter <- nrow(sims$mu)
  
  lambda_sims <- sims$lambda
  phi_sims <- 1/sims$lambda
  mu_sims <- sims$mu
  sigma_sims <- sqrt(sims$sigma2)
  gamma_sims <- get.gamma_sims(sims$gamma)
  
  indxs <- sample(1:niter, ndraw, replace = FALSE)
  i <- 1
  cat(" ... sampling from posterior predictive... \n")
  for (t in indxs) {
    if((i %% (ndraw/10)) == 0) {
      cat(" ...", as.integer((i/ndraw)*100), "% \n")
    }
    lambda_temp <- lambda_sims[t, ]
    phi_temp <- phi_sims[t, ]
    mu_temp <- mu_sims[t, ]
    sigma_temp <- sigma_sims[t, ]
    gamma_temp <- gamma_sims[, , t]
    
    z_star <- negBinom.HSMM.viterbi(obs, m, lambda_temp, phi_temp, 
                                    mu_temp, sigma_temp, gamma_temp, 
                                    draw = TRUE)
    y_hat[[i]] <- rnorm(N, mu_temp[z_star], sigma_temp[z_star])
    i <- i + 1
  }

  # plot: posterior predictive check
  par(mfrow = c(1, 1))
  lambda.hat <- colMeans(lambda_sims)
  phi.hat <- colMeans(phi_sims)
  mu.hat <- colMeans(mu_sims)
  sigma.hat <- colMeans(sigma_sims)
  gamma.hat <- apply(gamma_sims, c(1, 2), mean)
  
  z_star <- negBinom.HSMM.viterbi(obs, m, lambda.hat, phi.hat, mu.hat, 
                                  sigma.hat, gamma.hat, draw = FALSE)
  

  plot(obs, type = "n", ylab = "Observations", 
       xlab = "Time", ylim = c(min(obs) - 3, max(obs) + 3))
  for (i in 1:ndraw) {
    lines(y_hat[[i]], col = "gray85")
  }
  points(1:N, obs, pch = 20, cex = 0.6, col = z_star)
}




# B_ii_negativeBinomial;
B_ii_negativeBinomial <- function(m_i, lambda_i, phi_i){
  B_ii <- matrix(0, nrow = m_i, ncol = m_i);
  if(m_i > 1){
    for(i in 1:(m_i-1)){
      B_ii[i, i+1] = 1 - c_hazard_NegBinom(i, lambda_i, phi_i);
    }
  }
  B_ii[m_i, m_i] = 1 - c_hazard_NegBinom(m_i, lambda_i, phi_i);
  return(B_ii);
}

# B_ij_negativeBinomial;
B_ij_negativeBinomial <- function(m_i, m_j, a_ij, lambda_i, phi_i){
  B_ij <- matrix(0, nrow = m_i, ncol = m_j);
  for( i in 1:m_i){
    B_ij[i, 1] =  a_ij * c_hazard_NegBinom(i, lambda_i, phi_i);
  }
  return(B_ij);
}

# B_matrix_negativeBinomial;
B_matrix_negativeBinomial <- function(K, m_vect, a_mat, lambda_vec, phi_vec){
  sum_m = sum(m_vect);
  m_vect_temp <- rep(NA, K+1)
  B = matrix(0, nrow = sum_m, ncol = sum_m);
  m_vect_temp[1] = 0;
  for(i in 1:K){
    m_vect_temp[i+1] = m_vect[i];
  }
  for(i in 1:K){
    for(j in 1:K){
      if(i ==j){
        B[(sum(m_vect_temp[1:i])+1):sum(m_vect_temp[1:(i+1)]), (sum(m_vect_temp[1:j])+1):sum(m_vect_temp[1:(j+1)])] =
          B_ii_negativeBinomial(m_vect_temp[i+1], lambda_vec[i], phi_vec[i]);
      } 
      else{
        B[(sum(m_vect_temp[1:i])+1):sum(m_vect_temp[1:(i+1)]), (sum(m_vect_temp[1:j])+1):sum(m_vect_temp[1:(j+1)])] =
          B_ij_negativeBinomial(m_vect_temp[i+1], m_vect_temp[j+1], a_mat[i,j], lambda_vec[i], phi_vec[i]);
      }
    }
  }
  return(B);
}


# transform tpm gamma to vector of 
# unconstrained parameter (see Zucchini et al. p.49)
negBinom.gamma_mat.to.log_vec <- function(K, gamma) 
{
  lpar = c()
  for (i in 1:K) {
    temp = log(gamma[i, -i]/(tail(gamma[i, -i], 1) +1e-30))[1:(K-2)]
    lpar = c(lpar, temp)
  }
  return(lpar)
}

gamma_simplex_to_gamma_matrix <- function(gamma_simplex)
{
  p <- nrow(gamma_simplex)
  gamma_matrix <- matrix(0, nrow = p, ncol = p)
  for(i in 1:p){
    ind <- 0
    for(j in 1:p){
      if(i ==j){next}
      ind <- ind + 1
      gamma_matrix[i, j] <- gamma_simplex[i, ind]
    }
  }
  return(gamma_matrix)
}
# 
# geomDwellMean_mean_variance <- function(alpha_i, beta_i){
#   return(list("mean" = (alpha_i + beta_i -1)/(beta_i -1),
#               "var" = (alpha_i + beta_i -1)*(alpha_i + beta_i -2)/((beta_i -1)*(beta_i -2)) - ((alpha_i + beta_i -1)/(beta_i -1))^2))
# }
# 
# geomDwell_mean_var_target_error <- function(alpha_i, beta_i, mean_target, var_target){
#   temp <- geomDwellMean_mean_variance(alpha_i, beta_i)
#   return(abs(mean_target - temp$mean) + abs(var_target - temp$var))
# }
# 
# NBDwellMean_mean_variance <- function(a_0, b_0){
#   return(list("mean" = a_0/b_0, "var" = a_0/(b_0^2)))
# }
# 
# NBDwell_mean_var_target_error <- function(a_0, b_0, mean_target, var_target){
#   temp <- NBDwellMean_mean_variance(a_0, b_0)
#   return(abs(mean_target - temp$mean) + abs(var_target - temp$var))
# }



# transform vec of unconstrained parameter to tpm gamma
negBinom.log_vec.to.gamma_mat <- function(K, lpar) 
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





# -- For comparable priors: 


# mean and variance for neg-binom dwell 
nb.dwell.info <- function(a, b)
{
  list("mean" = a/b + 1, "var" = a/(b^2))
}

# error function for optim (neg-binom dwell)
nb.dwell.error<- function(a, b, mean.target, var.target)
{
  temp <- nb.dwell.info(a, b)
  abs(mean.target - temp$mean) + abs(var.target - temp$var)
}





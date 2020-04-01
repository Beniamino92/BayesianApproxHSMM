## Generate from the exact unstructured Geometric
runstructuredGeometric <- function(n, theta){
  d <- length(theta)
  samp <- sample(1:(d + 1), n, replace = TRUE, prob = c(theta, 1 - sum(theta)))
  samp[samp > d] <- d + 1 + rgeom(sum(samp > d), prob = sum(theta))
  return(samp)
}

## Density under the HSMM-HMM approximation to the geometric dwell
dunstructuredGeometric2 <- function(r, theta){
  if(sum(theta)>1){return("Error sum(theta)>1")}
  d <- length(theta)
  unstructured_start <- as.numeric(r <= d)*theta[r]
  unstructured_start[is.na(unstructured_start)==TRUE] <- 0
  geometric_tail <-  as.numeric(r > d)*(theta[d]*(1-theta[d]/(1-sum(theta[-d])))^(r-d))
  return(unstructured_start + geometric_tail)
}

## Sampling from the HSMM-HMM approximation to the geometric dwell
runstructuredGeometric2_approx <- function(n, theta, trun = 100){
  return(sample(1:trun, n, replace = TRUE, prob = (dunstructuredGeometric2(1:trun, theta))))
}

## Generating the HSMM with Gaussian dwell and the approximate unstructured geometric dwell distribution. 
HSMM_generate_unstructuredGeometric_approx <- function(N, K ,theta_list, a_mat, mu_vec, sigma_vec){
  ## N is the number of univariate observations
  ## K is the number of states
  ## theta_list is a list of K vectors each of different length giving the theta parameters of the unstructured Geometric dwell distribution  parameter of the dwell distribution for each state
  ## a_mat is the State change transition matrix with 0's on the diagnonal.
  ## mu_vec is the parameter of the emmision distribution for each state
  
  ## We start in state 1
  z_curr <- 1
  y <- c()
  z <- c()
  while(length(y) <= N){
    #r <- 1 + rpois(1, lambda = lambda_vec[z_curr])
    #r <- runstructuredGeometric(1, theta_list[[z_curr]])
    r <- runstructuredGeometric2_approx(1, theta_list[[z_curr]], trun = 100)
    y <- c(y, rnorm(r, mu_vec[z_curr], sigma_vec[z_curr]))
    z <- c(z, rep(z_curr, r))
    z_curr <- sample(1:K, 1, replace = TRUE, prob = a_mat[z_curr,])
  }
  return(list("y" = y[1:N], "z" = z[1:N]))
}


# EM algorithm for approx HSMM
#HSMM.EM <- function(K, obs, m, lambda, mu, sigma, gamma, niter = 1e3)
HSMM.ECM.geom <- function(K, m, obs, parms_init, niter = 1e3)
{
  
  mu <- parms_init$mu
  sigma <- parms_init$sigma
  theta_simplex_vec <- parms_init$theta_simplex_vec
  gamma <- parms_init$gamma
  
  N <- length(obs)
  M <- sum(m)
  
  mu.next <- mu
  sigma.next <- sigma
  #lambda.next <- lambda
  theta_simplex_vec.next <- theta_simplex_vec
  gamma.next <- gamma
  
  for (iter in 1:niter)  {
    
    if((iter %% (niter/10)) == 0) {
      cat(" ...", as.integer((iter/niter)*100), "% \n")
    }
    
    allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
    
    # - transitions
    # gamma.next <- update.transitions(K, m, gamma, lambda, allprobs_aggr)
    gamma.next <- update.transitions.geom(K, m, gamma, theta_simplex_vec, allprobs_aggr)
    # - emissions
    # fb <- forward_backwards.aggr(obs, m, allprobs_aggr, lambda,  gamma.next)
    fb <- forward_backwards.aggr.geom(obs, m, allprobs_aggr, theta_simplex_vec,  gamma.next)
    la <- fb$la
    lb <- fb$lb
    
    temp <- update.gauss.emission(m, obs, fb)
    mu.next <- temp$mu.next
    sigma.next <- temp$sigma.next
    
    # - durations
    #temp <- nlm(HSMM.mllk.lambda, log(lambda), 
    #            m = m, obs = obs, mu = mu.next,
    #            sigma = sigma.next, gamma = gamma.next)
    #lambda.next <- exp(temp$estimate)
    theta_raw_vec <- simplex_uncreate_vect(theta_simplex_vec, m_vect = m)
    temp <- nlm(HSMM.mllk.theta_simplex_vec, theta_raw_vec, 
                m = m, obs = obs, mu = mu.next,
                sigma = sigma.next, gamma = gamma.next)
    theta_simplex_vec.next <- simplex_create_vect(temp$estimate, m_vect = m)
    
    mu <- mu.next
    sigma <- sigma.next
    # lambda <- lambda.next
    theta_simplex_vec <- theta_simplex_vec.next
    gamma <- gamma.next
  }
  
  # z <- HSMM.viterbi(obs, m, lambda, mu, sigma, gamma)
  # lppd <- sum(dnorm(obs, mu[z], sigma[z], log = TRUE)) # log pointwise predictive density, BDA p. (168-169)
  #llk <- -HSMM.mllk(obs, m, lambda,  mu, sigma, gamma) 
  llk <- -HSMM.mllk.geom(obs, m, theta_simplex_vec,  mu, sigma, gamma) 
  
  #np <- if(K < 3) 3*K else K^2 - 2*K + 3*K
  # K x means, K x variances, K x (K-2) for the transition matrix, sum(m) for the theta's
  np <- K + K + K*(K-2) + sum(m)
  AIC <- -2*(llk - np)
  BIC <- -2*llk + np*log(N)
  
  #return(list(mu = mu, sigma = sigma, lambda = lambda,
  return(list(mu = mu, sigma = sigma, theta_simplex_vec = theta_simplex_vec,
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
update.transitions.geom <- function(K, m, gamma, theta_simplex_vec, allprobs_aggr)
{
  if (K < 3) {
    gamma.next <- matrix(0, K, K)
    gamma.next[!diag(K)] = 1
  } else {
    parvect <- gamma_mat.to.log_vec(K, gamma)
    #fit.mle <-nlm(HSMM.mllk.gamma, parvect, m = m, 
    #              lambda = lambda, allprobs_aggr = allprobs_aggr)
    fit.mle <-nlm(HSMM.mllk.gamma.geom, parvect, m = m, 
                  theta_simplex_vec = theta_simplex_vec, allprobs_aggr = allprobs_aggr)
    gamma.next <- log_vec.to.gamma_mat(K, fit.mle$estimate)
  }
  return(gamma.next)
}

# forward backward messages (aggregates)
# forward_backwards.aggr <- function(obs, m, allprobs_aggr, lambda, gamma) 
forward_backwards.aggr.geom <- function(obs, m, allprobs_aggr, theta_simplex_vec, gamma) 
{
  
  M <- sum(m)
  N <- length(obs)
  lalpha <- lbeta <- matrix(NA, M, N)
  
  #B <- B.matrix_poisson(K, m, gamma, lambda)
  B <- B_matrix_unstructuredGeometric(K, m_vect = m, a_mat = gamma, theta_simplex_vec)
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
  allprobs <- get.gauss.emission(obs, mu, sigma) 
  # - emission aggregates
  allprobs_aggr = c()
  for (i in 1:K) {
    temp = matrix(rep(allprobs[, i], times = m[i]), nrow = N, ncol = m[i])
    allprobs_aggr = cbind(allprobs_aggr, temp)
  }
  allprobs_aggr
}


#  (approx) HSSM likelihood, only tpm (gamma) as parameter
# HSMM.mllk.gamma <- function(parvect, m, lambda, allprobs_aggr) 
HSMM.mllk.gamma.geom <- function(parvect, m, theta_simplex_vec, allprobs_aggr)
{
  N <- nrow(allprobs_aggr)
  M <- sum(m)
  K <- length(m)
  gamma <- log_vec.to.gamma_mat(K, parvect)
  #B <- B.matrix_poisson(K, m, gamma, lambda)
  B <- B_matrix_unstructuredGeometric(K, m_vect = m, a_mat = gamma, theta_simplex_vec)
  
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
HSMM.mllk.theta_simplex_vec <- function(parvect, m, obs, mu, sigma, gamma) {  
  
  N <- length(obs)
  K <- length(m)
  M <- sum(m)
  
  #lambda <- exp(parvect)
  theta_simplex_vec <- simplex_create_vect(theta_raw = parvect, m_vect = m)
  #B <- B.matrix_poisson(K, m, gamma, lambda)
  B <- B_matrix_unstructuredGeometric(K, m_vect = m, a_mat = gamma, theta_simplex_vec)
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
HSMM.mllk.geom <- function(obs, m, theta_simplex_vec, mu, sigma, gamma) 
{
  N <- length(obs)
  M <- sum(m)
  K <- length(m)
  #B <- B.matrix_poisson(K, m, gamma, lambda)
  B <- B_matrix_unstructuredGeometric(K, m_vect = m, a_mat = gamma, theta_simplex_vec)
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
HSMM.pseudo_residuals.geom <- function(obs, m, theta_simplex_vec, mu, sigma, gamma,
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
  cdists <- HSMM.conditionals.geom(obs, m, theta_simplex_vec, mu, sigma, gamma,
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
  
  return(npsr)
}


# full conditionals observations (see Zucchini et al., p.76)
#HSMM.conditionals <- function(obs, m, lambda, mu, sigma, gamma, delta = NULL, xrange = NULL)
HSMM.conditionals.geom <- function(obs, m, theta_simplex_vec, mu, sigma, gamma, delta = NULL, xrange = NULL)
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
  B <- B_matrix_unstructuredGeometric(K, m_vect = m, a_mat = gamma, theta_simplex_vec)
  allprobs_aggr <- get.gauss.emission_aggr(obs, m, mu, sigma)
  #fb <- forward_backwards.aggr(obs, m, allprobs_aggr, lambda,  gamma)
  fb <- forward_backwards.aggr.geom(obs, m, allprobs_aggr, theta_simplex_vec,  gamma)
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
HSMM.viterbi.geom <- function(obs, m, theta_simplex_vec, mu, sigma, gamma, draw = FALSE, plt = FALSE) {
  
  N <- length(obs)
  K <- length(m)
  M <- sum(m)
  
  ## B <- B.matrix_poisson(K, m, gamma, lambda)
  B <- B_matrix_unstructuredGeometric(K, m_vect = m, a_mat = gamma, theta_simplex_vec)
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


# avoid over/under flow.
logsumexp <- function(vec) 
{
  maxvec = max(vec)
  return(maxvec + log(sum(exp(vec-maxvec))))
}


# The logit function
logit <- function(p){
  return(log(p/(1-p)))
}
# The invsere-logit function
inv_logit <- function(alpha){
  return(exp(alpha)/(exp(alpha) + 1))
}

# Creating a simplex from real numbers using centered stick breaking
simplex_create_stickBreaking <- function(theta_raw){
  ## Input: theta_raw is a vector d x 1 of real values 
  ## Output: theta_simplex is a vector (d + 1) x 1 simplex
  
  d <- length(theta_raw)
  theta_simplex <- rep(NA, d+1)
  stick_len <- 1
  for(k in 1:d){
    prop <- inv_logit(theta_raw[k] - log(d - k + 1))
    theta_simplex[k] <- stick_len * prop
    stick_len <- stick_len - theta_simplex[k]
  }
  theta_simplex[d + 1] <- stick_len 
  
  return(theta_simplex)
}

# Finding the real numbers used to create a simplex from centered stick breaking
simplex_uncreate_stickBreaking <- function(theta_simplex){
  ## Input: theta_simplex is a vector (d + 1) x 1 simplex
  ## Output: theta_raw is a vector d x 1 of real values 
  
  d <- length(theta_simplex) - 1
  theta_raw <- rep(NA, d)
  stick_len <- 0
  stick_len <- theta_simplex[d + 1]
  ####
  for(k in d:1){
    stick_len <- stick_len + theta_simplex[k]
    prop <- theta_simplex[k] / stick_len
    theta_raw[k] <- logit(prop) + log(d - k + 1)
  }
  
  return(theta_raw)
}


# Creating a vector of simplexes from a vector of real numbers
simplex_create_vect <- function(theta_raw, m_vect){
  ## Input: 
  ## m_vect is a K x 1 vector of the number of state aggregates
  ## theta_raw is a vector sum(m_vect) x 1 of real values to produce K simplexes
  ## Output: theta_simplex is a vector (sum(m_vect) + K) x 1 values which as K simplexs concatinated.
  
  K <- length(m_vect)
  theta_simplex <- rep(NA, sum(m_vect) + K)
  pos <- 0
  ## k = 1
  theta_simplex[1:(m_vect[1] + 1)] <- simplex_create_stickBreaking(theta_raw[1:m_vect[1]])
  pos = pos + (m_vect[1] + 1)
  ## k > 1
  for(k in 2:K){
    theta_simplex[(pos + 1):(pos + m_vect[k] + 1)] <- simplex_create_stickBreaking(theta_raw[(sum(m_vect[1:(k-1)]) + 1):sum(m_vect[1:k])])
    pos = pos + (m_vect[k] + 1)
  }
  return(theta_simplex)
}

# Finding the vector of real numbers used to create the vector of simplexes
simplex_uncreate_vect <- function(theta_simplex, m_vect){
  ## Input: 
  ## m_vect is a K x 1 vector of the number of state aggregates
  ## theta_simplex is a vector (sum(m_vect) + K) x 1 values which as K simplexs concatinated.
  ## Output: theta_raw is a vector sum(m_vect) x 1 of real values to produce K simplexes
  
  K <- length(m_vect)
  theta_raw <- rep(NA, sum(m_vect))
  pos <- 0
  ## k = 1
  theta_raw[1:m_vect[1]] <- simplex_uncreate_stickBreaking(theta_simplex[1:(m_vect[1] + 1)])
  pos = pos + (m_vect[1] + 1)
  ## k > 1
  for(k in 2:K){
    theta_raw[(sum(m_vect[1:(k-1)]) + 1):sum(m_vect[1:k])] <- simplex_uncreate_stickBreaking(theta_simplex[(pos + 1):(pos + m_vect[k] + 1)])
    pos = pos + (m_vect[k] + 1)
  }
  return(theta_raw)
}


## Test ##
#set.seed(8)
#m_vect <- c(3,5,4)
#theta_raw <- rnorm(sum(m_vect), 0, 1)
#theta_simplex <- simplex_create_vect(theta_raw, m_vect)
#theta_raw2 <- simplex_uncreate_vect(theta_simplex, m_vect)
#theta_raw  
#theta_simplex
#sum(theta_simplex[1:4])
#sum(theta_simplex[5:10])
#sum(theta_simplex[11:15])
#theta_raw2


# pmf of the unstructured geometric
unstructuredGeometric_pmf <- function(r, theta_simplex){
  ## r is the observation
  ## theta_simplex is a d+1 dimensional simplex
  
  d <- length(theta_simplex - 1)
  if( r <= d){
    pmf = theta_simplex[r];
  }
  else{
    pmf = (1-theta_simplex[d+1])*theta_simplex[d+1]^(r-d);
  }
  return(pmf);
}

# cmf of the unstructured geometric
unstructuredGeometric_cmf <- function(r, theta_simplex){
  ## r is the observation
  ## theta_simplex is a d+1 dimensional simplex
  
  d <- length(theta_simplex - 1)
  if( r <= 0){
    cmf = 0;
  } else{
    if( r <= d){
      cmf = sum(theta_simplex[1:r]);
    }
    else{
      cmf = (1-theta_simplex[d+1]) + (theta_simplex[d+1])*(1-theta_simplex[d+1]^(r-d));
    }
  }
  return(cmf)
}

# Hazard function of the unstructured geometric
c_hazard_unstructuredGeometric <- function(r, theta_simplex){
  ## r is the observation
  ## theta is a d+1 dimensional simplex
  
  return(unstructuredGeometric_pmf(r, theta_simplex)/(1-unstructuredGeometric_cmf(r-1, theta_simplex)))
}


# B_ii_unstructuredGeometric;
B_ii_unstructuredGeometric <- function(m_i, theta_simplex_i){
  B_ii <- matrix(0, nrow = m_i, ncol = m_i);
  if(m_i > 1){
    for(i in 1:(m_i-1)){
      B_ii[i, i+1] = 1 - c_hazard_unstructuredGeometric(i, theta_simplex_i);
    }
  }
  B_ii[m_i, m_i] = 1 - c_hazard_unstructuredGeometric(m_i, theta_simplex_i);
  return(B_ii);
}

# B_ij_unstructuredGeometric;
B_ij_unstructuredGeometric <- function(m_i, m_j, a_ij, theta_simplex_i){
  B_ij <- matrix(0, nrow = m_i, ncol = m_j);
  for( i in 1:m_i){
    B_ij[i, 1] =  a_ij * c_hazard_unstructuredGeometric(i, theta_simplex_i);
  }
  return(B_ij);
}

# B_matrix_unstructuredGeometric;
B_matrix_unstructuredGeometric <- function(K, m_vect, a_mat, theta_simplex_vec){
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
          B_ii_unstructuredGeometric(m_vect_temp[i+1], theta_simplex_vec[(sum(m_vect_temp[1:i]) + i):(sum(m_vect_temp[1:(i+1)]) + i)]);
      } 
      else{
        B[(sum(m_vect_temp[1:i])+1):sum(m_vect_temp[1:(i+1)]), (sum(m_vect_temp[1:j])+1):sum(m_vect_temp[1:(j+1)])] =
          B_ij_unstructuredGeometric(m_vect_temp[i+1], m_vect_temp[j+1], a_mat[i,j], theta_simplex_vec[(sum(m_vect_temp[1:i]) + i):(sum(m_vect_temp[1:(i+1)]) + i)]);
      }
    }
  }
  return(B);
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

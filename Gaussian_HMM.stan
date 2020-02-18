functions {
  
  /** Function: getEmissions; generate a N * K matrix with
  Gaussian emissions for t=1,..,N and j=1,..,K.
  */
  matrix getGaussianEmissions(int N, int K, real[] y, vector mu) {
    
    matrix[N, K] allprobs;
    for (n in 1:N) 
      for(k in 1:K) 
        allprobs[n, k] = exp(normal_lpdf(y[n] | mu[k], 1));
        // Excpetions about NaN?
    return allprobs;
    
  }
  
  
  
  /** Function: getForward; perform forward algorithm to 
  generate alpha_t (j) for t=1,..,N and j=1,..,K
  */
  matrix getForward(int N, int K, matrix emissions, 
                    matrix gamma_mat) {
    
    matrix[K, N] lalpha;     
    vector[K] foo;
    real sumfoo; 
    real lscale;
    
    // alpha_1
    for (k in 1:K) {
      foo[k] = emissions[1, k];
    }
    sumfoo = sum(foo);
    lscale = log(sumfoo);
    foo = foo/sumfoo; 
    lalpha[:, 1] = log(foo) + lscale;
    
    // alpha_t, t = 2, ..., N
    for (i in 2:N) {
      foo = (foo'*gamma_mat .* emissions[i, :])'; 
      sumfoo = sum(foo);
      lscale = lscale + log(sumfoo);
      foo = foo/sumfoo;
      lalpha[:, i] = log(foo) + lscale;
    }
    return lalpha;
  }
  
  /** Function: gammaToMat; convert simplex
  gamma to a K*K matrix. 
  */
  matrix gammaToMat(int K, vector[] gamma) {
    
    matrix[K, K] gamma_mat; 
    // Creating Gamma Matrix for Forward
    for(i in 1:K) 
      for(j in 1:K) 
        gamma_mat[i, j] = gamma[i][j];
    return gamma_mat;
  }
  
  
  /** Function: llk_lp; generate the likelihood p(y_1, .., y_T  | mu, gamma)
  via forward dynamic programming - p.243 (Zucchini et al.).
  */
  real llk_lp(int N, int K, real[] y, vector mu, vector[] gamma) {

    matrix[N, K] emissions = getGaussianEmissions(N, K, y, mu);
    matrix[K, K] gamma_mat = gammaToMat(K, gamma);
    matrix[K, N] lalpha =  getForward(N, K, emissions, gamma_mat); // need to store matrix ?
    real llk = log_sum_exp(lalpha[:, N]);

    return llk;
  }
  
    /** Function: ViterbiAlgorithm; generate the most likely state sequence
  via the Viterbi algorithm (with Gaussian emissions)
  */
  int[] ViterbiAlgorithm(int N, int K, real[] y, vector mu, vector[] gamma) {
    
    int z_star[N];
    int back_ptr[N, K];
    
    {
      real log_p_z_star;
      real best_logp[N, K];
      
      for (k in 1:K) 
        best_logp[1, k] = normal_lpdf(y[1] | mu[k], 1);
      for (t in 2:N) {
        for (k in 1:K) {
          best_logp[t, k] = negative_infinity();
          for (j in 1:K) {
            real logp;
            logp = best_logp[t-1, j] + log(gamma[j, k]) + normal_lpdf(y[t] | mu[k], 1);
            if (logp > best_logp[t, k]) {
              back_ptr[t, k] = j;
              best_logp[t, k] = logp;
            }
          }
        }
      }
      log_p_z_star = max(best_logp[N]);
      for (k in 1:K)
        if (best_logp[N, k] == log_p_z_star)
          z_star[N] = k;
      for (t in 1:(N-1))
        z_star[N-t] = back_ptr[N-t+1, z_star[N-t+1]];
    }
    return z_star;
  }
}


data {
  int<lower=0> N; // length of time series
  int<lower=0> K; // number of states
  real y[N]; // data
}


parameters {
  // transition probabilities
  simplex[K] gamma[K];
  // emission parameters
  positive_ordered[K] mu; 
}



model {
  // priors
  target+= normal_lpdf(mu[1] | 3, 1);
  target+= normal_lpdf(mu[2] | 10, 1);
  // likelihood
  target+= llk_lp(N, K, y, mu, gamma);
}

generated quantities {
  // generate hidden state sequence
  int<lower=1,upper=K> z_star[N] = ViterbiAlgorithm(N, K, y, mu, gamma);
}

functions {
  
    /** Function: getEmissions; generate a N * K matrix with
  Gaussian emissions for t=1,..,N and j=1,..,K.
  */
  matrix getGaussianEmissions(int N, int K, real[] y, vector[] mu, 
                              vector sigma, real omega) 
  {
    real mu_temp; 
    matrix[N, K] allprobs;
    for (n in 1:N) {
      for(k in 1:K) {
        mu_temp = mu[k, 1] + mu[k, 2]*cos(2*pi()*omega*n) + 
                  mu[k, 3]*sin(2*pi()*omega*n);
        allprobs[n, k] = exp(normal_lpdf(y[n] | mu_temp, sigma[k]));
      }
    }
    return allprobs;
  }
  
  
    /** perform forward algorithm to 
  generate alpha_t (j) for t=1,..,N and j=1,.., M = sum(m)
  via forward dynamic programming - p.243 (Zucchini et al.).
  */
  real forwardMessages(int N, int K, matrix emissions, matrix gamma_mat) 
  {
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
    // alpha_t, t = 2, ..., N
    for (i in 2:N) {
      foo = (foo'*gamma_mat .* emissions[i, :])'; 
      sumfoo = sum(foo);
      lscale = lscale + log(sumfoo);
      foo = foo/sumfoo;
    }
    return lscale;
  }
  
    /** Function: gammaToMat; convert simplex
  gamma to a K*K matrix. 
  */
  matrix gammaToMat(int K, vector[] gamma) 
  {
    matrix[K, K] gamma_mat; 
    // Creating Gamma Matrix for Forward
    for(i in 1:K) 
      for(j in 1:K) 
        gamma_mat[i, j] = gamma[i][j];
    return gamma_mat;
  }
  
    /** compute likelihood p(y_1, .., y_T  | )
  */
  real llk_lp(int N, int K, real[] y, vector[] mu, 
              vector sigma, vector[] gamma, real omega)
  {
    matrix[N, K] emissions = getGaussianEmissions(N, K, y, mu, sigma, omega);
    matrix[K, K] gamma_mat = gammaToMat(K, gamma);
    real llk =  forwardMessages(N, K, emissions, gamma_mat);
    return llk;
  } 
}


data {
  int<lower=0> N; // length time series
  int<lower=0> K; // number of states
  real y[N]; // data 
  real<lower=0, upper=0.5> omega; // frequency
  
  real mu_0[K]; // prior mean gauss emis
  real sigma_0; // prior sd      ""
  //real<lower = 0> alpha_0; // prior dirichlet probs
  vector<lower = 0>[K] alpha_0[K]; // prior dirichlet probs
}


parameters {
  simplex[K] gamma[K];// transition prob mat // what if K = 2?!
  ordered[K] mu_first; // mean gauss emission
  vector[2] mu_rest[K]; // mean gauss emission - rest of the dimension 
  vector<lower=0>[K] sigma2; // sd gauss emission
}

transformed parameters { 
  vector[3] mu[K];
  for(k in 1:K){
    mu[k] = append_row(mu_first[k], mu_rest[k]);
  }
}


model {
  // priors
  target += inv_gamma_lpdf(sigma2 | 2, 0.5);
  
  for(k in 1:K){
    target += dirichlet_lpdf(gamma[k] | alpha_0[k]);
    target += normal_lpdf(to_vector(mu[k]) | mu_0, sigma_0);
  }
  // likelihood
  target+= llk_lp(N, K, y, mu, sqrt(sigma2), gamma, omega);
}






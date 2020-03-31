
functions {
  
  // The unstructures Geometric pdf - theta is a d+1 simplex
  real unstructuredGeometric_pmf(int r, int d, vector theta){
    real pmf;
    if( r <= d){
      pmf = theta[r];
    }
    else{
      pmf = (1-theta[d+1])*theta[d+1]^(r-d);
    }
    return pmf;
  }
  
  // The unstructures Geometric cdf
  real unstructuredGeometric_cmf(int r, int d, vector theta){
    real cmf;
    if( r <= 0){
      cmf = 0;
    } else{
      if( r <= d){
        cmf = sum(theta[1:r]);
      }
      else{
        cmf = (1-theta[d+1]) + (theta[d+1])*(1-theta[d+1]^(r-d));
      }
    }
    return cmf;
  }
  
   /** c_harzard dwell unstrGeom
  */
  real c_hazard_dwell_unstructuredGeometric(int r, int d, vector theta){
    return unstructuredGeometric_pmf(r, d, theta)/(1-unstructuredGeometric_cmf(r-1, d, theta));
  }
  
  /* Functions to generate simplexes from real numbers
   - Takes a m vector of reals and turns this into a m+1 dim simplex
   - Uses centered stick-breaking process
  */

  // Creates one simplex
  vector simplex_create(vector theta_raw, int m){
    vector[m+1] theta_simplex;
    real stick_len = 1;
    real prop;
    for(j in 1:m){
      prop = inv_logit(theta_raw[j] - log(m - j + 1));
      theta_simplex[j] = stick_len * prop;
      stick_len = stick_len - theta_simplex[j];
    }
    theta_simplex[m + 1] = stick_len;
    
    return theta_simplex;
  }
  
  // Takes a sum(m_k) vector of reals and turns this into a vector of K m_k + 1 dim simplex
  vector simplex_create_vect(vector theta_raw, int[] m_vect, int K){
    vector[sum(m_vect) + K] theta_simplex;
    int pos = 0;
    // k = 1
    theta_simplex[1:(m_vect[1] + 1)] = simplex_create(theta_raw[1:m_vect[1]], m_vect[1]);
    pos  += (m_vect[1] + 1);
    // k > 1
    for(k in 2:K){
      theta_simplex[(pos + 1):(pos + m_vect[k] + 1)] = simplex_create(theta_raw[(sum(m_vect[1:(k-1)]) + 1):sum(m_vect[1:k])], m_vect[k]);
      pos += (m_vect[k] + 1);
    }
    return theta_simplex;
  }
  
  // To correctly specify we also need the log absolute Jacobian determinant of the simplex_create
  real simplex_create_lj(vector theta_raw, int m){
    real lj = 0;
    real stick_len = 1;
    real adj_theta_raw_j;
    for (j in 1:m) {
      adj_theta_raw_j = theta_raw[j] - log(m - j + 1);
      lj = lj + log(stick_len) - log1p_exp(-adj_theta_raw_j) - log1p_exp(adj_theta_raw_j);
      stick_len = stick_len * (1.0 - inv_logit(adj_theta_raw_j));
    }
    
    return lj;
  }
  
  /* Creating the B-Matrix
  */
  
  
  /** B_ii unstrGeom
  */
  matrix B_ii_unstructuredGeometric(int m_i, vector theta_simplex_i){
      matrix[m_i, m_i] B_ii;
      B_ii = rep_matrix(0, m_i, m_i);
      for(i in 1:(m_i-1)){
        B_ii[i, i+1] = 1 - c_hazard_dwell_unstructuredGeometric(i, m_i, theta_simplex_i);
      }
      B_ii[m_i, m_i] = 1 - c_hazard_dwell_unstructuredGeometric(m_i, m_i, theta_simplex_i);
      return B_ii;
  }
  
  /** B_ij  unstrGeom
  */
  matrix B_ij_unstructuredGeometric(int m_i, int m_j, real a_ij, vector theta_simplex_i){
      matrix[m_i, m_j] B_ij;
      B_ij = rep_matrix(0, m_i, m_j);
      for( i in 1:m_i){
        B_ij[i, 1] =  a_ij * c_hazard_dwell_unstructuredGeometric(i, m_i, theta_simplex_i);
      }
      return B_ij;
  }
  
  /** B_matrix unstrGeom
  */
  matrix B_matrix_unstructuredGeometric(int K, int[] m_vect, matrix a_mat, vector theta_simplex_vec){
      int sum_m = sum(m_vect);
      matrix[sum_m, sum_m] B; 
      int m_vect_temp[K+1];
      B = rep_matrix(0, sum_m, sum_m);
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
      return B;
    }
    
  /** convert (K-1) simplex gamma to a K*K matrix transition matrix.
  */
  matrix a_to_mat(int K, vector[] a)  
  {
    matrix[K, K] A = rep_matrix(0, K, K); 
    int count_j = 0;
    for (i in 1:K) {
      count_j = 0;
      for (j in 1:K) {
        if(i != j) {
          count_j += 1;
          A[i, j] = a[i][count_j]; 
        }
      }
    }
    return A;
  }
  
  /** generate a N * sum(m) matrix with Gaussian emissions
  for t=1,..,N and j=1,..,sum(m)..
 */
  matrix get_gauss_emissions(int N, int K, int[] m,  real[] y,
                             vector mu, vector sigma) 
  {
    int sum_m = sum(m);
    matrix[N, sum_m] allprobs;
    for (n in 1:N) {
      for(k in 1:K) {
        if(k ==1) {
          for (i in 1:m[k]) {
            allprobs[n, i] = exp(normal_lpdf(y[n] | mu[k], sigma[k]));
          }
        }
        else {
          for(i in 1:m[k]) {
            allprobs[n, sum(m[1:(k-1)]) + i] = exp(normal_lpdf(y[n] | mu[k], sigma[k]));
          }
        }
      }
    }
    return allprobs;
  }
  
  /** perform forward algorithm to 
  generate alpha_t (j) for t=1,..,N and j=1,.., M = sum(m)
  via forward dynamic programming - p.243 (Zucchini et al.).
  */
  real forward_messages(int N, int M, matrix emissions, matrix B) 
  {
    vector[M] foo;
    real sumfoo; 
    real lscale;
    // alpha_1
    for (k in 1:M) {
      foo[k] = emissions[1, k];
    }
    sumfoo = sum(foo);
    lscale = log(sumfoo);
    foo = foo/sumfoo; 
    // alpha_t, t = 2, ..., N
    for (i in 2:N) {
      foo = (foo'*B .* emissions[i, :])'; 
      sumfoo = sum(foo);
      lscale = lscale + log(sumfoo);
      foo = foo/sumfoo;
    }
    return lscale;
  }
  
  /** compute likelihood p(y_1, .., y_T  | )
  */
  real llk_lp(int N, int K, int[] m, real[] y, vector mu, 
              vector sigma, vector theta_simplex_vec, vector[] a)
  {
    int M = sum(m);
    matrix[K, K] A =  a_to_mat(K, a);
    matrix[M, M] B = B_matrix_unstructuredGeometric(K, m, A, theta_simplex_vec);
    matrix[N, M] allprobs = get_gauss_emissions(N, K, m, y, mu, sigma);
    real llk =  forward_messages(N, M, allprobs, B);
    return llk;
  } 
}

data {
  int<lower=0> N; // length time series
  int<lower=0> K; // number of states
  int m[K]; // size state aggr
  real y[N]; // data
  
  real mu_0[K]; // prior mean gauss emis
  real sigma_0; // prior sd      ""        ""
  real<lower = 0> alpha_0_gamma; // prior dirichlet transition probs
  real<lower = 0> alpha_0_theta; // prior dirichlet unstr-Geom params
}

parameters {
  simplex[K-1] gamma[K];// transition prob mat // what if K = 2?!
  vector[sum(m)] theta; // params of the unstr-Geom dwell distribution - each d+1 simplex only have d degrees of freedom
  ordered[K] mu; // mean gauss emission
  vector<lower=0>[K] sigma; // sd gauss emission
}

model {
  
  // First we turn our vector of theta's into a vector of simplexs
  vector[sum(m) + K] theta_simplex;
  theta_simplex = simplex_create_vect(theta, m, K);
  
  // priors
  // target += cauchy_lpdf(sigma | 0, 1);
  target += gamma_lpdf(sigma | 1, 0.5);
  // target += normal_lpdf(mu | mu_0, sigma_0*sigma);
  target += normal_lpdf(mu | mu_0, sigma_0);
  for(i in 1:K){
    target += dirichlet_lpdf(gamma[i]|rep_vector(alpha_0_gamma, K-1));
  }
  
  // Now the dirichelt prior on the simplexes not forgetting the jacobian as our parameters are the unconstrained!
  // k = 1
  target += dirichlet_lpdf(theta_simplex[1:(m[1]+1)] | rep_vector(alpha_0_theta, (m[1]+1))) + simplex_create_lj(theta[1:m[1]], m[1]);
  // k > 1
  for(k in 2:K){
    target += dirichlet_lpdf(theta_simplex[(sum(m[1:(k-1)]) + k):(sum(m[1:(k)]) + k)] | rep_vector(alpha_0_theta, (m[k]+1))) + simplex_create_lj(theta[(1+sum(m[1:(k-1)])):sum(m[1:k])], m[k]);
  }
  
  // likelihood
  target += llk_lp(N, K, m, y, mu, sigma, theta_simplex, gamma);
}







functions {
  
   /** c_harzard dwell pois
  */
  real c_hazard_dwell_poisson(int r, real lambda){
    // shifted 
    return exp(poisson_lpmf(r-1 | lambda))/(1-poisson_cdf(r-2, lambda));
  }
  
  /** B_ii poiss
  */
  matrix B_ii_poisson(int m_i, real lambda_i){
      matrix[m_i, m_i] B_ii;
      B_ii = rep_matrix(0, m_i, m_i);
      for( i in 1:(m_i-1)){
        B_ii[i, i+1] = 1 - c_hazard_dwell_poisson(i, lambda_i);
      }
      B_ii[m_i, m_i] = 1-c_hazard_dwell_poisson(m_i, lambda_i);
      return B_ii;
  }
  
  /** B_ij  poiss
  */
  matrix B_ij_poisson(int m_i, int m_j, real a_ij, real lambda_i){
      matrix[m_i, m_j] B_ij;
      B_ij = rep_matrix(0, m_i, m_j);
      for( i in 1:m_i){
        B_ij[i, 1] =  a_ij * c_hazard_dwell_poisson(i, lambda_i);
      }
      return B_ij;
  }
  
  /** B_matrix_poisson;
  */
  matrix B_matrix_poisson(int K, int[] m_vect, matrix a_mat, vector lambda_vec)
  {
      int sum_m = sum(m_vect);
      matrix[sum_m, sum_m] B; 
      int m_vect_temp[K+1];
      B = rep_matrix(0, sum_m, sum_m);
      m_vect_temp[1] = 0;
      for(i in 1:K){
        m_vect_temp[i+1] = m_vect[i];//Adding a 0 to this makes the indexing below easier
      }
      for(i in 1:K){
        for(j in 1:K){
          if(i ==j){
            B[(sum(m_vect_temp[1:i])+1):sum(m_vect_temp[1:(i+1)]), 
            (sum(m_vect_temp[1:j])+1):sum(m_vect_temp[1:(j+1)])] =
             B_ii_poisson(m_vect_temp[i+1], lambda_vec[i]);
          } 
          else{
            B[(sum(m_vect_temp[1:i])+1):sum(m_vect_temp[1:(i+1)]), 
            (sum(m_vect_temp[1:j])+1):sum(m_vect_temp[1:(j+1)])] =
             B_ij_poisson(m_vect_temp[i+1], m_vect_temp[j+1], a_mat[i,j], lambda_vec[i]);
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
              vector sigma, vector lambda, vector[] a)
  {
    int M = sum(m);
    matrix[K, K] A =  a_to_mat(K, a);
    matrix[M, M] B = B_matrix_poisson(K, m, A, lambda);
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
  real sigma_0; // prior sd      ""
  real<lower = 0> a_0;  // prior rate pois dwell
  real<lower = 0> b_0; //         ""
  real<lower = 0> alpha_0; // prior dirichlet probs
}

parameters {
  simplex[K-1] gamma[K];// transition prob mat // what if K = 2?!
  vector<lower=0>[K] lambda; // rate dwell distr
  ordered[K] mu; // mean gauss emission
  vector<lower=0>[K] sigma; // sd gauss emission
}

model {
  // priors
  // target += cauchy_lpdf(sigma | 0, 1);
  // target += normal_lpdf(mu | mu_0, sigma_0*sigma);
  target += gamma_lpdf(sigma | 2, 0.1);
  target += normal_lpdf(mu | mu_0, sigma_0);
  target += gamma_lpdf(lambda | a_0, b_0);
  for(i in 1:K){
    target += dirichlet_lpdf(gamma[i]|rep_vector(alpha_0, K-1));
  }
  // likelihood
  target+= llk_lp(N, K, m, y, mu, sigma, lambda, gamma);
}






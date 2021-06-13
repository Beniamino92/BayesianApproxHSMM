
functions {
  
  // phi = size of NB, phi = infty recovers the Poisson - Inverse gamma prior
  // phim1 = 1/size of NB, phim1 = 0 recovers the Poisson. - Gamma prior
  
   /** c_harzard dwell Negative-Binomial
  */
  real c_hazard_dwell_negativeBinomial(int r, real lambda, real phi){
    // shifted so that the minimum is 1
    return exp(neg_binomial_2_lpmf(r-1 | lambda, phi))/(1-neg_binomial_2_cdf(r-2, lambda, phi));
  }
  // E[NB] = lambda, var[NB] = lambda + lambda^2/phi
  
  /** B_ii Negative-Binomial
  */
  matrix B_ii_negativeBinomial(int m_i, real lambda_i, real phi_i){
      matrix[m_i, m_i] B_ii;
      B_ii = rep_matrix(0, m_i, m_i);
      for( i in 1:(m_i-1)){
        B_ii[i, i+1] = 1 - c_hazard_dwell_negativeBinomial(i, lambda_i, phi_i);
      }
      B_ii[m_i, m_i] = 1-c_hazard_dwell_negativeBinomial(m_i, lambda_i, phi_i);
      return B_ii;
  }
  
  /** B_ij  Negative-Binomial
  */
  matrix B_ij_negativeBinomial(int m_i, int m_j, real a_ij, real lambda_i, real phi_i){
      matrix[m_i, m_j] B_ij;
      B_ij = rep_matrix(0, m_i, m_j);
      for( i in 1:m_i){
        B_ij[i, 1] =  a_ij * c_hazard_dwell_negativeBinomial(i, lambda_i, phi_i);
      }
      return B_ij;
  }
  
  /** B_matrix_negativeBinomial;
  */
  matrix B_matrix_negativeBinomial(int K, int[] m_vect, matrix a_mat, vector lambda_vec, vector phi_vec)
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
             B_ii_negativeBinomial(m_vect_temp[i+1], lambda_vec[i], phi_vec[i]);
          } 
          else{
            B[(sum(m_vect_temp[1:i])+1):sum(m_vect_temp[1:(i+1)]), 
            (sum(m_vect_temp[1:j])+1):sum(m_vect_temp[1:(j+1)])] =
             B_ij_negativeBinomial(m_vect_temp[i+1], m_vect_temp[j+1], a_mat[i,j], lambda_vec[i], phi_vec[i]);
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
                             vector[] mu, vector sigma, real omega) 
  {
    int sum_m = sum(m);
    matrix[N, sum_m] allprobs;
    for (n in 1:N) {
      for(k in 1:K) {
        if(k ==1) {
          for (i in 1:m[k]) {
            allprobs[n, i] = exp(normal_lpdf(y[n] | mu[k, 1] + 
                 mu[k, 2]*cos(2*pi()*omega*n) + 
                 mu[k, 3]*sin(2*pi()*omega*n), sigma[k]));
          }
        }
        else {
          for(i in 1:m[k]) {
            allprobs[n, sum(m[1:(k-1)]) + i] = exp(normal_lpdf(y[n] | mu[k, 1] + 
                 mu[k, 2]*cos(2*pi()*omega*n) + 
                 mu[k, 3]*sin(2*pi()*omega*n), sigma[k]));
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
  real llk_lp(int N, int K, int[] m, real[] y, vector[] mu, 
              vector sigma, vector lambda, vector phi, vector[] a, real omega)
  {
    int M = sum(m);
    matrix[K, K] A =  a_to_mat(K, a);
    matrix[M, M] B = B_matrix_negativeBinomial(K, m, A, lambda, phi);
    matrix[N, M] allprobs = get_gauss_emissions(N, K, m, y, mu, sigma, omega);
    real llk =  forward_messages(N, M, allprobs, B);
    return llk;
  } 
}

data {
  int<lower=0> N; // length time series
  int<lower=0> K; // number of states
  int m[K]; // size state aggr
  real y[N]; // data
  real<lower=0, upper=0.5> omega; // frequency

  
  real mu_0[K]; // prior mean gauss emis
  real sigma_0; // prior sd      ""
  real<lower = 0> a_0_lambda[K];  // prior rate pois dwell
  real<lower = 0> b_0_lambda[K]; //         ""
  real<lower = 0> a_0_phi;  //           ""
  real<lower = 0> b_0_phi; //            ""
  //real<lower = 0> alpha_0; // prior dirichlet probs
  vector<lower = 0>[K-1] alpha_0[K]; // prior dirichlet probs
  
}

parameters {
  simplex[K-1] gamma[K];// transition prob mat // what if K = 2?!
  vector<lower=0>[K] lambda; // rate dwell distr
  // vector<lower=0>[K] phi; // overdispersion of the dwell distribution
  // vector<lower=0>[K] sqrt_phim1; // inverse of the square-root of the overdispersion of the dwell distribution
  vector<lower=0>[K] phim1;
  vector[3] mu[K]; // mean gauss emission
  positive_ordered[K] sigma2; // sd gauss emission
}


model {
  
  // transforming phi
  vector[K] phi;
  for(i in 1:K){
    // phi[i] = 1/(sqrt_phim1[i])^2;
    phi[i] = 1/(phim1[i]);
    target += normal_lpdf(to_vector(mu[i]) | mu_0, sigma_0);
  }

  // priors
  // target += cauchy_lpdf(sigma | 0, 1);
  //target += gamma_lpdf(sigma | 1, 0.5);
  target += inv_gamma_lpdf(sigma2 | 2, 0.5);
  // target += normal_lpdf(mu | mu_0, sigma_0*sigma);
  target += gamma_lpdf(lambda | a_0_lambda, b_0_lambda);
  // target += inv_gamma_lpdf(phi | a_0_phi, b_0_phi); // inverse gamma as 1/phi inflates the variance. 
  target += gamma_lpdf(phim1 | a_0_phi, b_0_phi); // inverse gamma as 1/phi inflates the variance.
  // target += normal_lpdf(sqrt_phim1| 0, 1);
  // target += normal_lpdf(phim1| 0, 1);
  for(i in 1:K){
    target += dirichlet_lpdf(gamma[i]|alpha_0[i]);
  }
  // likelihood
  target+= llk_lp(N, K, m, y, mu, sqrt(sigma2), lambda, phi, gamma, omega);
}






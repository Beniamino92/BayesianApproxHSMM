
functions {
  
   /** c_harzard dwell pois
  */
  real c_hazard_dwell_poisson(int r, real lambda){
    // shifted 
    return exp(poisson_lpmf(r-1 | lambda))/(1-poisson_cdf(r-2, lambda));
  }
  
  /** B_matrix_poisson - sparse representation of the transpose;
  */
  vector w_B_transpose_matrix_poisson_sparse(int K, int[] m_vect, matrix a_mat, vector lambda_vec)
  {
      int sum_m = sum(m_vect);
      vector[sum_m*K] w_B;// we can think about exactly getting the lengths right and indexing in if needs be
      int counter = 0;
      
      /*
      int m_vect_temp[K+1];
      m_vect_temp[1] = 0;
      for(i in 1:K){
        m_vect_temp[i+1] = m_vect[i];//Adding a 0 to this makes the indexing below easier
      }*/
      
      for(j in 1:K){
      // j_m is either 1, in which case we have the probabilities of arriving here from any other state 
        for(i in 1:K){
          if(i == j){
            // When m_j = 1 need to include the transition to itself - should this go in the id i == j?
            if(m_vect[j] == 1){
              counter += 1;
              w_B[counter] = (1 - c_hazard_dwell_poisson(m_vect[j], lambda_vec[j]));
            }
          } else{
            for(i_m in 1:m_vect[i]){
              counter += 1;
              w_B[counter] = (a_mat[i,j] * c_hazard_dwell_poisson(i_m, lambda_vec[i]));
            }
          }
        }
        if(m_vect[j] > 2){
        // j_m is in the middle, in which case we have it can only be reaced by the end state 
          for(j_m in 2:(m_vect[j] - 1)){
            counter += 1;
            w_B[counter] = (1 - c_hazard_dwell_poisson(j_m - 1, lambda_vec[j]));
          }
        }
        if(m_vect[j] > 1){
        // j_m = m_j and then it can be reached by the previous state or itself. 
          counter += 1;
          w_B[counter] =(1 - c_hazard_dwell_poisson(m_vect[j] - 1, lambda_vec[j]));
          counter += 1;
          w_B[counter] =(1 - c_hazard_dwell_poisson(m_vect[j], lambda_vec[j]));
        }
      }
      return w_B;
    }
    
  int[] v_B_transpose_matrix_poisson_sparse(int K, int[] m_vect, matrix a_mat, vector lambda_vec)
  {
      int sum_m = sum(m_vect);
      int v_B[K*sum_m];
      int counter = 0;
      
      int m_vect_temp[K+1];
      m_vect_temp[1] = 0;
      for(i in 1:K){
        m_vect_temp[i+1] = m_vect[i];//Adding a 0 to this makes the indexing below easier
      }
      
      for(j in 1:K){
      // j_m is either 1, in which case we have the probabilities of arriving here from any other state 
        for(i in 1:K){
          if(i == j){
            // When m_j = 1 need to include the transition to itself - should this go in the id i == j?
            if(m_vect_temp[j + 1] == 1){
              counter += 1;
              v_B[counter] = (sum(m_vect_temp[1:j]) + m_vect_temp[j + 1]);
            }
          } else{
            for(i_m in 1:m_vect_temp[i + 1]){
              counter += 1;
              v_B[counter] = (sum(m_vect_temp[1:i]) + i_m);
            }
          }
        }
        if(m_vect_temp[j + 1] > 2){
        // j_m is in the middle, in which case we have it can only be reaced by the end state 
          for(j_m in 2:(m_vect_temp[j + 1] - 1)){
            counter += 1;
            v_B[counter] = (sum(m_vect_temp[1:j]) + j_m - 1); 
          }
        }
        if(m_vect_temp[j + 1] > 1){
        // j_m = m_j and then it can be reached by the previous state or itself. 
          counter += 1;
          v_B[counter] =(sum(m_vect_temp[1:j]) + m_vect_temp[j + 1] - 1);
          counter += 1;
          v_B[counter] =(sum(m_vect_temp[1:j]) + m_vect_temp[j + 1]);
        }
      }
      return v_B;
    
    }
  
  int[] u_B_transpose_matrix_poisson_sparse(int K, int[] m_vect, matrix a_mat, vector lambda_vec)
  {
      int sum_m = sum(m_vect);
      int u_B[sum_m + 1];
      int counter = 0;
      int ind;
      
      counter += 1;
      u_B[counter] = 1;
      for(j in 1:K){
      // j_m is either 1, in which case we have the probabilities of arriving here from any other state 
        ind = 0;
        for(i in 1:K){
          if(i == j){
            // When m_j = 1 need to include the transition to itself - should this go in the id i == j?
            if(m_vect[j] == 1){
              ind += 1;
            }
          } else{
            for(i_m in 1:m_vect[i]){
              ind += 1;
            }
          }
        }
        counter += 1;
        u_B[counter] = (u_B[counter - 1] + ind);
      if(m_vect[j] > 2){
      // j_m is in the middle, in which case we have it can only be reaced by the end state 
        for(j_m in 2:(m_vect[j] - 1)){
          counter += 1;
          u_B[counter] = (u_B[counter - 1] + 1); // c[1:size(c)]
        }
      }
      if(m_vect[j] > 1){
        // j_m = m_j and then it can be reached by the previous state or itself. 
        counter += 1;
        u_B[counter] = (u_B[counter - 1] + 2);
      }
    }
      return u_B;
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
    real mu_temp;
    matrix[N, sum_m] allprobs;
    for (n in 1:N) {
      for(k in 1:K) {
        if(k ==1) {
          for (i in 1:m[k]) {
            mu_temp = mu[k, 1] + mu[k, 2]*cos(2*pi()*omega*n) + 
                      mu[k, 3]*sin(2*pi()*omega*n);
            allprobs[n, i] = exp(normal_lpdf(y[n] | mu_temp, sigma[k]));
          }
        }
        else {
          for(i in 1:m[k]) {
            mu_temp = mu[k, 1] + mu[k, 2]*cos(2*pi()*omega*n) + 
                      mu[k, 3]*sin(2*pi()*omega*n);
            allprobs[n, sum(m[1:(k-1)]) + i] = exp(normal_lpdf(y[n] | mu_temp, sigma[k]));
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
  real forward_messages_sparse(int N, int M, matrix emissions, vector w_B_transpose, int[] v_B_transpose, int[] u_B_transpose) 
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
      //foo = (foo'*B .* emissions[i, :])'; 
      // bA = (A^Tb^T)^T where A^T is sparse
      foo = ((csr_matrix_times_vector(M, M, w_B_transpose, v_B_transpose, u_B_transpose, foo))' .* emissions[i, :])';
      sumfoo = sum(foo);
      lscale = lscale + log(sumfoo);
      foo = foo/sumfoo;
    }
    return lscale;
  }
  
  /** compute likelihood p(y_1, .., y_T  | )
  */
  real llk_sparse_lp(int N, int K, int[] m, real[] y, vector[] mu, 
              vector sigma, vector lambda, vector[] a, real omega)
  {
    int M = sum(m);
    matrix[K, K] A =  a_to_mat(K, a);
    // bA = (A^Tb^T)^T where A^T is sparse
    vector[K*M] w_B_transpose = w_B_transpose_matrix_poisson_sparse(K, m, A, lambda);
    int v_B_transpose[K*M]= v_B_transpose_matrix_poisson_sparse(K, m, A, lambda);
    int u_B_transpose[M + 1] = u_B_transpose_matrix_poisson_sparse(K, m, A, lambda);
    matrix[N, M] allprobs = get_gauss_emissions(N, K, m, y, mu, sigma, omega);
    real llk =  forward_messages_sparse(N, M, allprobs, w_B_transpose, v_B_transpose, u_B_transpose);
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
  real<lower = 0> a_0[K];  // prior rate pois dwell
  real<lower = 0> b_0[K]; //         ""
  // real<lower = 0> alpha_0; // prior dirichlet probs
  vector<lower = 0>[K-1] alpha_0[K]; // prior dirichlet probs
}

parameters {
  simplex[K-1] gamma[K];// transition prob mat // what if K = 2?!
  vector<lower=0>[K] lambda; // rate dwell distr
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
  // target += cauchy_lpdf(sigma | 0, 1);
  // target += gamma_lpdf(sigma | 1, 0.5);
  target += inv_gamma_lpdf(sigma2 | 2, 0.5);
  // target += normal_lpdf(mu | mu_0, sigma_0*sigma);
  target += gamma_lpdf(lambda | a_0, b_0);
  for(i in 1:K){
    target += dirichlet_lpdf(gamma[i]|alpha_0[i]);
    target += normal_lpdf(to_vector(mu[i]) | mu_0, sigma_0);
  }
  // likelihood
  target+= llk_sparse_lp(N, K, m, y, mu, sqrt(sigma2), lambda, gamma, omega);
}






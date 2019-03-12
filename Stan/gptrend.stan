/*
  Stan implementation of the Trendiness of Trends
  This is the stationary fully Bayesian version. 
  The following hyper-parameters must be passed as data
  1) mean for the prior of mu
  2) shape and scale parameters for the prior of rho
*/

functions{
  matrix pred_rng(real[] xPred, vector y, real[] x, real m, real alpha, real rho, real sigma) {
    /* 
       Returns a 4 x length(xPred) matrix with the following rows:
       1) random sample from the posterior of f
       2) random sample from the posteriof of df
       3) random sample from the predictive posterior
       4) Trend Change Index
    */
    int n = rows(y);
    int nPred = size(xPred);
    matrix[2 * nPred, 2 * nPred] diag_delta = diag_matrix(rep_vector(1e-9, 2 * nPred));
    
    matrix[n, n] K;
    matrix[n, n] L;
    matrix[n, nPred] K1_f;
    matrix[nPred, nPred] K2_f;
    matrix[n, nPred] K1_df;
    matrix[nPred, nPred] K2_df;
    
    vector[n] L_div_y_f;
    matrix[n, nPred] L_div_K1_f;
    vector[n] L_div_y_df;
    matrix[n, nPred] L_div_K1_df;

    vector[nPred] mu_f;
    vector[nPred] mu_df;
    matrix[nPred, nPred] C11;
    matrix[nPred, nPred] C12;
    matrix[nPred, nPred] C21;
    matrix[nPred, nPred] C22;
    
    vector[2 * nPred] mu;
    matrix[2 * nPred, 2 * nPred] C;
    
    vector[2 * nPred] sim;
    matrix[nPred, 4] result;
    
    //Get C(x, x) + sigma^2*I
    K = cov_exp_quad(x, alpha, rho);
    for (i in 1:n) {
      K[i, i] = K[i, i] + square(sigma);
    }
    L = cholesky_decompose(K);
    
    /* 
    Stuff for f
    */
    //Get C(x, xPred)
    K1_f = cov_exp_quad(x, xPred, alpha, rho);
    //Get C(xPred, xPred)
    K2_f = cov_exp_quad(xPred, alpha, rho);      
    
    L_div_y_f = mdivide_left_tri_low(L, y - m);
    L_div_K1_f = mdivide_left_tri_low(L, K1_f);
    mu_f = m + L_div_K1_f' * L_div_y_f;
    C11 = K2_f - L_div_K1_f' * L_div_K1_f;
    
    /* 
    Stuff for df 
    */
    //Get \partial_2 C(x, xPred)
    for(i in 1:n) {
      for(j in 1:nPred) {
        K1_df[i, j] = K1_f[i, j] * (x[i] - xPred[j]) / square(rho);
      }
    }
    //Get \partial_1 \partial_2 C(xPred, xPred)
    for(i in 1:nPred) {
      for(j in 1:nPred) {
       K2_df[i, j] = K2_f[i, j] * (pow(rho, -2) - pow(rho, -4) * square(xPred[i] - xPred[j]));
      }
    }
  
    L_div_y_df = mdivide_left_tri_low(L, y);
    L_div_K1_df = mdivide_left_tri_low(L, K1_df);
    mu_df = L_div_K1_df' * L_div_y_df;
    C22 = K2_df - L_div_K1_df' * L_div_K1_df;
    
    /* 
    Get cross-covariances - fix this!
    */
    for(i in 1:nPred) {
      for(j in 1:nPred) {
        C12[i, j] = 0;
        C21[i, j] = 0;
      }
    }
    
    /* 
    Collect the joint moments
    */
    mu = append_row(mu_f, mu_df);    
    C = append_row(append_col(C11, C12), append_col(C21, C22));

    /* 
    Simulate from the joint distribution
    */
    sim = multi_normal_rng(mu, C + diag_delta);

    /* 
    Gather the results
    */
    result[1:nPred, 1] = sim[1:nPred];
    result[1:nPred, 2] = sim[(nPred + 1):(2*nPred)];
    for(i in 1:nPred) {
      result[i, 3] = normal_rng(sim[i], sigma);
      result[i, 4] = 1 - normal_cdf(0, mu_df[i], sqrt(C22[i, i]));
    }

    return result;
  }
}

data {
  int<lower=1> n;
  real x[n];
  vector[n] y;
  
  int<lower = 1> nPred;
  real xPred[nPred];
  
  //Hyper-parameters
  real mu_mu;
  real<lower = 0> rho_shape;
  real<lower = 0> rho_scale;
}

parameters {
  real mu;
  real<lower=0> alpha;
  real<lower=0> rho;
  
  real<lower=0> sigma;
}

model {
  matrix[n, n] L;
  {
    matrix[n, n] K = cov_exp_quad(x, alpha, rho);
    for (i in 1:n) {
      K[i, i] = K[i, i] + square(sigma);
    }
    L = cholesky_decompose(K);
  }
  
  mu ~ normal(mu_mu, 3);
  alpha ~ normal(0, 3);
  rho ~ inv_gamma(rho_shape, rho_scale);
  
  sigma ~ normal(0, 3);

  y ~ multi_normal_cholesky(rep_vector(mu, n), L);
}

generated quantities {
  matrix[nPred, 4] pred;
  pred = pred_rng(xPred, y, x, mu, alpha, rho, sigma);
}


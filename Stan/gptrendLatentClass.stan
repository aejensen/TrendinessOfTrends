/*
  Stan implementation of the Trendiness of Trends
  This is the nonstationary fully Bayesian version (with K = 2). 
  
  The following must be passed as data:
  1) number of degrees of freedom for the spline basis
  2) the basis evaluated at x (basis)
  3) the basis evaluted at xPred (basisPred)
  4) the basis derivatives evaluated at x (basisD)
  5) the basis derivaties evaluated at xPred (basisPredD)
  
  As well ad the following hyper-parameters:
  6) mean for the prior of mu (mu_mu)
  7) shape and scale parameters for the prior of rho (rho_shape, rho_scale)
*/

functions{
  vector inv_logit_D(vector f) {
    /*
      Returns the derivative of the inverse logit function
    */
    int n = rows(f);
    vector[n] out = exp(f) ./ square(1 + exp(f));
    return out;
  }
  
  matrix pred_rng(real[] xPred, vector y, real[] x, real m, 
                  vector w, vector wPred, vector wD, vector wPredD,
                  vector alpha, vector rho, real sigma) {
    /* 
       Returns a length(xPred) x 4 matrix with the following rows:
       1) random sample from the posterior of f
       2) random sample from the posteriof of df
       3) random sample from the predictive posterior
       4) Trend Direction Index
    */
    int n = rows(y);
    int nPred = size(xPred);
    matrix[nPred, nPred] diag_delta = diag_matrix(rep_vector(1e-9, nPred));
    
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
    {
      matrix[n, n] K_1 = cov_exp_quad(x, alpha[1], rho[1]);
      matrix[n, n] K_2 = cov_exp_quad(x, alpha[2], rho[2]);
      for(i in 1:n) {
        for(j in 1:n) {
          K[i, j] = (w[i] * w[j]) * K_1[i, j] + (1 - w[i]) * (1 - w[j]) * K_2[i, j];
        }
      }
    }
    for (i in 1:n) {
      K[i, i] = K[i, i] + square(sigma);
    }
    L = cholesky_decompose(K);
    
    /* 
    Stuff for f
    */
    //Get C(x, xPred)
    {
      matrix[n, nPred] K1_f_1 = cov_exp_quad(x, xPred, alpha[1], rho[1]);
      matrix[n, nPred] K1_f_2 = cov_exp_quad(x, xPred, alpha[2], rho[2]);
      for(i in 1:n) {
        for(j in 1:nPred) {
          K1_f[i, j] = (w[i] * wPred[j]) * K1_f_1[i, j] + (1 - w[i]) * (1 - wPred[j]) * K1_f_2[i, j];
        }
      }
    }
    
    //Get C(xPred, xPred)
    {
      matrix[nPred, nPred] K2_f_1 = cov_exp_quad(xPred, alpha[1], rho[1]);
      matrix[nPred, nPred] K2_f_2 = cov_exp_quad(xPred, alpha[2], rho[2]);
      for(i in 1:nPred) {
        for(j in 1:nPred) {
          K2_f[i, j] = (wPred[i] * wPred[j]) * K2_f_1[i, j] + (1 - wPred[i]) * (1 - wPred[j]) * K2_f_2[i, j];
        }
      }
    }
    
    //Get posterior moments of f
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
        matrix[n, nPred] K1_df_1[i, j] = K1_f[i, j] * (x[i] - xPred[j]) / square(rho[1]);
        matrix[n, nPred] K1_df_2[i, j] = K1_f[i, j] * (x[i] - xPred[j]) / square(rho[2]);
        //K1_df[i, j] = ?
      }
    }
    //Get \partial_1 \partial_2 C(xPred, xPred)
    for(i in 1:nPred) {
      for(j in 1:nPred) {
       matrix[nPred, nPred] K2_df_1[i, j] = K2_f[i, j] * (pow(rho[1], -2) - pow(rho[1], -4) * square(xPred[i] - xPred[j]));
       matrix[nPred, nPred] K2_df_2[i, j] = K2_f[i, j] * (pow(rho[2], -2) - pow(rho[2], -4) * square(xPred[i] - xPred[j]));
       //K2_df[i, j] = ?
      }
    }
  
    //Get posterior moments of df
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
  
  //Spline data
  int<lower = 1> df;
  matrix[df, n] basis;
  matrix[df, nPred] basisPred;
  matrix[df, n] basisD;
  matrix[df, nPred] basisPredD;
  
  //Hyper-parameters
  real mu_mu;
  real<lower = 0> rho_shape;
  real<lower = 0> rho_scale;
}

parameters {
  real mu;
  vector<lower = 0>[2] alpha;
  positive_ordered[2] rho;

  real<lower=0> sigma;
  
  row_vector[df] a_raw; 
  real a0; 
  real<lower =0> tau;   
}

transformed parameters { 
  row_vector[df] a; 
  vector[n] w; 
 
  a[1] = a_raw[1];
  for (i in 2:df) {
    a[i] = a[i-1] + tau * a_raw[i];
  }
    
  w = inv_logit(a0 + to_vector(a * basis)); 
} 

model {
  matrix[n, n] K1;
  matrix[n, n] K2;
  matrix[n, n] K;
  matrix[n, n] L;
  
  K1 = cov_exp_quad(x, alpha[1], rho[1]);
  K2 = cov_exp_quad(x, alpha[2], rho[2]);
  for(i in 1:n) {
    for(j in 1:n) {
      K[i, j] = (w[i] * w[j]) * K1[i, j] + (1 - w[i]) * (1 - w[j]) * K2[i, j];
    }
  }
  for (i in 1:n) {
    K[i, i] = K[i, i] + square(sigma);
  }
  L = cholesky_decompose(K);

  a0 ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  tau ~ normal(0, 1); 
  
  mu ~ normal(mu_mu, 3);
  
  alpha ~ normal(0, 3);
  rho ~ inv_gamma(rho_shape, rho_scale);
  
  sigma ~ normal(0, 3);
  
  y ~ multi_normal_cholesky(rep_vector(mu, n), L);
}

generated quantities {
  vector[nPred] wD;
  vector[nPred] wPred;
  vector[nPred] wPredD;
  matrix[nPred, 2] pred;

  wD = to_vector(a * basisD) .* inv_logit_D(a0 + to_vector(a * basis));
  wPred = inv_logit(a0 + to_vector(a * basisPred)); 
  wPredD = to_vector(a * basisPredD) .* inv_logit_D(a0 + to_vector(a * basisPred)); 

  pred = pred_rng(xPred, y, x, mu, w, wPred, wD, wPredD, alpha, rho, sigma);
}

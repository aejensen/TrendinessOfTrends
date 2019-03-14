/*
  Stan implementation of the Trendiness of Trends
  This is the stationary fully Bayesian version. 
  The following hyper-parameters must be passed as data
  1) mean for the prior of mu
  2) shape and scale parameters for the prior of rho
*/

functions{
  real generalized_inverse_gaussian_lpdf(real x, int p, real a, real b) {
    return p * 0.5 * log(a / b) - log(2 * modified_bessel_second_kind(p, sqrt(a * b))) + (p - 1) * log(x) - (a * x + b / x) * 0.5;
  }  
  
  real cov_matern_52(real s, real t, real alpha, real rho) {
    real c;
    real d = fabs(s - t);
    if(d < 1e-9) {
      c = 1;
    } else {
      c = exp(-sqrt(5) * d / rho) * (5 * square(s - t) + 3 * rho * (sqrt(5) * d + rho)) / (3 * square(rho));
    }
    return square(alpha) * c;
  }
  
  real cov_matern_52_D(real s, real t, real alpha, real rho) {
    real c;
    real d = fabs(s - t);
    if(d < 1e-9) {
      c = 0;
    } else {
      c = 5 * exp(-sqrt(5) * d / rho) * (s - t) * (sqrt(5) * d  + rho) / (3 * pow(rho, 3));
    }
    return square(alpha) * c;
  }
  
  real cov_matern_52_DD(real s, real t, real alpha, real rho) {
    real c;
    real d = fabs(s - t);
    if(d < 1e-9) {
      c = 5 / (3 * square(rho));
    } else {
      c = 5 * exp(-sqrt(5) * d / rho) * (rho * (sqrt(5) * d + rho) - 5 * square(s - t)) / (3 * pow(rho, 4));
    }
    return square(alpha) * c;
  }
  
  matrix pred_rng(real[] xPred, vector y, real[] x, real m0, real m1, real alpha, real rho, real sigma) {
    /* 
       Returns a 4 x length(xPred) matrix with the following rows:
       1) random sample from the posterior of f
       2) random sample from the posteriof of df
       3) random sample from the predictive posterior
       4) Trend Change Index
    */
    int n = rows(y);
    int nPred = size(xPred);
    matrix[2 * nPred, 2 * nPred] diag_delta = diag_matrix(rep_vector(1e-6, 2 * nPred));
    //matrix[nPred, nPred] diag_delta = diag_matrix(rep_vector(1e-9, nPred));
    
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
    //vector[nPred] sim;
    matrix[nPred, 4] result;
    
    //Get C(x, x) + sigma^2*I
    //K = cov_exp_quad(x, alpha, rho);
    for(i in 1:n) {
      for(j in 1:n) {
        K[i, j] = cov_matern_52(x[i], x[j], alpha, rho);
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
    //K1_f = cov_exp_quad(x, xPred, alpha, rho);
    for(i in 1:n) {
      for(j in 1:nPred) {
        K1_f[i, j] = cov_matern_52(x[i], xPred[j], alpha, rho);
      }
    }
    //Get C(xPred, xPred)
    //K2_f = cov_exp_quad(xPred, alpha, rho);      
    for(i in 1:nPred) {
      for(j in 1:nPred) {
        K2_f[i, j] = cov_matern_52(xPred[i], xPred[j], alpha, rho);
      }
    }
    
    L_div_y_f = mdivide_left_tri_low(L, y - (m0 + m1 * to_vector(x)));
    L_div_K1_f = mdivide_left_tri_low(L, K1_f);
    mu_f = (m0 + m1*to_vector(xPred)) + L_div_K1_f' * L_div_y_f;
    C11 = K2_f - L_div_K1_f' * L_div_K1_f;
    
    /* 
    Stuff for df 
    */
    //Get \partial_2 C(x, xPred)
    for(i in 1:n) {
      for(j in 1:nPred) {
        //K1_df[i, j] = K1_f[i, j] * (x[i] - xPred[j]) / (sqrt(3) * fabs(x[i] - xPred[j]) * rho + square(rho));
        K1_df[i, j] = cov_matern_52_D(x[i], xPred[j], alpha, rho);
      }
    }
    
    //Get \partial_1 \partial_2 C(xPred, xPred)
    for(i in 1:nPred) {
      for(j in 1:nPred) {
        K2_df[i, j] = cov_matern_52_DD(xPred[i], xPred[j], alpha, rho);
      }
    }
    
    L_div_y_df = mdivide_left_tri_low(L, y - m1);
    L_div_K1_df = mdivide_left_tri_low(L, K1_df);
    mu_df = m1 + L_div_K1_df' * L_div_y_df;
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
    //sim = multi_normal_rng(mu_f, C11);

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
  real<lower = 0> rho_shape;
  real<lower = 0> rho_scale;
}

parameters {
  real mu0;
  real mu1;
  real<lower=0> alpha;
  real<lower=0> rho;
  
  real<lower=0> sigma;
}

transformed parameters {
  vector[n] mu = mu0 + mu1 * to_vector(x);
}

model {
  matrix[n, n] L;
  {
    matrix[n, n] K;
    for(i in 1:n) {
      for(j in 1:n) {
        K[i, j] = cov_matern_52(x[i], x[j], alpha, rho);
      }
    }
    for (i in 1:n) {
      K[i, i] = K[i, i] + square(sigma);
    }
    L = cholesky_decompose(K);
  }
  
  mu0 ~ normal(26.8350000, 3);
  mu1 ~ normal(-0.6941479, 3);
  alpha ~ normal(0, 3);
  rho ~ inv_gamma(rho_shape, rho_scale);
  //rho ~ generalized_inverse_gaussian(-1, 2, 1);
  
  sigma ~ normal(0, 3);

  y ~ multi_normal_cholesky(mu, L);
}

generated quantities {
  matrix[nPred, 4] pred;
  pred = pred_rng(xPred, y, x, mu0, mu1, alpha, rho, sigma);
}


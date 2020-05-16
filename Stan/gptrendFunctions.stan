/*
  Stan implementation of the Trendiness of Trends
  This file contains functions doing most of the work
  AKJ, 2019
  
  Addendum 2020: Adding SE covariance functions
*/

functions{
  real cov_se(real s, real t, real alpha, real rho) {
    //Squared Exponential (SE) covariance function
    return square(alpha) * exp(-square(s - t) / (2 * square(rho)));
  }
  
  real cov_se_D2(real s, real t, real alpha, real rho) {
    //d_2 C(s,t)
    return cov_se(s, t, alpha, rho) * (s - t) / square(rho);
  }

  real cov_se_D2_D2(real s, real t, real alpha, real rho) {
    //d_2^2 C(s,t)
    return cov_se(s, t, alpha, rho) * (square(s-t) - square(rho)) / pow(rho, 4);
  }

  real cov_se_D1_D2(real s, real t, real alpha, real rho) {
    //d_1 d_2 C(s,t)
    return cov_se(s, t, alpha, rho) * (square(rho) - square(s - t)) / pow(rho, 4);
  }

  real cov_se_D1_D2_D2(real s, real t, real alpha, real rho) {
    //d_1 d_2^2 C(s,t)
    real k = 3 * (s - t) * square(rho) - pow(s - t, 3);
    return cov_se(s, t, alpha, rho) * k / pow(rho, 6);
  }

  real cov_se_D1_D1_D2_D2(real s, real t, real alpha, real rho) {
    //d_1^2 d_2^2 C(s,t)
    real k = pow(s - t, 4) - 6 * square(s - t) * square(rho) + 3 * pow(rho, 4);
    return cov_se(s, t, alpha, rho) * k / pow(rho, 8);
  }

  real cov_rq(real s, real t, real alpha, real rho, real nu) {
    //Rational Quadratic (RQ) covariance function
    return square(alpha) * pow(1 + square(s - t) / (2*nu*square(rho)), -nu);
  }

  real cov_rq_D2(real s, real t, real alpha, real rho, real nu) {
    //d_2 C(s,t)
    real k = (2 * (s-t) * nu) / (square(s-t) + 2 * nu * square(rho));
    return cov_rq(s, t, alpha, rho, nu) * k; 
  }

  real cov_rq_D2_D2(real s, real t, real alpha, real rho, real nu) {
    //d_2^2 C(s,t)
    real d = square(s - t);
    real num = 2 * nu * (d * (1 + 2 * nu) - 2 * nu * square(rho));
    real den = square(d + 2 * nu * square(rho));
    return cov_rq(s, t, alpha, rho, nu) * (num / den); 
  }

  real cov_rq_D1_D2(real s, real t, real alpha, real rho, real nu) {
    //d_1 d_2 C(s,t)
    real num = 4 * square(nu) * square(rho) - 2 * square(s-t) * nu * (1 + 2 * nu);
    real den = square(square(s-t) + 2 * nu * square(rho));
    return cov_rq(s, t, alpha, rho, nu) * (num / den);
  }

  real cov_rq_D1_D2_D2(real s, real t, real alpha, real rho, real nu) {
    //d_1 d_2^2 C(s,t)
    real d = square(s - t);
    real num = 4 * (s - t) * nu * (1 + nu) * (-d * (1 + 2 * nu) + 6 * nu * square(rho));
    real den = pow(d + 2 * nu * square(rho), 3);
    return cov_rq(s, t, alpha, rho, nu) * (num / den);
  }

  real cov_rq_D1_D1_D2_D2(real s, real t, real alpha, real rho, real nu) {
    //d_1^2 d_2^2 C(s,t)
    real d = square(s -t);
    real den = pow(d + 2 * nu * square(rho), 4);
    real num1 = 4 * square(d) * nu * (1 + nu) * (3 + 8 * nu + 4 * square(nu));
    real num2 = 48 * d * square(nu) * (1 + nu) * (3 + 2 * nu) * square(rho);
    real num3 = 48 * pow(nu, 3) * (1 + nu) * pow(rho, 4);
    real k = num1/den - num2/den + num3/den;
    return cov_rq(s, t, alpha, rho, nu) * k;
  }
  
  vector dnorm(vector x) {
    //Vectorized standard normal density function
    return exp(-0.5 * square(x)) / sqrt(2 * pi());
  }

  vector pnorm(vector x) {
    //Vectorized standard normal cumulative density function
    return 0.5 * (erf(x / sqrt(2)) + 1);
  }

  vector calcTDI(int p, vector mu, matrix cov) {
    /*
    Returns a p-vector of the Trend Direction Index
    mu, cov are the joint posterior covariances matrices of (f, df, d^2f)
    */
    vector[p] teddy;
    for(i in 1:p) {
      teddy[i] = 1 - normal_cdf(0, mu[p + i], sqrt(cov[p + i, p + i]));
    }
    return teddy;
  }
  
  vector calcETI(int p, vector mu, matrix cov) {
    /*
    Returns a p-vector of the local Expected Trend Instability
    mu, cov are the joint posterior covariances matrices of (f, df, d^2f)
    */
    vector[p] mu_df = mu[(p+1):(2*p)];
    vector[p] mu_ddf = mu[(2*p+1):(3*p)];
    vector[p] sigma_df = sqrt(diagonal(cov[(p+1):(2*p), (p+1):(2*p)]));
    vector[p] sigma_ddf = sqrt(diagonal(cov[(2*p+1):(3*p), (2*p+1):(3*p)]));
    vector[p] sigma_df_ddf = diagonal(cov[(p+1):(2*p), (2*p+1):(3*p)]);
  
    vector[p] omega = sigma_df_ddf ./ (sigma_df .* sigma_ddf);
    vector[p] eta = (mu_ddf - sigma_ddf .* omega .* mu_df ./ sigma_df) ./ (sigma_ddf .* sqrt(1 - square(omega)));
  
    vector[p] t1 = sigma_ddf ./ sigma_df .* sqrt(1 - square(omega));
    vector[p] t2 = dnorm(mu_df ./ sigma_df) .* (2 * dnorm(eta) + eta .* (2* pnorm(eta) - 1));
    
    vector[p] eddy = t1 .* t2;
    
    return eddy;
  }  

  matrix calcPostMoments(real[] tPred, real[] t, vector y, 
                         vector mY, vector m, vector dm, vector ddm,
                         real alpha, real rho, real nu, real sigma) {
    /*
    Returns a (3*p) x (3*p + 1) matrix of posterior moments of
    (f, df, ddf). The first column is the 3*p vector of posterior means
    */
    int n = rows(y);     //Number of observations
    int p = size(tPred); //Number of prediction points

    matrix[n, n] L;
    matrix[n, 3*p] L_div_K1;
    matrix[3*p, 3*p] K2;
    
    vector[3*p] mu_joint;
    matrix[3*p, 3*p] cov_joint;

    /*
    Calculate the posterior moments
    */
    {
      matrix[n, n] K; //C(t, t) + sigma^2 * I
      for(i in 1:n) {
        for(j in 1:n) {
          K[i, j] = cov_rq(t[i], t[j], alpha, rho, nu);
        }
      }
      for (i in 1:n) {
        K[i, i] = K[i, i] + square(sigma);
      }
      L = cholesky_decompose(K);
    }
    
    //Calculate posterior means
    {
      matrix[n, p] K1_f;   //C(t, tPred)
      matrix[n, p] K1_df;  //d_2 C(t, tPred)
      matrix[n, p] K1_ddf; //d_2^2 C(t, tPred)
      matrix[n, 3*p] K1;
      vector[3*p] m_vec;      
      for(i in 1:n) {
        for(j in 1:p) {
          K1_f[i, j] = cov_rq(t[i], tPred[j], alpha, rho, nu);
          K1_df[i, j] = cov_rq_D2(t[i], tPred[j], alpha, rho, nu);
          K1_ddf[i, j] = cov_rq_D2_D2(t[i], tPred[j], alpha, rho, nu);
        }
      }
      K1 = append_col(append_col(K1_f, K1_df), K1_ddf);
      L_div_K1 = mdivide_left_tri_low(L, K1);
      
      m_vec = append_row(append_row(m, dm), ddm);
      mu_joint = m_vec + L_div_K1' * mdivide_left_tri_low(L, y - mY); 
    }
    
    //Calculate posterior covariance
    {
      matrix[p, p] K2_11; //C(tPred, tPred)
      matrix[p, p] K2_12; //d_2 C(tPred, tPred)
      matrix[p, p] K2_13; //d_2^2 C(tPred, tPred)
      matrix[p, p] K2_22; //d_1 d_2 C(tPred, tPred)
      matrix[p, p] K2_23; //d_1 d_2^2 C(tPred, tPred)
      matrix[p, p] K2_33; //d_1^2 d_2^2 C(tPred, tPred)
      for(i in 1:p) {
        for(j in 1:p) {
          K2_11[i, j] = cov_rq(tPred[i], tPred[j], alpha, rho, nu);
          K2_12[i, j] = cov_rq_D2(tPred[i], tPred[j], alpha, rho, nu);
          K2_13[i, j] = cov_rq_D2_D2(tPred[i], tPred[j], alpha, rho, nu);
          K2_22[i, j] = cov_rq_D1_D2(tPred[i], tPred[j], alpha, rho, nu);
          K2_23[i, j] = cov_rq_D1_D2_D2(tPred[i], tPred[j], alpha, rho, nu);
          K2_33[i, j] = cov_rq_D1_D1_D2_D2(tPred[i], tPred[j], alpha, rho, nu);
        }
      }
      K2 = append_row(append_row(append_col(append_col(K2_11,  K2_12),  K2_13),
                                 append_col(append_col(K2_12', K2_22),  K2_23)),
                                 append_col(append_col(K2_13', K2_23'), K2_33));
    }
    cov_joint = K2 - L_div_K1' * L_div_K1;
    
    return append_col(mu_joint, cov_joint);
  }
  
  matrix gpFit_rng(real[] tPred, real[] t, vector y, 
                   vector mY, vector m, vector dm, vector ddm,
                   real alpha, real rho, real nu, 
                   real sigma) {
                     
    int n = rows(y);     //Number of observations
    int p = size(tPred); //Number of prediction points
    
    vector[3*p] mu_joint;
    matrix[3*p, 3*p] cov_joint;
    vector[3*p] sim;
    matrix[p, 6] result;
    
    /* 
    Get the posterior moments
    */
    {
      matrix[3*p, 3*p + 1] pm = calcPostMoments(tPred, t, y, mY, m, dm, ddm, alpha, rho, nu, sigma);
      mu_joint = pm[, 1];
      cov_joint = pm[, 2:(3*p + 1)];
    }
    
    /* 
    Simulate from the joint posterior distribution
    */
    {
      //Add a small constant to the diagonal for numerical stability
      matrix[3*p, 3*p] diag_delta = diag_matrix(rep_vector(1e-9, 3*p));
      sim = multi_normal_rng(mu_joint, cov_joint + diag_delta);
    }
    
    /*
    Gather results
    */
    result[1:p, 1] = sim[1:p];
    for(i in 1:p) {
      result[i, 2] = normal_rng(sim[i], sigma);
    }
    result[1:p, 3] = sim[(p+1):(2*p)];
    result[1:p, 4] = sim[(2*p+1):(3*p)];
    result[1:p, 5] = calcTDI(p, mu_joint, cov_joint);
    result[1:p, 6] = calcETI(p, mu_joint, cov_joint);
    
    return result;  
  }  
}

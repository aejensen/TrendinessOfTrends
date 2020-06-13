functions {
  real cov_se(real s, real t, real alpha, real rho) {
    //Squared Exponential (SE) covariance function
    return square(alpha) * exp(-square(s - t) / (2 * square(rho)));
  }
  
  real cov_se_D2(real s, real t, real alpha, real rho) {
    //d_2 C(s,t)
    return cov_se(s, t, alpha, rho) * (s - t) / square(rho);
  }
  
  real cov_se_D1_D2(real s, real t, real alpha, real rho) {
    //d_1 d_2 C(s,t)
    return cov_se(s, t, alpha, rho) * (square(rho) - square(s - t)) / pow(rho, 4);
  }
}

data {
  int<lower = 0> n;
  real t[n];
  vector[n] y;
  
  real<lower = 0, upper = 1> quantile;
}

parameters {
  //Parameters for (f, df)
  real m;
  real<lower = 0> alpha;
  real<lower = 0> rho;
  vector[2*n] eta;
  
  real<lower = 0> sigma;
  vector<lower = 0>[n] w;
}

transformed parameters {
  real<lower = 0> tau;
  vector[n] me;
  vector[n] pe;
  //vector[n] pe2;
  
  vector[n] f;
  vector[n] df;
  
  {
    vector[2*n] latent;
    matrix[2*n, 2*n] K;
    matrix[n, n] K11;
    matrix[n, n] K12;
    matrix[n, n] K22;
    
    for(i in 1:n) {
      for(j in 1:n) {
        K11[i, j] = cov_se(t[i], t[j], alpha, rho);
        K12[i, j] = cov_se_D2(t[i], t[j], alpha, rho);
        K22[i, j] = cov_se_D1_D2(t[i], t[j], alpha, rho);
      }
    }
    K = append_row(append_col(K11,  K12),
                   append_col(K12', K22));
    for (i in 1:(2*n)) {
      K[i, i] = K[i, i] + 1e-9; //numerical stability
    }
    latent = cholesky_decompose(K) * eta;
    
    f = m + latent[1:n];
    df = latent[(n+1):(2*n)];
  }  
  
  tau = pow(sigma, -2);
  me  = (1 - 2 * quantile) / (quantile * (1 - quantile)) * w + f;
  pe  = (quantile * (1 - quantile) * tau) ./ (2 * w);
  
  /*
  for(i in 1:n) {
    pe2[i] = inv_sqrt(pe[i]);
  }
  */
  //pe2 = inv_sqrt(pe);
}

model {
  m ~ normal(0, 10);
  alpha ~ normal(0, 1);
  rho ~ inv_gamma(5, 5);
  eta ~ normal(0, 1);  
  
  sigma ~ cauchy(0, 2.5);
  
  w ~ exponential(tau);
  
  y ~ normal(me, inv_sqrt(pe));
}

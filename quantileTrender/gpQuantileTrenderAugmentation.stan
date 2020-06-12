functions {
  real cov_se(real s, real t, real alpha, real rho) {
    //Squared Exponential (SE) covariance function
    return square(alpha) * exp(-square(s - t) / (2 * square(rho)));
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
  vector[n] eta;
  
  real<lower = 0> sigma;
  vector<lower = 0>[n] w; //data augmentation
}

transformed parameters {
  real<lower = 0> tau;
  vector[n] me;
  vector[n] pe;
  vector[n] pe2;
  vector[n] f;
  vector[n] mu;
  
  {
    matrix[n, n] K;
    for(i in 1:n) {
      for(j in 1:n) {
        K[i, j] = cov_se(t[i], t[j], alpha, rho);
      }
    }
    for (i in 1:n) {
      K[i, i] = K[i, i] + 1e-9;
    }
    f = cholesky_decompose(K) * eta;
  }  
  
  mu = m + f;
  tau = pow(sigma, -2);
  me  = (1 - 2 * quantile) / (quantile * (1 - quantile)) * w + mu;
  pe  = (quantile * (1 - quantile) * tau) ./ (2 * w);
  
  for(i in 1:n) {
    pe2[i] = inv_sqrt(pe[i]);
  }
}

model {
  m ~ normal(0, 10);
  alpha ~ normal(0, 1);
  rho ~ inv_gamma(5, 5);
  eta ~ normal(0, 1);  
  
  sigma ~ cauchy(0, 2.5);
  
  w ~ exponential(tau);
  
  y ~ normal(me, pe2);
}

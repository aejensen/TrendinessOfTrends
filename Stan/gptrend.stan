#include /gptrendFunctions.stan

data {
  int<lower = 1> n;  // Number of observations
  real t[n];         // Vector of sampling times
  vector[n] y;       // Vector of outcomes
  
  int<lower = 1> p;  // Number of prediction points
  real tPred[p];     // Vector of prediction points

  real mu_mu;
  real<lower = 0> alpha_mu;  
  real<lower = 0> rho_mu;
  real<lower = 0> nu_mu;
  real<lower = 0> sigma_mu;
}

parameters {
  real mu;
  
  real<lower = 0> alpha;  
  real<lower = 0> rho;
  real<lower = 0> nu;
  real<lower = 0> sigma;
}

model {
  matrix[n, n] L;
  
  {
    matrix[n, n] K;
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
  
  mu ~ normal(mu_mu, 1);
  
  alpha ~ gamma(square(alpha_mu), alpha_mu);
  
  rho ~ gamma(square(rho_mu), rho_mu);
  
  nu ~ normal(nu_mu, 0.001);
  
  sigma ~ gamma(square(sigma_mu) / 1.0, sigma_mu / 1.0);
  
  y ~ multi_normal_cholesky(rep_vector(mu, n), L);
}

generated quantities {
  vector[n] mY = rep_vector(mu, n);
  vector[p] m = rep_vector(mu, p);
  vector[p] dm = rep_vector(0, p);
  vector[p] ddm = rep_vector(0, p);
  
  matrix[p, 6] pred = gpFit_rng(tPred, t, y, mY, m, dm, ddm, alpha, rho, nu, sigma);
}
/*
  Stan implementation of the Trendiness of Trends
  This is the Bayesian version
  AKJ, 2019
*/

#include /gptrendFunctions.stan

data {
  int<lower = 1> n;  // Number of observations
  real t[n];         // Vector of sampling times
  vector[n] y;       // Vector of outcomes
  
  int<lower = 1> p;  // Number of prediction points
  real tPred[p];     // Vector of prediction points

  real mu_mu;
  
  real alpha_mu;  
  real rho_mu;
  real nu_mu;
  
  real sigma_mu;
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
  
  mu ~ student_t(3, mu_mu, 3);
  
  alpha ~ student_t(3, alpha_mu, 3);
  rho ~ normal(rho_mu, 1);
  nu ~ student_t(3, nu_mu, 3);
  
  sigma ~ student_t(3, sigma_mu, 3);
  
  y ~ multi_normal_cholesky(rep_vector(mu, n), L);
}

generated quantities {
  matrix[p, 6] pred;
  {
    vector[n] mY = rep_vector(mu, n); //mu_beta(t)
    vector[p] m = rep_vector(mu, p);  //mu_beta(tPred)
    vector[p] dm = rep_vector(0, p);  //d mu_beta(tPred)
    vector[p] ddm = rep_vector(0, p); //d^2 mu_beta(tPred)
  
    pred = gpFit_rng(tPred, t, y, mY, m, dm, ddm, alpha, rho, nu, sigma);
  }
}

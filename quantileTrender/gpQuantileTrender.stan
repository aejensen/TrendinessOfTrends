functions {
   real rho_quantile(real y, real quantile) {
     if (y < 0) {
       return y * (quantile - 1);
     } else {
       return y * quantile;
     }
   }
   
   real asym_laplace_lpdf(real y, real mu, real sigma, real quantile) { 
     return log(quantile * (1 - quantile)) - log(sigma) - 
            rho_quantile((y - mu) / sigma, quantile); 
   }
}

data {
  int<lower = 1> n;
  real t[n];
  vector[n] y;
  
  int<lower = 1> p;
  real tPred[p];
  
  real<lower=0, upper=1> quantile;
}

transformed data {
  int<lower=1> n_tot = n + p;
  
  real t_tot[n_tot];
  for (i in 1:n) {
    t_tot[i] = t[i];
  }
  for (i in 1:p) {
    t_tot[n + i] = tPred[i];
  }
}

parameters {
  real mu;
  real<lower = 0> alpha;  
  real<lower = 0> rho;
  real<lower = 0> sigma;
  
  vector[n_tot] eta;
}

transformed parameters {
  /*
    In the next iteration f should be matrix[n_tot, 2] to get
    the joint posterior of (f, df).
  */
  vector[n_tot] f;
  
  {
    matrix[n_tot, n_tot] K = cov_exp_quad(t_tot, alpha, rho);
    for (i in 1:n_tot) {
      K[i, i] = K[i, i] + 1e-9; //numerical stability
    }
    
    f = cholesky_decompose(K) * eta;
  }
}

model {
  mu ~ normal(0, 10);
  alpha ~ normal(0, 1);
  rho ~ inv_gamma(5, 5);
  sigma ~ cauchy(0, 2.5);
  
  eta ~ normal(0, 1);
  
  for (i in 1:n) {
    target += asym_laplace_lpdf(y[i] | mu + f[i], sigma, quantile);
  }
}

generated quantities {
  vector[p] fPred = mu + f[(n+1):n_tot];
}

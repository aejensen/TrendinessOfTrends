library(mvtnorm)
library(DEoptimR)
library(rstan)
library(pracma)
library(parallel)
library(rootSolve)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####################################
# Global variables
####################################
tPred <- seq(0, 1, length.out=1000)
rho <- sqrt(3) / (pi * 2)
ETI_true <- sqrt(3) / (pi * rho)

m <- stan_model("gptrendFixedSE.stan")

####################################
# Functions
####################################
seCov <- function(s, t, alpha, rho) {
  alpha^2 * exp(-(s-t)^2 / (2*rho^2))
}

optPar <- function(t, y) {
  op <- JDEoptim(function(par) {
    mu <- rep(par[1], length(t))
    cMat <- outer(t, t, seCov, par[2], par[3]) + diag(par[4] ,length(t))
    -mvtnorm::dmvnorm(y, mu, cMat, log=TRUE)
  }, lower = c(-10, 0, 0, 0), upper = c(10, 2, 2, 1), 
  trace=TRUE, triter = 1000, maxiter = 500 * 10^3, tol = 1e-10)
  op$par
}

getStan <- function(t, y, tPred, par) {
  sDat <- list(n = length(t), t = t, y = y, p = length(tPred), tPred = tPred)
  sDat$mu <- par[1]
  sDat$alpha <- par[2]
  sDat$rho <- par[3]
  sDat$sigma <- par[4]
  
  fit <- sampling(m, data = sDat, iter = 1, seed = 12345, 
                  chains = 1, algorithm = "Fixed_param", refresh = -1)
  extract(fit, "pred")$pred
}

doSim <- function(n, sigma) {
  tSeq <- seq(0, 1, length.out=n)
  
  cMat <- outer(tPred, tPred, seCov, alpha = 1, rho = rho)
  mu <- rep(0, length(tPred))
  
  f <- as.numeric(rmvnorm(1, mu, cMat))
  df <- predict(smooth.spline(tPred, f, all.knots=TRUE), deriv=1)$y
  
  y <- approxfun(tPred, f)(tSeq) + rnorm(length(tSeq), 0, sigma)
  
  opt <- optPar(tSeq, y)
  post <- getStan(tSeq, y, tPred, opt)
  
  list(t = tSeq, tPred = tPred, f = f, df = df, y = y, 
       parOpt = opt, fMu = post[1,,7], dfMu = post[1,,8], ddfMu = post[1,,9],
       TDI = post[1,,5], ETI = post[1,,6])
}

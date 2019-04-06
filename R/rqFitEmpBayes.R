rm(list=ls())

library(DEoptim)
library(mvtnorm)
library(numDeriv)
library(fields)

########################################################
# Load data
########################################################
load("../Data/smoking.RData")
tPred <- seq(1998, 2018, 0.2)
dat <- data.frame(t = smoking$year, y = smoking$p)

dat$t <- dat$t - mean(smoking$year)
tPred <- seq(1998, 2018, length.out = 500)
tPred <- tPred - mean(smoking$year)
tPredCalendar <- seq(1998, 2018, length.out = 500)

########################################################
# Define RQ covariance function and its derivatives
########################################################
rqCov <- function(s, t, alpha, rho, nu) {
  alpha^2 * (1 + (s - t)^2 / (2 * nu * rho^2))^(-nu)
}

rqCov_D1 <- function(s, t, alpha, rho, nu) {
  k <- (2 * (t - s) * nu) / ((s - t)^2 + 2 * nu * rho^2)
  rqCov(s, t, alpha, rho, nu) * k
}

rqCov_D1_D1 <- function(s, t, alpha, rho, nu) {
  d <- (s - t)^2
  k <- 2 * nu * (d * (1 + 2*nu) - 2*nu*rho^2) / (d + 2 * nu * rho^2)^2
  rqCov(s, t, alpha, rho, nu) * k
}

rqCov_D1_D2 <- function(s, t, alpha, rho, nu) {
  d <- (s - t)^2
  k <- (4 * nu^2 * rho^2 - 2 * d * nu * (1 + 2 * nu)) / (d + 2 * nu * rho^2)^2
  rqCov(s, t, alpha, rho, nu) * k
}

rqCov_D1_D1_D2 <- function(s, t, alpha, rho, nu) {
  d <- (s - t)^2
  k <- 4 * (s - t) * nu * (1 + nu) * (d * (1 + 2 * nu) - 6 * nu * rho^2) / (d + 2 * nu * rho^2)^3
  rqCov(s, t, alpha, rho, nu) * k
}

rqCov_D1_D1_D2_D2 <- function(s, t, alpha, rho, nu) {
  d <- (s-t)^2
  den <- (d + 2 * nu * rho^2)^4
  k1 <- 4 * d^2 * nu * (1 + nu) * (3 + 8 * nu + 4 * nu^2)
  k2 <- 48 * d * nu^2 * (1 + nu) * (3 + 2 * nu) * rho^2
  k3 <- 48 * nu^3 * (1 + nu) * rho^4
  k <- k1/den - k2/den + k3/den
  rqCov(s, t, alpha, rho, nu) * k
}

########################################################
# Empirical Bayes parameter estimates
########################################################
ctl <- DEoptim.control(itermax = 1000, trace = 100)
opt.rq <- DEoptim(function(par) {
  mu <- rep(par[1], nrow(dat))
  cMat <- outer(dat$t, dat$t, rqCov, par[2], par[3], par[4]) + diag(par[5]^2 , nrow(dat))
  -mvtnorm::dmvnorm(dat$y, mu, cMat, log=TRUE)
}, lower = c(0,0,0,0,0), upper = c(50,50,50,50,50), control = ctl)
par.rq <- opt.rq$optim$bestmem
par.rq

########################################################
# Get posterior moments
########################################################
getPosterior <- function(tPred, tObs, y, par) {
  n <- length(y)
  p <- length(tPred)
  
  mu <- par[1]
  alpha <- par[2]
  rho <- par[3]
  nu <- par[4]
  sigma <- par[5]
  
  K <- outer(tObs, tObs, rqCov, alpha, rho, nu) + diag(sigma^2, n)
  
  K1_f <- outer(tPred, tObs, rqCov, alpha, rho, nu)
  K1_df <- outer(tPred, tObs, rqCov_D1, alpha, rho, nu)
  K1_ddf <- outer(tPred, tObs, rqCov_D1_D1, alpha, rho, nu)
  K1 <- rbind(K1_f, K1_df, K1_ddf)
  
  #Get posterior mean
  m_vec <- matrix(c(rep(mu, p), rep(0, p), rep(0, p)), 3*p, 1)
  mu_joint <- as.numeric(m_vec + K1 %*% solve(K) %*% as.matrix(y - mu, n, 1))
  mu_f <- mu_joint[1:p]
  mu_df <- mu_joint[(p+1):(2*p)]
  mu_ddf <- mu_joint[(2*p+1):(3*p)]
  
  
  #Get posterior covariance
  K2_11 <- outer(tPred, tPred, rqCov, alpha, rho, nu)
  K2_21 <- outer(tPred, tPred, rqCov_D1, alpha, rho, nu)
  K2_22 <- outer(tPred, tPred, rqCov_D1_D2, alpha, rho, nu)
  K2_31 <- outer(tPred, tPred, rqCov_D1_D1, alpha, rho, nu)
  K2_32 <- outer(tPred, tPred, rqCov_D1_D1_D2, alpha, rho, nu)
  K2_33 <- outer(tPred, tPred, rqCov_D1_D1_D2_D2, alpha, rho, nu)
  K2_12 <- t(K2_21)
  K2_13 <- t(K2_31)
  K2_23 <- t(K2_32)
  K2 <- rbind(cbind(K2_11, K2_12, K2_13),
              cbind(K2_21, K2_22, K2_23),
              cbind(K2_31, K2_32, K2_33))
  
  v <-  rbind(K1_f, K1_df, K1_ddf)
  cov_joint <- K2 - v %*% solve(K) %*% t(v)
  
  list(mu_joint = mu_joint, 
       mu_f = mu_f, mu_df = mu_df, mu_ddf = mu_ddf,
       cov_joint = cov_joint)
}

post <- getPosterior(tPred, dat$t, dat$y, par.rq)

plot(tPred, post$mu_f, type="l")
plot(tPred, post$mu_df, type="l", ylim=c(-2,2))
abline(h = 0, lty=2)
plot(tPred, post$mu_ddf, type="l", ylim=c(-1,1))
abline(h = 0, lty=2)

#Sanity check
mean((post$mu_df - sapply(tPred, function(t) grad(splinefun(tPred, post$mu_f), t)))^2)
mean((post$mu_ddf - sapply(tPred, function(t) grad(splinefun(tPred, post$mu_df), t)))^2)

########################################################
# Simulate from the posterior
########################################################
doSim <- function(post, r = 5000) {
  p <- length(tPred)
  sim <- rmvnorm(r, post$mu_joint, post$cov_joint)
  f <- sim[,1:p]
  df <- sim[,(p+1):(2*p)]
  ddf <- sim[,(2*p+1):(3*p)]
  list(f = t(f), df = t(df), ddf = t(ddf))
}

TDI <- function(post) {
  p <- length(tPred)
  sapply(1:p, function(i) 1 - pnorm(0, post$mu_joint[p + i], sqrt(post$cov_joint[p + i, p + i])))
}

ETI <- function(post) {
  p <- length(tPred)
  mu_df <- post$mu_joint[(p+1):(2*p)]
  mu_ddf <- post$mu_joint[(2*p+1):(3*p)]
  sigma_df <- diag(post$cov_joint[(p+1):(2*p), (p+1):(2*p)])
  sigma_ddf <- diag(post$cov_joint[(2*p+1):(3*p), (2*p+1):(3*p)])
  sigma_df_ddf <- diag(post$cov_joint[(p+1):(2*p), (2*p+1):(3*p)])
  
  omega <- sigma_df_ddf / (sqrt(sigma_df) * sqrt(sigma_ddf))
  eta <- mu_ddf / (sigma_ddf * sqrt(1 - omega^2)) - (mu_df * omega) / (sigma_df * sqrt(1 - omega^2))
  
  ETI <- sqrt(sigma_ddf)/sqrt(sigma_df) * sqrt(1 - omega^2) * dnorm(mu_df/sqrt(sigma_df)) * (2*dnorm(eta) + eta*(2*pnorm(eta)-1))
  ETI
}

ETI.empirical <- function(y, thres=0) {
  x1 <- sum(y[1:(length(y) - 1)] > thres & y[2:length(y)] < thres)
  x2 <- sum(y[1:(length(y) - 1)] < thres & y[2:length(y)] > thres)
  x1 + x2
}

sim <- doSim(post)
teddy <- TDI(post)
eddy <- ETI(post)


nPlot <- 200
matplot(tPredCalendar, sim$f[,1:nPlot], type="l", lty = 1, col="lightgray")
lines(tPredCalendar, post$mu_f)
lines(tPredCalendar, apply(sim$f, 1, mean), col = 2)

matplot(tPredCalendar, sim$df[,1:nPlot], type="l", lty = 1, col="lightgray")
lines(tPredCalendar, post$mu_df)
lines(tPredCalendar, apply(sim$df, 1, mean), col = 2)
abline(h = 0, lty = 2)

matplot(tPredCalendar, sim$ddf[,1:nPlot], type="l", lty = 1, col="lightgray")
lines(tPredCalendar, post$mu_ddf)
lines(tPredCalendar, apply(sim$ddf, 1, mean), col = 2)

plot(tPredCalendar, apply(sim$df > 0, 1, mean), type="l")
lines(tPredCalendar, teddy, col = 2)

plot(tPredCalendar, eddy, type="l")
mean(sapply(1:ncol(sim$df), function(i) ETI.empirical(sim$df[,i])))
integrate(approxfun(tPredCalendar, eddy), 1998, 2018)



####
par(mfrow=c(2,1), bty="n")
matplot(tPredCalendar, sim$df[,1:nPlot], type="l", lty = 1, col="lightgray")
lines(tPredCalendar, post$mu_df)
lines(tPredCalendar, apply(sim$df, 1, mean), col = 2)
abline(h = 0, lty = 2)

plot(tPredCalendar, eddy, type="l")
mean(sapply(1:ncol(sim$df), function(i) ETI.empirical(sim$df[,i])))
integrate(approxfun(tPredCalendar, eddy), 1998, 2018)


plot(tPredCalendar, sim$df[, 4], type="l")
abline(h = 0, lty = 2)

library(rootSolve)
uniroot.all(approxfun(tPredCalendar, sim$df[,4]), interval = c(1998, 2018))

roots <- sapply(1:5000, function(i) uniroot.all(approxfun(tPredCalendar, sim$df[,i]), interval = c(1998, 2018)))
mean(sapply(roots, length))
quantile(sapply(roots, length), c(0.025, 0.5, 0.975))
plot(density(unlist(roots), bw="SJ"))


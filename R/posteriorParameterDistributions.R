rm(list=ls())

library(DEoptim)
library(mvtnorm)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

########################################################
# Load data
########################################################
load("../Data/smoking.RData")
dat <- data.frame(t = smoking$year, y = smoking$p)
tPred <- seq(1998, 2018, length.out = 5)

########################################################
# Estimate parameters form marginal likelihood
########################################################
rqCov <- function(s, t, alpha, rho, nu) {
  alpha^2 * (1 + (s-t)^2 / (2 * nu * rho^2))^(-nu)
}

ctl <- DEoptim.control(itermax = 1000, trace = 100)
set.seed(1234)
opt.rq <- DEoptim(function(par) {
  mu <- rep(par[1], nrow(dat))
  cMat <- outer(dat$t, dat$t, rqCov, par[2], par[3], par[4]) + diag(par[5]^2 , nrow(dat))
  -mvtnorm::dmvnorm(dat$y, mu, cMat, log=TRUE)
}, lower = c(0,0,0,0,0), upper = c(50,50,50,50,50), control = ctl)
par.rq <- opt.rq$optim$bestmem

########################################################
# Run Stan script
########################################################
sDat <- list(n = nrow(dat), t = dat$t, y = dat$y, p = length(tPred), tPred = tPred)
sDat$mu_mu <- par.rq[1]
sDat$alpha_mu <- par.rq[2]
sDat$rho_mu <- par.rq[3]
sDat$nu_mu <- par.rq[4]
sDat$sigma_mu <- par.rq[5]

iter <- 100000
seed <- 12345

m <- stan_model("gptrend.stan")
fit <- sampling(m, data = sDat, iter = iter, seed = seed)
summary(fit, c("mu", "alpha", "rho", "nu", "sigma"))$summary

par(mfrow=c(1,3), bty="n")
plot(density(log(extract(fit, "alpha")$alpha)), xlab=expression(log(alpha)), ylab="Posterior density", 
     main=expression("Posterior density of " ~ log(alpha)))
abline(v = log(par.rq[2]))

plot(density(log(extract(fit, "rho")$rho)), xlab=expression(log(rho)), ylab="Posterior density", 
     main=expression("Posterior density of " ~ log(rho)), xlim = c(0.5, 2))
abline(v = log(par.rq[3]))

plot(density(log(extract(fit, "nu")$nu)), xlab=expression(log(nu)), ylab="Posterior density",
     main=expression("Posterior density of " ~ log(nu)), xlim=c(-4,4))
abline(v = log(par.rq[4]))

quantile(extract(fit, "alpha")$alpha, c(0.025, 0.5, 0.975))
par.rq[2]

quantile(extract(fit, "rho")$rho, c(0.025, 0.5, 0.975))
par.rq[3]

quantile(extract(fit, "nu")$nu, c(0.025, 0.5, 0.975))
par.rq[4]

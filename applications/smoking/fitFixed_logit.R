rm(list=ls())

library(DEoptim)
library(mvtnorm)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

########################################################
# Load data
########################################################
logit <- function(p) log(p / (1-p))

load("smoking.RData")
dat <- data.frame(t = smoking$year, y = logit(smoking$p/100))
tPred <- seq(1998, 2018, length.out = 200) #length = 500 used for manuscript

########################################################
# Estimate Empirical Bayes parameters
########################################################
rqCov <- function(s, t, alpha, rho, nu) {
  alpha^2 * (1 + (s-t)^2 / (2 * nu * rho^2))^(-nu)
}

ctl <- DEoptim.control(itermax = 5000, trace = 100)
set.seed(1234)
opt.rq <- DEoptim(function(par) {
  mu <- rep(par[1], nrow(dat))
  cMat <- outer(dat$t, dat$t, rqCov, par[2], par[3], par[4]) + diag(par[5]^2 , nrow(dat))
  -mvtnorm::dmvnorm(dat$y, mu, cMat, log=TRUE)
}, lower = c(-50,0,0,0,0), upper = c(50,50,50,50,50), control = ctl)
par.rq <- opt.rq$optim$bestmem

########################################################
# Run Stan script
########################################################
sDat <- list(n = nrow(dat), t = dat$t, y = dat$y, p = length(tPred), tPred = tPred)
sDat$mu <- par.rq[1]
sDat$alpha <- par.rq[2]
sDat$rho <- par.rq[3]
sDat$nu <- par.rq[4]
sDat$sigma <- par.rq[5]

iter <- 10000
seed <- 12345

m <- stan_model("../../Stan/gptrendFixed.stan")
fit <- sampling(m, data = sDat, iter = iter, seed = seed, algorithm = "Fixed_param")
pred <- extract(fit, "pred")$pred

########################################################
# Plot fit
########################################################
band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

pdf("../../figures/smoking_gpFixed_logit.pdf", width = 8, height = 3)
par(mfrow=c(1,3), bty="n", mar = c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0))
plot(dat$t, dat$y, pch = 19, ylim=c(-1.4, -0.4), xaxt="n", xlab="Year", ylab="f | logit(Y)", type="n")
axis(1, 1998:2018)
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.005), apply(pred[,,1], 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.025), apply(pred[,,1], 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.25), apply(pred[,,1], 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,1], 2, mean), lwd = 2)
points(dat$t, dat$y, pch = 19, cex=0.8)
legend("topleft", c("Mean", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"), 
       lwd = 2, bty="n", cex=0.7, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5)

plot(tPred, apply(pred[,,3], 2, mean), lwd = 2, type="n", ylim = c(-0.15, 0.1), xaxt="n", xlab="Year", ylab="df | logit(Y)")
axis(1, 1998:2018)
band(tPred, apply(pred[,,3], 2, quantile, prob = 0.005), apply(pred[,,3], 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,3], 2, quantile, prob = 0.025), apply(pred[,,3], 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,3], 2, quantile, prob = 0.25), apply(pred[,,3], 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,3], 2, mean), lwd = 2)
abline(h = 0, lty = 2)
legend("topleft", c("Mean", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"), 
       lwd = 2, bty="n", cex=0.7, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5)

plot(tPred, t(pred[1,,5])*100, type="l", lty = 1, lwd = 2, xaxt="n", xlab="Year", 
     ylab="Trend Direction Index [%]", ylim=c(0,100))
axis(1, 1998:2018)
abline(h = 50, lty = 2)
title("TDI(2018, Year - 2018)", font.main=1)
dev.off()

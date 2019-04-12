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
tPred <- seq(1998, 2018, length.out = 200) #length = 500 used for manuscript

########################################################
# Estimate Empirical Bayes parameters
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
par.rq
signif(par.rq, 3)

########################################################
# Run Stan script
########################################################
sDat <- list(n = nrow(dat), t = dat$t, y = dat$y, p = length(tPred), tPred = tPred)
sDat$mu <- par.rq[1]
sDat$alpha <- par.rq[2]
sDat$rho <- par.rq[3]
sDat$nu <- par.rq[4]
sDat$sigma <- par.rq[5]

iter <- 10000 #25000 used for manuscript
seed <- 12345

m <- stan_model("gptrendFixed.stan")
fit <- sampling(m, data = sDat, iter = iter, seed = seed, algorithm = "Fixed_param")
pred <- extract(fit, "pred")$pred
save(tPred, dat, pred, file="predFixed.RData")

########################################################
# Get some summary statistics
########################################################
load("predFixed.RData")

#TDI summary at selected years
round(approxfun(tPred, pred[1,,5]*100)(2013:2018), 2)

uniroot(function(x) approxfun(tPred, pred[1,,5]*100)(x) - 50, c(2010, 2018))$root

#ETI summary
round(integrate(splinefun(tPred, pred[1,,6]), 1998, 2018)$value, 2)
round(integrate(splinefun(tPred, pred[1,,6]), 2008, 2018)$value, 2)

########################################################
# Plot fit
########################################################
band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

pdf("../figures/gpFit_fixed.pdf", width = 8, height = 6)
par(mfrow=c(2,2), bty="n", mar = c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0))
plot(dat$t, dat$y, pch = 19, xaxt="n", ylim=c(20,40), xlab="Year", ylab="f | Y", type="n")
axis(1, 1998:2018)
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.005), apply(pred[,,1], 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.025), apply(pred[,,1], 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.25), apply(pred[,,1], 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,1], 2, mean), lwd = 2)
points(dat$t, dat$y, pch = 19, cex=0.8)
legend("topleft", c("Mean", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"), 
       lwd = 2, bty="n", cex=0.7, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5)

plot(tPred, apply(pred[,,3], 2, mean), lwd = 2, type="n", ylim = c(-3, 3), xaxt="n", xlab="Year", ylab="df | Y")
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

plot(tPred, t(pred[1,,6]), type="l", lty = 1, lwd = 2, xaxt="n", xlab="Year", 
     ylab="Expected Trend Instability", ylim=c(0,1.8), yaxt="n")
title(bquote(ETI[Year] ~ "([1998; 2018])"), font.main=1)
axis(1, 1998:2018)
axis(2, seq(0, 1.8, length.out=7))
dev.off()

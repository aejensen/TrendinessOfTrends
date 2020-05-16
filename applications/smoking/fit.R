rm(list=ls())

library(DEoptim)
library(mvtnorm)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

########################################################
# Load data
########################################################
load("smoking.RData")
dat <- data.frame(t = smoking$year, y = smoking$p)
tPred <- seq(1998, 2018, length.out = 200) #length = 500 used for manuscript

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
par.rq

########################################################
# Run Stan script
########################################################
sDat <- list(n = nrow(dat), t = dat$t, y = dat$y, p = length(tPred), tPred = tPred)
sDat$mu_mu <- par.rq[1]
sDat$alpha_mu <- par.rq[2]
sDat$rho_mu <- par.rq[3]
sDat$nu_mu <- par.rq[4]
sDat$sigma_mu <- par.rq[5]

iter <- 10000 #25000 used for manuscript
seed <- 12345

m <- stan_model("../../Stan/gptrend.stan")
fit <- sampling(m, data = sDat, iter = iter, seed = seed)
summary(fit, c("mu", "alpha", "rho", "nu", "sigma"))$summary
pred <- extract(fit, "pred")$pred
save(tPred, dat, pred, file="pred.RData")

########################################################
# Trace plots
########################################################
pdf("../../figures/smoking_traceplot.pdf", width = 8, height = 6)
traceplot(fit, pars=c("mu", "alpha", "rho", "nu", "sigma")) + theme(legend.position = "top")
dev.off()

#color_scheme_set("viridis")
#mcmc_trace(posterior, pars = "mu")

########################################################
# Get some summary statistics
########################################################
load("pred.RData")

#TDI summary
round(cbind(2018-(0:5), 
            approxfun(tPred, apply(pred[,,5]*100, 2, median))(2018-(0:5)),
            approxfun(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.025))(2018-(0:5)),
            approxfun(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.975))(2018-(0:5))), 2)

uniroot(function(x) approxfun(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.025))(x) - 50, c(2010, 2018))$root
uniroot(function(x) approxfun(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.5))(x) - 50, c(2010, 2018))$root
uniroot(function(x) approxfun(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.975))(x) - 50, c(2010, 2018))$root


#ETI summary
ETI1 <- sapply(1:dim(pred)[1], function(i) integrate(splinefun(tPred, pred[i,,6]), 1998, 2018)$value)
ETI2 <- sapply(1:dim(pred)[1], function(i) integrate(splinefun(tPred, pred[i,,6]), 2008, 2018)$value)
round(quantile(ETI1, c(0.025, 0.5, 0.975)), 2)
round(quantile(ETI2, c(0.025, 0.5, 0.975)), 2)

########################################################
# Plot 1
########################################################
band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

pdf("../figures/gpFit.pdf", width = 8, height = 6)
par(mfrow=c(2,2), bty="n", mar =  c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0))
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

plot(tPred, apply(pred[,,5]*100, 2, median), type="n", lty = 1, lwd = 2, xaxt="n", xlab="Year", ylab="TDI [%]", ylim=c(0,100))
band(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.005), apply(pred[,,5]*100, 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.025), apply(pred[,,5]*100, 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.25), apply(pred[,,5]*100, 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,5]*100, 2, median), lwd = 2)
axis(1, 1998:2018)
abline(h = 50, lty = 2)
legend("topleft", c("Median", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"),
       lwd = 2, bty="n", cex=0.7, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5)

plot(tPred, t(pred[1,,6]), type="n", lty = 1, lwd = 2, xaxt="n", xlab="Year", ylab="Local ETI", ylim=c(0,1.8), yaxt="n")
axis(1, 1998:2018)
axis(2, seq(0, 1.8, length.out=7))
band(tPred, apply(pred[,,6], 2, quantile, prob = 0.005), apply(pred[,,6], 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,6], 2, quantile, prob = 0.025), apply(pred[,,6], 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,6], 2, quantile, prob = 0.25), apply(pred[,,6], 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,6], 2, mean), lwd = 2)
legend("topleft", c("Median", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"), 
       lwd = 2, bty="n", cex=0.7, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5)
dev.off()

########################################################
# Plot 2
########################################################
ETI1_d <- density(ETI1)
ETI1_d_f <- approxfun(ETI1_d$x, ETI1_d$y)
ETI1_seq1 <- seq(quantile(ETI1, 0.025), quantile(ETI1, 0.975), length.out=100)
ETI1_seq2 <- seq(quantile(ETI1, 0.25), quantile(ETI1, 0.75), length.out=100)

ETI2_d <- density(ETI2)
ETI2_d_f <- approxfun(ETI2_d$x, ETI2_d$y)
ETI2_seq1 <- seq(quantile(ETI2, 0.025), quantile(ETI2, 0.975), length.out=100)
ETI2_seq2 <- seq(quantile(ETI2, 0.25), quantile(ETI2, 0.75), length.out=100)

pdf("../figures/gpFit_ETI.pdf", width = 8, height = 4)
par(mfrow=c(1,2), bty="n", mar =  c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0))
plot(0, 0, type="n", xlim=c(0,6), ylim=c(0,0.8), xlab="ETI", ylab="Posterior density", xaxt="n", yaxt="n")
axis(1, seq(0, 6, 1), cex.axis = 0.8)
axis(2, cex.axis = 0.8)
band(ETI1_seq1, rep(0, length(ETI1_seq1)), ETI1_d_f(ETI1_seq1), col="gray80")
band(ETI1_seq2, rep(0, length(ETI1_seq2)), ETI1_d_f(ETI1_seq2), col="gray65")
lines(rep(quantile(ETI1, 0.5), 2), c(0, ETI1_d_f(quantile(ETI1, 0.5))), lwd = 1, lty=2)
curve(ETI1_d_f(x), 0, 6, ylim=c(0,1), add=TRUE, lwd = 2)
title("ETI([1998; 2018])", font.main=1)
legend("topleft", c("Median", "50%", "95%"), col = c("black", "gray65", "gray80"), 
       lwd = 2, bty="n", cex=0.9, lty=c(2,NA, NA), pch = c(NA, 15, 15), pt.cex=1.5)

plot(0, 0, type="n", xlim=c(0,3), ylim=c(0,2.5), xlab="ETI", ylab="Posterior density", xaxt="n", yaxt="n")
axis(1, seq(0, 3, 0.5), cex.axis = 0.8)
axis(2, cex.axis = 0.8)
band(ETI2_seq1, rep(0, length(ETI2_seq1)), ETI2_d_f(ETI2_seq1), col="gray80")
band(ETI2_seq2, rep(0, length(ETI2_seq2)), ETI2_d_f(ETI2_seq2), col="gray65")
lines(rep(quantile(ETI2, 0.5), 2), c(0, ETI2_d_f(quantile(ETI2, 0.5))), lwd = 1, lty=2)
curve(ETI2_d_f(x), 0, 3, ylim=c(0,2), add=TRUE, lwd = 2)
title("ETI([2008; 2018])", font.main=1)
legend("topleft", c("Median", "50%", "95%"), col = c("black", "gray65", "gray80"), 
       lwd = 2, bty="n", cex=0.9, lty=c(2,NA, NA), pch = c(NA, 15, 15), pt.cex=1.5)
dev.off()


rm(list=ls())

library(DEoptim)
library(mvtnorm)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

########################################################
# Load data
########################################################
#Data downloaded May 23th from 
#https://github.com/pcm-dpc/COVID-19/tree/master/dati-andamento-nazionale
d <- read.csv("dpc-covid19-ita-andamento-nazionale.csv")

#We look at the number of new positives
dat <- data.frame(t = 0:(nrow(d)-1), y = d$nuovi_positivi)
dat$yScaled <- dat$y / max(dat$y)
tPred <- seq(0, nrow(d)-1, length.out = 200)

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
  -mvtnorm::dmvnorm(dat$yScaled, mu, cMat, log=TRUE)
}, lower = c(0,0,0,0,0), upper = c(1,10,50,50,50), control = ctl)
par.rq <- opt.rq$optim$bestmem
par.rq

########################################################
# Run Stan script
########################################################
sDat <- list(n = nrow(dat), t = dat$t, y = dat$yScaled, p = length(tPred), tPred = tPred)
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

save.image("~/Desktop/covid19Trendiness.RData")

########################################################
# Plot fit
########################################################
load("~/Desktop/covid19Trendiness.RData")

band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

pdf("../../figures/covid19_gpFit.pdf", width = 8, height = 3)
par(mfrow=c(1,3), bty="n", mar = c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0))
plot(dat$t, dat$y, pch = 19, xlab="Number of days since Feb 24th", ylab="f | Y", 
     type="n", ylim=c(0, 8000), xaxt="n", xlim=c(0, 90))
axis(1, seq(0, 90, 10))
band(tPred, apply(pred[,,2]*max(dat$y), 2, quantile, prob = 0.025), 
     apply(pred[,,2]*max(dat$y), 2, quantile, prob = 0.975), col = "gray85")
band(tPred, apply(pred[,,1]*max(dat$y), 2, quantile, prob = 0.025), 
     apply(pred[,,1]*max(dat$y), 2, quantile, prob = 0.975), col = "gray65")
lines(tPred, apply(pred[,,1]*max(dat$y), 2, mean), lwd = 2)
points(dat$t, dat$y, pch = 19, cex=0.6)
legend("topleft", 
       c("Mean", "95% credible interval", "95% posterior prediction interval"), 
       col = c("black", "gray65", "gray85"), lwd = 2, bty="n", cex=0.7, 
       lty = c(1, NA, NA), pch = c(NA, 15, 15), pt.cex=1.5)

plot(tPred, apply(pred[,,3]*max(dat$y), 2, mean), lwd = 2, type="n", yaxt="n", xlim=c(0,90),
     ylim=c(-500, 500), xlab="Number of days since Feb 24th", ylab="df | Y", xaxt="n")
axis(2, seq(-500, 500, length.out=5))
axis(1, seq(0, 90, 10))
band(tPred, apply(pred[,,3]*max(dat$y), 2, quantile, prob = 0.025), 
     apply(pred[,,3]*max(dat$y), 2, quantile, prob = 0.975), col = "gray65")
lines(tPred, apply(pred[,,3]*max(dat$y), 2, mean), lwd = 2)
abline(h = 0, lty = 2)
legend("topleft", 
       c("Mean", "95% credible interval"), col = c("black", "gray65"), 
       lwd = 2, bty="n", cex=0.7, lty = c(1, NA), pch = c(NA, 15), pt.cex=1.5)

plot(tPred, t(pred[1,,5])*100, type="l", lty = 1, lwd = 2, xlab="Number of days since Feb 24th", 
     ylab="Trend Direction Index [%]", ylim=c(0,100), xaxt="n", xlim=c(0,90))
axis(1, seq(0, 90, 10))
abline(h = 50, lty = 2)
dev.off()

#Summaries
max(pred[1,,5]*100)
rootSolve::uniroot.all(function(x) approxfun(tPred, t(pred[1,,5])*100)(x) - 95, c(0, 90))
rootSolve::uniroot.all(function(x) approxfun(tPred, t(pred[1,,5])*100)(x) - 50, c(0, 90))

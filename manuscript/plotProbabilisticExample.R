rm(list=ls())
library(mvtnorm)
library(fields)

add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

seCov <- function(s, t, alpha, l) {
  alpha^2 * exp(-(s-t)^2 / (2*l^2))
}

seCovD1 <- function(s, t, alpha, l) {
  seCov(s, t, alpha, l) * (t-s) / l^2
}

seCovD2 <- function(s, t, alpha, l) {
  seCov(s, t, alpha, l) * (s-t) / l^2
}

seCovD12 <- function(s, t, alpha, l) {
  seCov(s, t, alpha, l) * (l^(-2) - l^(-4) * (s-t)^2)
}

getPost <- function(tObs, yObs, tPred, alpha, l, sigmaResid = 1e-9) {
  k1 <- outer(tPred, tObs, seCov, alpha = alpha, l = l)
  k2 <- outer(tObs, tObs, seCov, alpha = alpha, l = l) + diag(sigmaResid^2, length(tObs))
  k3 <- outer(tPred, tPred, seCov, alpha = alpha, l = l)
  k4 <- outer(tObs, tPred, seCov, alpha = alpha, l = l)
  
  mu.f <- as.numeric(k1 %*% solve(k2) %*% yObs)
  cov.f <- k3 - k1 %*% solve(k2) %*% k4
  sim.f <- t(mvtnorm::rmvnorm(k, mu.f, cov.f))
  list(sim = sim.f, mu = mu.f)
}

getPostD <- function(tObs, yObs, tPred, alpha, l, sigmaResid = 1e-9, k = k) {
  k1 <- outer(tPred, tObs, seCovD1, alpha = alpha, l = l)
  k2 <- outer(tObs, tObs, seCov, alpha = alpha, l = l) + diag(sigmaResid^2, length(tObs))
  k3 <- outer(tPred, tPred, seCovD12, alpha = alpha, l = l)
  k4 <- outer(tObs, tPred, seCovD2, alpha = alpha, l = l)
  
  mu.f <- as.numeric(k1 %*% solve(k2) %*% yObs)
  cov.f <- k3 - k1 %*% solve(k2) %*% k4
  sim.f <- t(mvtnorm::rmvnorm(k, mu.f, cov.f))
  list(sim = sim.f, mu = mu.f)
}


n <- 100
k <- 150
tPred <- seq(0, 1, length.out = n)

doPlotF <- function(obsT, obsY, ...) {
  set.seed(12345)
  matplot(tPred, getPost(obsT, obsY, tPred, alpha = 2, l = 0.1)$sim, type="l", lty=1, ylim = c(-7, 7),
          xlab = "t", ylab="f(t)", col = add.alpha(1:6, 0.3), xaxt="n", yaxt="n", ...)
  lines(tPred, getPost(obsT, obsY, tPred, alpha = 2, l = 0.1)$mu, lwd = 3)
  points(obsT, obsY, pch=19, cex=1.5)
  axis(1, cex.axis = 0.69)
  axis(2, cex.axis = 0.69)
}

doPlotDF <- function(obsT, obsY) {
  set.seed(12345)
  matplot(tPred, getPostD(obsT, obsY, tPred, alpha = 2, l = 0.1, k = 150)$sim, type="l", lty=1, ylim = c(-100, 100),
          xlab = "t", ylab="df(t)", col = add.alpha(1:6, 0.3), xaxt="n", yaxt="n")
  lines(tPred, getPostD(obsT, obsY, tPred, alpha = 2, l = 0.1, k = 150)$mu, lwd = 3)
  lines(c(0,1), c(0, 0), lty=2)
  axis(1, cex.axis = 0.69)
  axis(2, cex.axis = 0.69)
}

doPlotTCI <- function(obsT, obsY) {
  set.seed(12345)
  dfSim <- getPostD(obsT, obsY, tPred, alpha = 2, l = 0.1, k = 10^6)$sim
  prob <- apply(dfSim > 0, 1, mean)
  plot(tPred, prob, type="l", ylim=c(0,1), xaxt="n", yaxt="n", lwd=2, xlab="t", ylab="TDI")
  lines(c(0,1), c(0.5, 0.5), lty=2)
  axis(1, cex.axis = 0.69)
  axis(2, cex.axis = 0.69)
}

pdf("probabilisticExample.pdf", width = 8, height = 5)
par(mfrow=c(3,3), bty="n", mar =  c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0))
doPlotF(c(0.1), c(0), main = "One observation", font.main = 1)
doPlotF(c(0.1, 0.3), c(0, 2), main = "Two observations", font.main = 1)
doPlotF(c(0.1, 0.3, 0.5, 0.9), c(0, 2, -3, 2),  main = "Four observations", font.main = 1)

doPlotDF(c(0.1), c(0))
doPlotDF(c(0.1, 0.3), c(0, 2))
doPlotDF(c(0.1, 0.3, 0.5, 0.9), c(0, 2, -3, 2))

doPlotTCI(c(0.1), c(0))
doPlotTCI(c(0.1, 0.3), c(0, 2))
doPlotTCI(c(0.1, 0.3, 0.5, 0.9), c(0, 2, -3, 2))
dev.off()


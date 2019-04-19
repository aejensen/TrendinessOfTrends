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


n <- 251
k <- 150
tPred <- seq(0, 1, length.out = n)

doPlotF <- function(obsT, obsY, ...) {
    beforeT <- which(tPred <= max(obsT))
    afterT <- which(tPred > max(obsT))
  
  set.seed(12345)
  postDat <- getPost(obsT, obsY, tPred, alpha = 1.5, l = 0.1)
  
  matplot(tPred, postDat$sim, type="n", lty=1, ylim = c(-6, 6),
          xlab = "t", ylab="f(t)", col = add.alpha(1:6, 0.3), xaxt="n", yaxt="n", ...)
  
  matplot(tPred[beforeT], postDat$sim[beforeT,], type="l", lty=1, col = add.alpha(1:6, 0.6), add=TRUE)  
  if(length(afterT) > 0) {
    matplot(tPred[afterT], postDat$sim[afterT,], type="l", lty=1, col = add.alpha(1:6, 0.2), add=TRUE)  
  }
  
  lines(tPred[beforeT], postDat$mu[beforeT], lwd = 3)
  if(length(afterT) > 0) {
    lines(tPred[afterT], postDat$mu[afterT], lwd = 3, col = add.alpha("black", 0.5))
  }
  
  points(obsT, obsY, pch=19, cex=1.5)
  axis(1, cex.axis = 0.9)
  axis(2, cex.axis = 0.9)
  
  if(length(afterT) > 0) {
    lines(rep(obsT[length(obsT)], 2), c(-6, 6), lty = 3)
  }
}

doPlotDF <- function(obsT, obsY) {
  beforeT <- which(tPred <= max(obsT))
  afterT <- which(tPred > max(obsT))
  
  set.seed(12345)
  postDat <- getPostD(obsT, obsY, tPred, alpha = 1.5, l = 0.1, k = 150)
  
  matplot(tPred, postDat$sim, type="n", lty=1, ylim = c(-60, 60),
          xlab = "t", ylab="df(t)", col = add.alpha(1:6, 0.3), xaxt="n", yaxt="n")

  matplot(tPred[beforeT], postDat$sim[beforeT,], type="l", lty=1, col = add.alpha(1:6, 0.6), add=TRUE)  
  if(length(afterT) > 0) {
    matplot(tPred[afterT], postDat$sim[afterT,], type="l", lty=1, col = add.alpha(1:6, 0.2), add=TRUE)  
  }
  
  lines(tPred[beforeT], postDat$mu[beforeT], lwd = 3)
  if(length(afterT) > 0) {
    lines(tPred[afterT], postDat$mu[afterT], lwd = 3, col = add.alpha("black", 0.5))
  }
  
  lines(c(0,1), c(0, 0), lty=3)
  axis(1, cex.axis = 0.9)
  axis(2, seq(-60, 60, length.out=5), cex.axis = 0.9)
  if(length(afterT) > 0) {
    lines(rep(obsT[length(obsT)], 2), c(-60, 60), lty = 3)
  }
}

doPlotTCI <- function(obsT, obsY) {
  beforeT <- which(tPred <= max(obsT))
  afterT <- which(tPred > max(obsT))
  
  set.seed(12345)
  dfSim <- getPostD(obsT, obsY, tPred, alpha = 1.5, l = 0.1, k = 500000)$sim
  prob <- apply(dfSim > 0, 1, mean)*100
  
  plot(tPred, prob, type="n", ylim=c(0,100), xaxt="n", yaxt="n", lwd=2, xlab="t", ylab="Trend Direction Index [%]")

  lines(tPred[beforeT], prob[beforeT], lwd = 2)
  if(length(afterT) > 0) {
    lines(tPred[afterT], prob[afterT], lwd = 2, col = add.alpha("black", 0.6))
  }
  
  lines(c(0,1), c(50, 50), lty=3)
  axis(1, cex.axis = 0.9)
  axis(2, cex.axis = 0.9)
  if(length(afterT) > 0) {
    lines(rep(obsT[length(obsT)], 2), c(0, 100), lty = 3)
  }
}

pdf("../figures/TDIexample.pdf", width = 8, height = 5)
par(mfrow=c(3,3), bty="n", mar =  c(2.3, 2.3, 0.8, 0), mgp=c(1.3,0.4,0))
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


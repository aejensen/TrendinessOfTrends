rm(list=ls())
library(mvtnorm)
library(fields)

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

doPlotPrior <- function(tPred) {
  cMat11 <- outer(tPred, tPred, seCov, alpha = 2, l = 0.1)
  cMat12 <- outer(tPred, tPred, seCovD1, alpha = 2, l = 0.1)
  cMat21 <- outer(tPred, tPred, seCovD2, alpha = 2, l = 0.1)
  cMat22 <- outer(tPred, tPred, seCovD12, alpha = 2, l = 0.1)
  cMat <- rbind(cbind(cMat11, cMat12), cbind(cMat21, cMat22))
  
  set.seed(12345)
  #sim <- t(rmvnorm(k, rep(0, 2*length(tPred)), cMat))
  sim.f <- t(rmvnorm(k, rep(0, length(tPred)), cMat11))

  #sim.f <- sim[1:n,]
  #sim.df <- sim[(n+1):(2*n),]

  par(mar =  c(5.1, 4.1, 4.1, 2.1), bty="n", bty="n", bg=NA, col.lab="white", col.axis="white")
  matplot(tPred, sim.f, type = "l", lty=1, ylim = c(-7, 7),
          xlab = "t", ylab="f(t)", col = 1:6, xaxt="n", yaxt="n")
  axis(1, col = "white")
  axis(2, col = "white")
}

doPlot <- function(obsT, obsY) {
  set.seed(12345)
  par(mfrow=c(1,2), bty="n", bg=NA, col.lab="white", col.axis="white", mar =  c(5.1, 4.1, 4.1, 0))
  matplot(tPred, getPost(obsT, obsY, tPred, alpha = 2, l = 0.1)$sim, type="l", lty=1, ylim = c(-8, 8),
          xlab = "t", ylab="f(t)", col = 1:6, xaxt="n", yaxt="n")
  lines(tPred, getPost(obsT, obsY, tPred, alpha = 2, l = 0.1)$mu, lwd = 4, col = "white")
  points(obsT, obsY, pch=19, cex=1.5, col = "white")
  axis(1, col = "white")
  axis(2, seq(-8, 8, 4), col = "white")
  lines(rep(obsT[length(obsT)], 2), c(-8, 8), lty = 2, col = "white", lwd = 2)
  
  matplot(tPred, getPostD(obsT, obsY, tPred, alpha = 2, l = 0.1, k = 150)$sim, type="l", lty=1, ylim = c(-100, 100),
          xlab = "t", ylab="df(t)", col = 1:6, xaxt="n", yaxt="n")
  lines(tPred, getPostD(obsT, obsY, tPred, alpha = 2, l = 0.1, k = 150)$mu, lwd = 4, col = "white")
  abline(h = 0, lty=2, col = "white")
  axis(1, col = "white")
  axis(2, col = "white")
  lines(rep(obsT[length(obsT)], 2), c(-100, 100), lty = 2, col = "white", lwd = 2)
}

doPlot2 <- function(obsT, obsY) {
  set.seed(12345)
  par(mfrow=c(1,2), bty="n", bg=NA, col.lab="white", col.axis="white", mar =  c(5.1, 4.1, 4.1, 0))
  matplot(tPred, getPost(obsT, obsY, tPred, alpha = 2, l = 0.1)$sim, type="l", lty=1, ylim = c(-8, 8),
          xlab = "t", ylab="f(t)", col = 1:6, xaxt="n", yaxt="n")
  lines(tPred, getPost(obsT, obsY, tPred, alpha = 2, l = 0.1)$mu, lwd = 4, col = "white")
  points(obsT, obsY, pch=19, cex=1.5, col = "white")
  axis(1, col = "white")
  axis(2, seq(-8, 8, 4), col = "white")
  lines(rep(obsT[length(obsT)], 2), c(-8, 8), lty = 2, col = "white", lwd = 2)
  
  dfSim <- getPostD(obsT, obsY, tPred, alpha = 2, l = 0.1, k = 10^6)$sim 
  prob <- apply(dfSim > 0, 1, mean)
  plot(tPred, prob, type="l", col="white", ylim=c(0,1), xaxt="n", yaxt="n", lwd=4,
       xlab="t", ylab="Trend Direction Index")
  abline(h = 0.5, lty=2, col = "white")
  #abline(h = 0, lty=2, col = "white")
  axis(1, col = "white")
  axis(2, col = "white")
  lines(rep(obsT[length(obsT)], 2), c(0, 1), lty = 2, col = "white", lwd = 2)
}

pdf(file = "postAni%02d.pdf", width = 8, height = 6, onefile = FALSE)
doPlotPrior(tPred)
doPlot(c(0.1), c(0))
doPlot(c(0.1, 0.3), c(0, 2))
doPlot(c(0.1, 0.3, 0.35), c(0, 2, 0))
doPlot(c(0.1, 0.3, 0.35, 0.5), c(0, 2, 0, -3))
doPlot(c(0.1, 0.3, 0.35, 0.5, 0.8), c(0, 2, 0, -3, 2))
doPlot(c(0.1, 0.3, 0.35, 0.5, 0.8, 0.9), c(0, 2, 0, -3, 2, 0))
#doPlot(c(0.1, 0.3, 0.35, 0.5, 0.8, 0.9, 1), c(0, 2, 0, -3, 2, 0, 0))
dev.off()

pdf(file = "probAni%02d.pdf", width = 8, height = 6, onefile = FALSE)
doPlot2(c(0.1), c(0))
doPlot2(c(0.1, 0.3), c(0, 2))
doPlot2(c(0.1, 0.3, 0.35), c(0, 2, 0))
doPlot2(c(0.1, 0.3, 0.35, 0.5), c(0, 2, 0, -3))
doPlot2(c(0.1, 0.3, 0.35, 0.5, 0.8), c(0, 2, 0, -3, 2))
doPlot2(c(0.1, 0.3, 0.35, 0.5, 0.8, 0.9), c(0, 2, 0, -3, 2, 0))
#doPlot2(c(0.1, 0.3, 0.35, 0.5, 0.8, 0.9, 1), c(0, 2, 0, -3, 2, 0, 0))
dev.off()


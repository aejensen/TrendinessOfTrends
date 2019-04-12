rm(list=ls())
library(fields)
library(mvtnorm)

seCov <- function(s, t, alpha, rho) {
  alpha^2 * exp(-(s-t)^2 / (2 * rho^2))
}

seCov_d2 <- function(s, t, alpha, rho) {
  seCov(s, t, alpha, rho) * (s-t) / rho^2
} 

seCov_d2_d2 <- function(s, t, alpha, rho) {
  seCov(s, t, alpha, rho) * ((s-t)^2 - rho^2) / rho^4
} 

seCov_d1_d2 <- function(s, t, alpha, rho) {
  seCov(s, t, alpha, rho) * (-(s-t)^2 + rho^2) / rho^4
} 

seCov_d1_d2_d2 <- function(s, t, alpha, rho) {
  seCov(s, t, alpha, rho) * (-(s-t)^3 + 3*(s-t)*rho^2) / rho^6
} 

seCov_d1_d1_d2_d2 <- function(s, t, alpha, rho) {
  seCov(s, t, alpha, rho) * ((s-t)^4 - 6*(s-t)^2 * rho^2 + 3*rho^4) / rho^8
} 

simData <- function(t, alpha, rho, k = 100) {
  n <- length(t)
  
  c11 <- outer(t, t, seCov, alpha, rho)
  c12 <- outer(t, t, seCov_d2, alpha, rho)
  c13 <- outer(t, t, seCov_d2_d2, alpha, rho)
  c21 <- t(c12)
  c22 <- outer(t, t, seCov_d1_d2, alpha, rho)
  c23 <- outer(t, t, seCov_d1_d2_d2, alpha, rho)
  c31 <- t(c13)
  c32 <- t(c23)
  c33 <- outer(t, t, seCov_d1_d1_d2_d2, alpha, rho)
  cMat <- rbind(cbind(c11, c12, c13), cbind(c21, c22, c23), cbind(c31, c32, c33))
  
  sim <- mvtnorm::rmvnorm(k, rep(0, 3*n), cMat)
  f <- t(sim[,1:n])
  df <- t(sim[,(n+1):(2*n)])
  ddf <- t(sim[,(2*n+1):(3*n)])
  
  list(f = f, df = df, ddf = ddf, cov = cMat)
}

empiricalEstimate <- function(y, thres=0) {
  x1 <- sum(y[1:(length(y) - 1)] > thres & y[2:length(y)] < thres)
  x2 <- sum(y[1:(length(y) - 1)] < thres & y[2:length(y)] > thres)
  x1 + x2
}

analyticEstimate <- function(t, alpha, rho) {
  c1 <- seCov_d1_d2(t, t, alpha, rho)
  c2 <- seCov_d1_d1_d2_d2(t, t, alpha, rho)
  c2 / (sqrt(c1 * c2) * pi)
}

parameterEstimate <- function(rho) {
  sqrt(3) / (pi * rho)
}

getRho <- function(D) {
  sqrt(3) / (D * pi)
}

########################################################################
n <- 100
tEval <- seq(0, 1, length.out = n)
alpha <-  0.6
rho1 <- getRho(0.25)
rho2 <- getRho(0.5)
rho3 <- getRho(1)

set.seed(1234)
d1 <- simData(tEval, alpha, rho1, k = 25)
d2 <- simData(tEval, alpha, rho2, k = 25)
d3 <- simData(tEval, alpha, rho3, k = 25)

which1 <- sapply(1:ncol(d1$df), function(k) !all(d1$df[,k] >= 0) & !all(d1$df[,k] <= 0))
which2 <- sapply(1:ncol(d2$df), function(k) !all(d2$df[,k] >= 0) & !all(d2$df[,k] <= 0))
which3 <- sapply(1:ncol(d3$df), function(k) !all(d3$df[,k] >= 0) & !all(d3$df[,k] <= 0))

col <- c("navy", "darkgoldenrod")

pdf("../figures/ETIexample.pdf", width = 8, height = 4.5)
par(mfrow=c(2,3), bty="n", mar =  c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0), font.main = 1)
matplot(tEval, d1$f, type="l", ylim=c(-2,2), col = col[which1 + 1], lty=c(1,5)[which1 + 1], 
        xlab="t", ylab="f(t)", main = "ETI([0;1]) = 0.25")
legend("topleft", c("Stable", "Unstable"), col = col, bty="n", lty=c(1,5), lwd = 2, seg.len=3)
matplot(tEval, d2$f, type="l",  ylim=c(-2,2), col = col[which2 + 1], lty=c(1,5)[which2 + 1], 
        xlab="t", ylab="f(t)", main = "ETI([0;1]) = 0.5")
matplot(tEval, d3$f, type="l", ylim=c(-2,2), col = col[which3 + 1], lty=c(1, 5)[which3 + 1], 
        xlab="t", ylab="f(t)", main = "ETI([0;1]) = 1")

matplot(tEval, d1$df, type="l", ylim=c(-1,1), col = col[which1 + 1], lty=c(1,5)[which1 + 1], 
        xlab="t", ylab="df(t)")
matplot(tEval, d2$df, type="l", ylim=c(-2,2), col = col[which2 + 1], lty=c(1,5)[which2 + 1], 
        xlab="t", ylab="df(t)")
matplot(tEval, d3$df, type="l",  ylim=c(-4,4), col = col[which3 + 1], lty=c(1,5)[which3 + 1], 
        xlab="t", ylab="df(t)")
dev.off()

### Sanity checks
c(analyticEstimate(0, alpha, rho1), parameterEstimate(rho1))
c(analyticEstimate(0, alpha, rho2), parameterEstimate(rho2))
c(analyticEstimate(0, alpha, rho3), parameterEstimate(rho3))




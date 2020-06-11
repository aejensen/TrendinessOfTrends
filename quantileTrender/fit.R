rm(list=ls())

library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

tObs <- seq(0, 1, length.out=50)
tPred <- seq(0, 1, length.out=100)

set.seed(12345)
mu <- sin(2*pi*tObs)
sigma <- 0.5*(tObs^2 + 0.1*tObs)
#plot(tObs, sigma)
y <- mu + rnorm(length(tObs), rep(0, length(tObs)), sigma)

m <- rstan::stan_model("gpQuantileTrender.stan")
iter <- 5000
seed <- 12345

quants <- c(0.1, 0.25, 0.5, 0.75, 0.9)
fits <- mclapply(quants, function(q) {
  sDat <- list(n = length(tObs), t = tObs, y = y,
               p = length(tPred), tPred = tPred,
               quantile = q)
  rstan::sampling(m, data = sDat, iter = iter, seed = seed, chains=4)  
}, mc.cores=length(quants))

save.image("quantileTrenderTest.RData")

############################################################

postMeans <- sapply(fits, function(fit) {
  apply(rstan::extract(fit, "fPred")$fPred, 2, mean) 
})

postLower <- sapply(fits, function(fit) {
  apply(rstan::extract(fit, "fPred")$fPred, 2, quantile, prob=0.025) 
})

postUpper <- sapply(fits, function(fit) {
  apply(rstan::extract(fit, "fPred")$fPred, 2, quantile, prob=0.975) 
})


col <- fields::tim.colors(5)

plot(tObs, y, pch=19, cex=0.8, bty="n", ylim=c(-2,2), xlab="t", ylab="y(t)")
lines(tPred, postMeans[,1], col=col[1], lwd=2)
lines(tPred, postLower[,1], col=col[1], lty=2, lwd=0.7)
lines(tPred, postUpper[,1], col=col[1], lty=2, lwd=0.7)

lines(tPred, postMeans[,2], col=col[2], lwd=2)
lines(tPred, postLower[,2], col=col[2], lty=2, lwd=0.7)
lines(tPred, postUpper[,2], col=col[2], lty=2, lwd=0.7)

lines(tPred, postMeans[,3], col=col[3], lwd=2)
lines(tPred, postLower[,3], col=col[3], lty=2, lwd=0.7)
lines(tPred, postUpper[,3], col=col[3], lty=2, lwd=0.7)

lines(tPred, postMeans[,4], col=col[4], lwd=2)
lines(tPred, postLower[,4], col=col[4], lty=2, lwd=0.7)
lines(tPred, postUpper[,4], col=col[4], lty=2, lwd=0.7)

lines(tPred, postMeans[,5], col=col[5], lwd=2)
lines(tPred, postLower[,5], col=col[5], lty=2, lwd=0.7)
lines(tPred, postUpper[,5], col=col[5], lty=2, lwd=0.7)

legend("topleft", paste(quants*100, "%", sep=""), col = col, lwd=2, bty="n", cex=0.7)
dev.off()


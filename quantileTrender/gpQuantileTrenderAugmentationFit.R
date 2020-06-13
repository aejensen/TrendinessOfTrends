rm(list=ls())

library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

tObs <- seq(0, 1, length.out=25)

set.seed(12345)
mu <- sin(2*pi*tObs)
sigma <- 0.5*(tObs^2 + 0.1*tObs)
#plot(tObs, sigma)
y <- mu + rnorm(length(tObs), rep(0, length(tObs)), sigma)

m <- rstan::stan_model("gpQuantileTrenderAugmentation.stan")

iter <- 10000
seed <- 12345

quants <- c(0.1, 0.25, 0.5, 0.75, 0.9)
fits <- mclapply(quants, function(q) {
  sDat <- list(n = length(tObs), t = tObs, y = y, quantile = q)
  rstan::sampling(m, data = sDat, iter = iter, seed = seed, chains=4, refresh = 100)  
}, mc.cores=length(quants))

save.image(file="gpQuantAug.RData")

###################
f <- sapply(fits, function(fit) apply(extract(fit, "f")$f, 2, mean))
df <- sapply(fits, function(fit) apply(extract(fit, "df")$df, 2, mean))
QTDI <- sapply(fits, function(fit) apply(extract(fit, "df")$df, 2, function(q) mean(q > 0)))

band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}
col <- fields::tim.colors(length(quants))

pdf("cutiepie.pdf", width = 10, height=5)
par(mfrow=c(1,3), bty="n")
plot(tObs, y, xlab="t", ylab="E[Q | Y]", pch=19)
for(i in 1:5) {
  lines(tObs, f[,i], col=col[i], lwd=2)
}
legend("topright", paste(quants*100, "%", sep=""), col = col, lwd=2, bty="n", cex=0.7)

plot(tObs, rep(0, 25), type="n", ylim=c(-10,10), xlab="t", ylab="E[dQ | Y]")
for(i in 1:5) {
  lines(tObs, df[,i], col=col[i], lwd=2)
}
abline(h = 0, lty=2)

plot(0, 0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="t", ylab="Q-TDI")
for(i in 1:5) {
  lines(tObs, QTDI[,i], col=col[i], lwd=2)
}
abline(h = 0.5, lty=2)
dev.off()


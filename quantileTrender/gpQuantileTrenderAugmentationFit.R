rm(list=ls())

library(rstan)
library(parallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

m <- rstan::stan_model("gpQuantileTrenderAugmentation.stan")

tObs <- seq(0, 1, length.out=50)
set.seed(12345)
mu <- sin(2*pi*tObs)
sigma <- 0.5*(tObs^2 + 0.1*tObs)
y <- mu + rnorm(length(tObs), rep(0, length(tObs)), sigma)

plot(tObs, y, pch=19, cex=0.8)
lines(tObs, mu)
iter <- 2000
seed <- 12345

quants <- c(0.1, 0.25, 0.5, 0.75, 0.9)
fits <- mclapply(quants, function(q) {
  sDat <- list(n = length(tObs), t = tObs, y = y, quantile = q)
  rstan::sampling(m, data = sDat, iter = iter, seed = seed, chains=4, refresh = 100)  
}, mc.cores=length(quants))
save.image(file="gpQuantAug2.RData")

quants <- seq(0.1, 0.9, length.out=17)
fits <- mclapply(quants, function(q) {
  sDat <- list(n = length(tObs), t = tObs, y = y, quantile = q)
  rstan::sampling(m, data = sDat, iter = iter, seed = seed, chains=4, refresh = 100)  
}, mc.cores=length(quants))
save.image(file="gpQuantAug3.RData")

###################
rm(list=ls())
load("gpQuantAug2.RData")

band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

col <- fields::tim.colors(length(quants))

par(mfrow=c(1,3), bty="n")
plot(tObs, y, pch=19, ylim=c(-1.5, 1.5), xlab="Time", ylab="Q | Y")
for(i in 1:length(quants)) {
  f <- extract(fits[[i]], "f")$f

  band(tObs, apply(f, 2, quantile, prob=0.9),
       apply(f, 2, quantile, prob=0.1), 
       col=adjustcolor(col[i], alpha.f=0.2))
  #lines(tObs, apply(f, 2, mean), col=col[i], lwd=2)
}

plot(tObs, y, pch=19, type="n", ylim=c(-10, 10), xlab="Time", ylab="dQ | Y")
for(i in 1:length(quants)) {
  df <- extract(fits[[i]], "df")$df
  
  band(tObs, apply(df, 2, quantile, prob=0.9),
       apply(df, 2, quantile, prob=0.1), 
       col=adjustcolor(col[i], alpha.f=0.2))
  #lines(tObs, apply(df, 2, mean), col=col[i], lwd=2)
}
abline(h=0, lty=2)

plot(tObs, y, pch=19, type="n", ylim=c(0, 1), xlab="Time", ylab="Q-TDI")
for(i in 1:length(quants)) {
  df <- extract(fits[[i]], "df")$df
  lines(tObs, apply(df, 2, function(q) mean(q > 0)), col=col[i], lwd=2)
}
abline(h=0.5, lty=2)
#############################

rm(list=ls())
load("gpQuantAug3.RData")
length(fits)
heatm <- matrix(NA, length(tObs), length(quants))
heatm
fit <- fits[[1]]
test <- sapply(fits, function(fit) apply(extract(fit, "df")$df, 2, median))
test <- sapply(fits, function(fit) apply(extract(fit, "df")$df, 2, function(q) mean(q > 0)))
fields::image.plot(tObs, quants, test)


##############
f <- t(extract(fit, "f")$f)
df <- t(extract(fit, "df")$df)
test <- sapply(1:4000, function(i) predict(smooth.spline(tObs, f[,i], all.knots=TRUE), deriv=1)$y)

matplot(tObs, f, type="l", lty=1)
matplot(tObs, df, type="l", lty=1)
matplot(tObs, test, type="l", lty=1)
matplot(tObs, test - df, type="l", lty=1)

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


rm(list=ls())
library(pracma)
library(rootSolve)

getStats <- function(sim) {
  cat("Calculating for f\n")
  f <- sapply(1:length(sim), function(r) {
    pracma::trapz(sim[[r]]$tPred, sim[[r]]$f - sim[[r]]$fMu)  
  })
  
  cat("Calculating for df\n")
  df <- sapply(1:length(sim), function(r) {
    pracma::trapz(sim[[r]]$tPred, sim[[r]]$df - sim[[r]]$dfMu)  
  })
  
  cat("Calculating for TDI\n")
  TDI <- sapply(1:length(sim), function(r) {
    pracma::trapz(sim[[r]]$tPred, as.numeric(sim[[r]]$df > 0) - sim[[r]]$TDI)
  })
  
  cat("Calculating for ETI\n")
  ETI <- sapply(1:length(sim), function(r) {
    #Roots of true df
    roots <- rootSolve::uniroot.all(approxfun(sim[[r]]$tPred, sim[[r]]$df), c(0,1))
    if(length(roots) == 0) { #special case with no crossings
      ETI_true <- function(t) 0*t 
    } else {
      #True counting process
      ETI_true <- approxfun(c(roots), 1:(length(roots)), 
                            yleft = 0, yright = length(roots), method="constant")
    }
    #Difference between true counting process and cumulative estimated dETI
    dif <- ETI_true(sim[[r]]$tPred) - pracma::cumtrapz(sim[[r]]$tPred, sim[[r]]$ETI)[,1]
    #Integrate the difference on [0;1]
    pracma::trapz(sim[[r]]$tPred, dif)
  })
  
  list(f = f, df = df, TDI = TDI, ETI = ETI)
}

getSamples <- function(sim, r = 100) {
  list(t = sim[[1]]$t, tPred = sim[[1]]$tPred,
       f = sapply(sim[1:r], function(q) q$f),
       df = sapply(sim[1:r], function(q) q$df),
       y = sapply(sim[1:r], function(q) q$y))
}

##########
# n = 25
##########
load("sim25_res.RData")
sim25_1 <- sim25_1[!sapply(sim25_1, function(q) is.null(q$fMu))]
sim25_2 <- sim25_2[!sapply(sim25_2, function(q) is.null(q$fMu))]
sim25_3 <- sim25_3[!sapply(sim25_3, function(q) is.null(q$fMu))]
sim25_4 <- sim25_4[!sapply(sim25_4, function(q) is.null(q$fMu))]

stat_25_1 <- getStats(sim25_1)
stat_25_2 <- getStats(sim25_2)
stat_25_3 <- getStats(sim25_3)
stat_25_4 <- getStats(sim25_4)

samps_25_1 <- getSamples(sim25_1)
samps_25_2 <- getSamples(sim25_2)
samps_25_3 <- getSamples(sim25_3)
samps_25_4 <- getSamples(sim25_4)

##########
# n = 50
##########
load("sim50_res.RData")
sim50_1 <- sim50_1[!sapply(sim50_1, function(q) is.null(q$fMu))]
sim50_2 <- sim50_2[!sapply(sim50_2, function(q) is.null(q$fMu))]
sim50_3 <- sim50_3[!sapply(sim50_3, function(q) is.null(q$fMu))]
sim50_4 <- sim50_4[!sapply(sim50_4, function(q) is.null(q$fMu))]

stat_50_1 <- getStats(sim50_1)
stat_50_2 <- getStats(sim50_2)
stat_50_3 <- getStats(sim50_3)
stat_50_4 <- getStats(sim50_4)

samps_50_1 <- getSamples(sim50_1)
samps_50_2 <- getSamples(sim50_2)
samps_50_3 <- getSamples(sim50_3)
samps_50_4 <- getSamples(sim50_4)

##########
# n = 100
##########
load("sim100_res.RData")
sim100_1 <- sim100_1[!sapply(sim100_1, function(q) is.null(q$fMu))]
sim100_2 <- sim100_2[!sapply(sim100_2, function(q) is.null(q$fMu))]
sim100_3 <- sim100_3[!sapply(sim100_3, function(q) is.null(q$fMu))]
sim100_4 <- sim100_4[!sapply(sim100_4, function(q) is.null(q$fMu))]

stat_100_1 <- getStats(sim100_1)
stat_100_2 <- getStats(sim100_2)
stat_100_3 <- getStats(sim100_3)
stat_100_4 <- getStats(sim100_4)

samps_100_1 <- getSamples(sim100_1)
samps_100_2 <- getSamples(sim100_2)
samps_100_3 <- getSamples(sim100_3)
samps_100_4 <- getSamples(sim100_4)



################################################################################
plotSamps <- function(l1, l2, l3, l4) {
  col <- fields::tim.colors(8)
  par(mfrow=c(4, 3), bty="n", mar=c(3,3,1,0), mgp=c(1.8,1,0))
  matplot(l1$tPred, l1$f, type="l", lty=1, ylim=c(-4, 4), col=col,
          xlab="t", ylab="f", main=expression(plain("f") ~~ sigma == 0.025 ~~ n == 25))
  matplot(l1$tPred, l1$df, type="l", lty=1, xlab="t", ylab="df", ylim=c(-15, 15),
          main=expression(plain("df") ~~ sigma == 0.025 ~~ n == 25), col=col)
  matplot(l1$t, l1$y, type="l", lty=1, xlab="t", ylab="Y", ylim=c(-4,4),
          main=expression(plain("Y") ~~ sigma == 0.025 ~~ n == 25), col=col)
  
  matplot(l2$tPred, l2$f, type="l", lty=1, ylim=c(-4,4), col=col,
          xlab="t", ylab="f", main=expression(plain("f") ~~ sigma == 0.05 ~~ n == 25))
  matplot(l2$tPred, l2$df, type="l", lty=1, xlab="t", ylab="df", ylim=c(-15,15),
          main=expression(plain("df") ~~ sigma == 0.05 ~~ n == 25), col=col)
  matplot(l2$t, l2$y, type="l", lty=1, xlab="t", ylab="Y", ylim=c(-4,4), col=col)
  
  matplot(l3$tPred, l3$f, type="l", lty=1, ylim=c(-4,4), col=col,
          xlab="t", ylab="f", main=expression(plain("f") ~~ sigma == 0.1 ~~ n == 25))
  matplot(l3$tPred, l3$df, type="l", lty=1, xlab="t", ylab="df", ylim=c(-15,15),
          main=expression(plain("df") ~~ sigma == 0.1 ~~ n == 25), col=col)
  matplot(l3$t, l3$y, type="l", lty=1, xlab="t", ylab="Y", ylim=c(-4,4), col=col)
  
  matplot(l4$tPred, l4$f, type="l", lty=1, ylim=c(-4,4), col=col,
          xlab="t", ylab="f", main=expression(plain("f") ~~ sigma == 0.15 ~~ n == 25))
  matplot(l4$tPred, l4$df, type="l", lty=1, xlab="t", ylab="df", ylim=c(-15, 15),
          main=expression(plain("df") ~~ sigma == 0.15 ~~ n == 25), col=col)
  matplot(l4$t, l4$y, type="l", lty=1, xlab="t", ylab="Y", ylim=c(-4, 4), col=col)
}

plotSamps(samps_25_1, samps_25_2, samps_25_3, samps_25_4)
plotSamps(samps_50_1, samps_50_2, samps_50_3, samps_50_4)
plotSamps(samps_100_1, samps_100_2, samps_100_3, samps_100_4)
################



par(mfrow=c(2,2), bty="n", mgp=c(2,0.5,0), mar=c(2,2,2,0.5))
plot(0, 0, type="n", xlim=c(0, 16), ylim=c(-0.2, 0.2), xaxt="n",
     xlab="", ylab="", main="Integrated residual of E[f | Y]")
abline(h = 0, lty=2)
vioplot(stat_25_1$f, at = 1, add = TRUE, col="gray80")
vioplot(stat_50_1$f, at = 2, add = TRUE, col="gray60")
vioplot(stat_100_1$f, at = 3, add = TRUE, col="gray40")

vioplot(stat_25_2$f, at = 5, add = TRUE, col="gray80")
vioplot(stat_50_2$f, at = 6, add = TRUE, col="gray60")
vioplot(stat_100_2$f, at = 7, add = TRUE, col="gray40")

vioplot(stat_25_3$f, at = 9, add = TRUE, col="gray80")
vioplot(stat_50_3$f, at = 10, add = TRUE, col="gray60")
vioplot(stat_100_3$f, at = 11, add = TRUE, col="gray40")

vioplot(stat_25_4$f, at = 13, add = TRUE, col="gray80")
vioplot(stat_50_4$f, at = 14, add = TRUE, col="gray60")
vioplot(stat_100_4$f, at = 15, add = TRUE, col="gray40")

axis(1, c(0, 2, 6, 10, 14, 16),
     c("", expression(sigma == 0.025), expression(sigma == 0.5), 
       expression(sigma == 0.1), expression(sigma == 0.15), ""))
legend("topleft", c(expression(n == 25), expression(n == 50), expression(n == 100)),
       col = c("gray80", "gray60", "gray40"),
       bty="n", pch = c(15, 15, 15), pt.cex=1.5)

plot(0, 0, type="n", xlim=c(0, 16), ylim=c(-1, 1), xaxt="n",
     xlab="", ylab="", main="Integrated residual of E[df | Y]")
abline(h = 0, lty=2)
vioplot(stat_25_1$df, at = 1, add = TRUE, col="gray80")
vioplot(stat_50_1$df, at = 2, add = TRUE, col="gray60")
vioplot(stat_100_1$df, at = 3, add = TRUE, col="gray40")

vioplot(stat_25_2$df, at = 5, add = TRUE, col="gray80")
vioplot(stat_50_2$df, at = 6, add = TRUE, col="gray60")
vioplot(stat_100_2$df, at = 7, add = TRUE, col="gray40")

vioplot(stat_25_3$df, at = 9, add = TRUE, col="gray80")
vioplot(stat_50_3$df, at = 10, add = TRUE, col="gray60")
vioplot(stat_100_3$df, at = 11, add = TRUE, col="gray40")

vioplot(stat_25_4$df, at = 13, add = TRUE, col="gray80")
vioplot(stat_50_4$df, at = 14, add = TRUE, col="gray60")
vioplot(stat_100_4$df, at = 15, add = TRUE, col="gray40")
axis(1, c(0, 2, 6, 10, 14, 16),
     c("", expression(sigma == 0.025), expression(sigma == 0.5), 
       expression(sigma == 0.1), expression(sigma == 0.15), ""))

plot(0, 0, type="n", xlim=c(0, 16), ylim=c(-0.6, 0.6), xaxt="n",
     xlab="", ylab="", main="Integrated residual of TDI")
abline(h = 0, lty=2)
vioplot(stat_25_1$TDI, at = 1, add = TRUE, col="gray80")
vioplot(stat_50_1$TDI, at = 2, add = TRUE, col="gray60")
vioplot(stat_100_1$TDI, at = 3, add = TRUE, col="gray40")

vioplot(stat_25_2$TDI, at = 5, add = TRUE, col="gray80")
vioplot(stat_50_2$TDI, at = 6, add = TRUE, col="gray60")
vioplot(stat_100_2$TDI, at = 7, add = TRUE, col="gray40")

vioplot(stat_25_3$TDI, at = 9, add = TRUE, col="gray80")
vioplot(stat_50_3$TDI, at = 10, add = TRUE, col="gray60")
vioplot(stat_100_3$TDI, at = 11, add = TRUE, col="gray40")

vioplot(stat_25_4$TDI, at = 13, add = TRUE, col="gray80")
vioplot(stat_50_4$TDI, at = 14, add = TRUE, col="gray60")
vioplot(stat_100_4$TDI, at = 15, add = TRUE, col="gray40")
axis(1, c(0, 2, 6, 10, 14, 16),
     c("", expression(sigma == 0.025), expression(sigma == 0.5), 
       expression(sigma == 0.1), expression(sigma == 0.15), ""))

plot(0, 0, type="n", xlim=c(0, 16), ylim=c(-30, 30), xaxt="n",
     xlab="", ylab="", main="Integrated residual of ETI")
abline(h = 0, lty=2)
vioplot(stat_25_1$ETI, at = 1, add = TRUE, col="gray80")
vioplot(stat_50_1$ETI, at = 2, add = TRUE, col="gray60")
vioplot(stat_100_1$ETI, at = 3, add = TRUE, col="gray40")

vioplot(stat_25_2$ETI, at = 5, add = TRUE, col="gray80")
vioplot(stat_50_2$ETI, at = 6, add = TRUE, col="gray60")
vioplot(stat_100_2$ETI, at = 7, add = TRUE, col="gray40")

vioplot(stat_25_3$ETI, at = 9, add = TRUE, col="gray80")
vioplot(stat_50_3$ETI, at = 10, add = TRUE, col="gray60")
vioplot(stat_100_3$ETI, at = 11, add = TRUE, col="gray40")

vioplot(stat_25_4$ETI, at = 13, add = TRUE, col="gray80")
vioplot(stat_50_4$ETI, at = 14, add = TRUE, col="gray60")
vioplot(stat_100_4$ETI, at = 15, add = TRUE, col="gray40")
axis(1, c(0, 2, 6, 10, 14, 16),
     c("", expression(sigma == 0.025), expression(sigma == 0.5), 
       expression(sigma == 0.1), expression(sigma == 0.15), ""))

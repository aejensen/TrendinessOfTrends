rm(list=ls())

load("~/Dropbox/igangværende projekter/TrendinessOfTrends/Stan/pred.RData")

band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

pdf("fitBayes.pdf", width = 8, height = 6)
par(mfrow=c(2,2), bty="n", mar = c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0), col.lab="white", col.axis="white", bg=NA)
plot(dat$t, dat$y, pch = 19, xaxt="n", ylim=c(20,40), xlab="Year", ylab="f | Y", type="n", yaxt="n")
axis(1, 1998:2018, col = "white")
axis(2, col = "white")
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.005), apply(pred[,,1], 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.025), apply(pred[,,1], 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.25), apply(pred[,,1], 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,1], 2, mean), lwd = 2)
points(dat$t, dat$y, pch = 19, cex=0.8)
legend("topleft", c("Mean", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"), 
       lwd = 2, bty="n", cex=0.7, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5, text.col="white")

plot(tPred, apply(pred[,,3], 2, mean), lwd = 2, type="n", ylim = c(-3, 3), xaxt="n", xlab="Year", ylab="df | Y", yaxt="n")
axis(1, 1998:2018, col = "white")
axis(2, col = "white")
band(tPred, apply(pred[,,3], 2, quantile, prob = 0.005), apply(pred[,,3], 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,3], 2, quantile, prob = 0.025), apply(pred[,,3], 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,3], 2, quantile, prob = 0.25), apply(pred[,,3], 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,3], 2, mean), lwd = 2)
abline(h = 0, lty = 2, col="white")
legend("topleft", c("Mean", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"), 
       lwd = 2, bty="n", cex=0.7, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5, text.col="white")

plot(tPred, t(pred[1,,5])*100, type="n", lty = 1, lwd = 2, xaxt="n", xlab="Year", ylab="TDI [%]", col="white", yaxt="n", ylim=c(0,100))
band(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.005), apply(pred[,,5]*100, 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.025), apply(pred[,,5]*100, 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.25), apply(pred[,,5]*100, 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,5]*100, 2, median), lwd = 2)
axis(1, 1998:2018, col = "white")
axis(2, col = "white")
abline(h = 50, lty = 2, col="white")
legend("topleft", c("Median", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"), 
       lwd = 2, bty="n", cex=0.7, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5, text.col="white")

plot(tPred, t(pred[1,,6]), type="n", lty = 1, lwd = 2, xaxt="n", xlab="Year", ylab="Local ETI", col="white", yaxt="n", ylim=c(0,2))
band(tPred, apply(pred[,,6], 2, quantile, prob = 0.005), apply(pred[,,6], 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,6], 2, quantile, prob = 0.025), apply(pred[,,6], 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,6], 2, quantile, prob = 0.25), apply(pred[,,6], 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,6], 2, median), lwd = 2)
axis(1, 1998:2018, col = "white")
axis(2, col = "white")
legend("topleft", c("Median", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"), 
       lwd = 2, bty="n", cex=0.7, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5, text.col="white")
dev.off()


########################################################
# Plot 2
########################################################
ETI1 <- sapply(1:dim(pred)[1], function(i) integrate(splinefun(tPred, pred[i,,6]), 1998, 2018)$value)
ETI2 <- sapply(1:dim(pred)[1], function(i) integrate(splinefun(tPred, pred[i,,6]), 2008, 2018)$value)

ETI1_d <- density(ETI1)
ETI1_d_f <- approxfun(ETI1_d$x, ETI1_d$y)
ETI1_seq1 <- seq(quantile(ETI1, 0.025), quantile(ETI1, 0.975), length.out=100)
ETI1_seq2 <- seq(quantile(ETI1, 0.25), quantile(ETI1, 0.75), length.out=100)

ETI2_d <- density(ETI2)
ETI2_d_f <- approxfun(ETI2_d$x, ETI2_d$y)
ETI2_seq1 <- seq(quantile(ETI2, 0.025), quantile(ETI2, 0.975), length.out=100)
ETI2_seq2 <- seq(quantile(ETI2, 0.25), quantile(ETI2, 0.75), length.out=100)

pdf("ETIplot.pdf", width = 8, height = 4)
par(mfrow=c(1,2), bty="n", mar =  c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0), bg=NA, col.lab="white", col.axis="white")
plot(0, 0, type="n", xlim=c(0,6), ylim=c(0,0.8), xlab="ETI", ylab="Posterior density", xaxt="n", yaxt="n")
axis(1, seq(0, 6, 1), cex.axis = 0.8, col = "white")
axis(2, cex.axis = 0.8, col = "white")
band(ETI1_seq1, rep(0, length(ETI1_seq1)), ETI1_d_f(ETI1_seq1), col="gray65")
band(ETI1_seq2, rep(0, length(ETI1_seq2)), ETI1_d_f(ETI1_seq2), col="gray45")
lines(rep(quantile(ETI1, 0.5), 2), c(0, ETI1_d_f(quantile(ETI1, 0.5))), lwd = 2, lty=2, col="black")
curve(ETI1_d_f(x), 0, 6, ylim=c(0,1), add=TRUE, lwd = 2, col="white")
title("ETI([1998; 2018])", font.main=1, col.main = "white")
legend("topleft", c("Median", "50%", "95%"), col = c("black", "gray45", "gray65"), 
       lwd = 2, bty="n", cex=0.9, lty=c(2,NA, NA), pch = c(NA, 15, 15), pt.cex=1.5, text.col="white")

plot(0, 0, type="n", xlim=c(0,3), ylim=c(0,2.5), xlab="ETI", ylab="Posterior density", xaxt="n", yaxt="n")
axis(1, seq(0, 3, 0.5), cex.axis = 0.8, col="white")
axis(2, cex.axis = 0.8, col="white")
band(ETI2_seq1, rep(0, length(ETI2_seq1)), ETI2_d_f(ETI2_seq1), col="gray65")
band(ETI2_seq2, rep(0, length(ETI2_seq2)), ETI2_d_f(ETI2_seq2), col="gray45")
lines(rep(quantile(ETI2, 0.5), 2), c(0, ETI2_d_f(quantile(ETI2, 0.5))), lwd = 2, lty=2)
curve(ETI2_d_f(x), 0, 3, ylim=c(0,2.5), add=TRUE, lwd = 2, col="white")
title("ETI([2008; 2018])", font.main=1, col.main="white")
legend("topleft", c("Median", "50%", "95%"), col = c("black", "gray45", "gray65"), 
       lwd = 2, bty="n", cex=0.9, lty=c(2,NA, NA), pch = c(NA, 15, 15), pt.cex=1.5, text.col="white")
dev.off()


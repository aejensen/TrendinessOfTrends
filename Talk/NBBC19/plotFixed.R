rm(list=ls())

load("~/Dropbox/igangværende projekter/TrendinessOfTrends/Stan/predFixed.RData")

band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

pdf("fitLogLik.pdf", width = 8, height = 6)
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

plot(tPred, t(pred[1,,5])*100, type="l", lty = 1, lwd = 2, xaxt="n", xlab="Year", ylab="TDI [%]", col="white", yaxt="n", ylim=c(0,100))
axis(1, 1998:2018, col = "white")
axis(2, col = "white")
abline(h = 50, lty = 2, col="white")

plot(tPred, t(pred[1,,6]), type="l", lty = 1, lwd = 2, xaxt="n", xlab="Year", ylab="Local ETI", col="white", yaxt="n", ylim=c(0,2))
axis(1, 1998:2018, col = "white")
axis(2, col = "white")
dev.off()

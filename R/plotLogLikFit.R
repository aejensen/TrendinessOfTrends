rm(list=ls())

load("../Stan/predFixed.RData") #not in repository due to file size

#TDI summary at selected years
rbind(2018-(0:5), round(approxfun(tPred, pred[1,,5]*100)(2018-(0:5)), 2))
uniroot(function(x) approxfun(tPred, pred[1,,5]*100)(x) - 50, c(2010, 2018))$root

optimize(approxfun(tPred, pred[1,,5]*100), c(2004, 2008), maximum=TRUE)

#ETI summary
round(integrate(splinefun(tPred, pred[1,,6]), 2008, 2018)$value, 2)
round(integrate(splinefun(tPred, pred[1,,6]), 1998, 2018)$value, 2)


########################################################
# Plot fit
########################################################
band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

pdf("../figures/likFitPlot.pdf", width = 8, height = 6)
par(mfrow=c(2,2), bty="n", mar = c(2.3, 2.6, 1.5, 0), mgp=c(1.3,0.4,0))
plot(dat$t, dat$y, pch = 19, xaxt="n", ylim=c(20,40), xlab="Year", 
     ylab=expression(f ~ "|" ~ Y~","~widehat(Theta^ML)), type="n", yaxt="n")
title("Latent function", font.main=1)
axis(1, 1998:2018, cex.axis = 0.7)
axis(2, cex.axis = 0.7)
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.005), apply(pred[,,1], 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.025), apply(pred[,,1], 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,1], 2, quantile, prob = 0.25), apply(pred[,,1], 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,1], 2, mean), lwd = 2)
points(dat$t, dat$y, pch = 19, cex=0.8)
legend(1997.5, 41.5, c("Mean", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"), 
       lwd = 2, bty="n", cex=0.9, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5, seg.len=1.3)

plot(tPred, apply(pred[,,3], 2, mean), lwd = 2, type="n", ylim = c(-3, 3), xaxt="n", xlab="Year", 
     ylab=expression(df ~ "|" ~ Y~","~widehat(Theta^ML)), yaxt="n")
title("Latent trend", font.main=1)
axis(1, 1998:2018, cex.axis = 0.7)
axis(2, cex.axis = 0.7)
band(tPred, apply(pred[,,3], 2, quantile, prob = 0.005), apply(pred[,,3], 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,3], 2, quantile, prob = 0.025), apply(pred[,,3], 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,3], 2, quantile, prob = 0.25), apply(pred[,,3], 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,3], 2, mean), lwd = 2)
abline(h = 0, lty = 2)
legend(1997.5, 3.4, c("Mean", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"), 
       lwd = 2, bty="n", cex=0.9, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5, seg.len=1.3)

plot(tPred, t(pred[1,,5])*100, type="l", lty = 1, lwd = 2, xaxt="n", xlab="Year", 
     ylab="Trend Direction Index [%]", ylim=c(0,100), yaxt="n")
axis(1, 1998:2018, cex.axis = 0.7)
axis(2, cex.axis = 0.7)
abline(h = 50, lty = 2)
#title("TDI(2018, Year - 2018)", font.main=1)
title(bquote(TDI ~ ("2018, Year - 2018" ~ "|" ~ widehat(Theta^ML))), font.main=1)

plot(tPred, t(pred[1,,6]), type="l", lty = 1, lwd = 2, xaxt="n", xlab="Year", 
     ylab="Local Expected Trend Instability", ylim=c(0,1.8), yaxt="n")
#title(bquote(dETI[Year] ~ "([1998; 2018])"), font.main=1)
title(bquote(dETI ~ "(Year" ~ "|" ~ widehat(Theta^ML) ~ ")"), font.main=1)
axis(1, 1998:2018, cex.axis = 0.7)
axis(2, seq(0, 1.8, length.out=7), cex.axis = 0.7)
dev.off()

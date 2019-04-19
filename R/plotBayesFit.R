rm(list=ls())

load("../Stan/pred.RData") #not in repository due to file size
  
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

#TDI summary
round(cbind(2018-(0:5), 
            approxfun(tPred, apply(pred[,,5]*100, 2, median))(2018-(0:5)),
            approxfun(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.025))(2018-(0:5)),
            approxfun(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.975))(2018-(0:5))), 2)

optimize(approxfun(tPred, apply(pred[,,5]*100, 2, median)), c(2004, 2008), maximum=TRUE)
approxfun(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.025))(2005.896)
approxfun(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.975))(2005.896)

fiftyCrossYear <- apply(pred[,,5]*100, 1, function(q) {
  uniroot(function(x) approxfun(tPred, q)(x) - 50, c(2010, 2018))$root
})
quantile(fiftyCrossYear, c(0.025, 0.5, 0.975))

#ETI summary
round(quantile(ETI2, c(0.025, 0.5, 0.975)), 2)
round(quantile(ETI1, c(0.025, 0.5, 0.975)), 2)

round(integrate(splinefun(tPred, pred[1,,6]), 2008, 2018)$value, 2)
round(integrate(splinefun(tPred, pred[1,,6]), 1998, 2018)$value, 2)

########################################################
# Plot fit
########################################################
band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}
  
pdf("../figures/bayesFitPlot.pdf", width = 8, height = 6)
par(mfrow=c(2,2), bty="n", mar =  c(2.3, 2.3, 1.5, 0), mgp=c(1.3,0.4,0))
plot(tPred, apply(pred[,,5]*100, 2, median), type="n", lty = 1, lwd = 2, 
     xaxt="n", xlab="Year", ylab="Trend Direction Index [%]", ylim=c(0,100), yaxt="n")
axis(1, 1998:2018, cex.axis = 0.7)
axis(2, cex.axis = 0.7)
#title("TDI(2018, Year - 2018)", font.main=1)
title(bquote(TDI ~ ("2018, Year - 2018" ~ "|" ~ widetilde(Theta))), font.main=1)
band(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.005), apply(pred[,,5]*100, 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.025), apply(pred[,,5]*100, 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,5]*100, 2, quantile, prob = 0.25), apply(pred[,,5]*100, 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,5]*100, 2, median), lwd = 2)
abline(h = 50, lty = 2)
legend("topleft", c("Median", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"),
       lwd = 2, bty="n", cex=0.9, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5, seg.len=1.3)
  
plot(tPred, t(pred[1,,6]), type="n", lty = 1, lwd = 2, xaxt="n", xlab="Year", 
       ylab="Local Expected Trend Instability", ylim=c(0,1.8), yaxt="n")
#title(bquote(dETI[Year] ~ "([1998; 2018])"), font.main=1)
title(bquote(dETI ~ "(Year" ~ "|" ~ widetilde(Theta) ~ ")"), font.main=1)
axis(1, 1998:2018, cex.axis = 0.7)
axis(2, seq(0, 1.8, length.out=7), cex.axis = 0.7)
band(tPred, apply(pred[,,6], 2, quantile, prob = 0.005), apply(pred[,,6], 2, quantile, prob = 0.995), col = "gray80")
band(tPred, apply(pred[,,6], 2, quantile, prob = 0.025), apply(pred[,,6], 2, quantile, prob = 0.975), col = "gray65")
band(tPred, apply(pred[,,6], 2, quantile, prob = 0.25), apply(pred[,,6], 2, quantile, prob = 0.75), col = "gray45")
lines(tPred, apply(pred[,,6], 2, mean), lwd = 2)
legend("topleft", c("Median", "50%", "95%", "99%"), col = c("black", "gray45", "gray65", "gray85"), 
       lwd = 2, bty="n", cex=0.9, lty = c(1, NA, NA, NA), pch = c(NA, 15, 15, 15), pt.cex=1.5, seg.len=1.3)
  
plot(0, 0, type="n", xlim=c(0,6), ylim=c(0,0.8), xlab="Expected Trend Instability", 
     ylab="Density", xaxt="n", yaxt="n")
axis(1, seq(0, 6, 1), cex.axis = 0.7)
axis(2, cex.axis = 0.7)
band(ETI1_seq1, rep(0, length(ETI1_seq1)), ETI1_d_f(ETI1_seq1), col="gray80")
band(ETI1_seq2, rep(0, length(ETI1_seq2)), ETI1_d_f(ETI1_seq2), col="gray65")
lines(rep(quantile(ETI1, 0.5), 2), c(0, ETI1_d_f(quantile(ETI1, 0.5))), lwd = 2, lty=1)
curve(ETI1_d_f(x), 0, 6, ylim=c(0,1), add=TRUE, lwd = 2, n = 500)
#title("ETI([1998; 2018])", font.main=1)
title(bquote(ETI ~ "([1998; 2018]" ~ "|" ~ widetilde(Theta) ~ ")"), font.main=1)
legend("topleft", c("Median", "50%", "95%"), col = c("black", "gray65", "gray80"), 
       lwd = 2, bty="n", cex=0.9, lty=c(1,NA, NA), pch = c(NA, 15, 15), pt.cex=1.5, seg.len=1.3)
  
plot(0, 0, type="n", xlim=c(0,3), ylim=c(0,2.5), xlab="Expected Trend Instability", 
     ylab="Density", xaxt="n", yaxt="n")
axis(1, seq(0, 3, 0.5), cex.axis = 0.7)
axis(2, cex.axis = 0.7)
band(ETI2_seq1, rep(0, length(ETI2_seq1)), ETI2_d_f(ETI2_seq1), col="gray80")
band(ETI2_seq2, rep(0, length(ETI2_seq2)), ETI2_d_f(ETI2_seq2), col="gray65")
lines(rep(quantile(ETI2, 0.5), 2), c(0, ETI2_d_f(quantile(ETI2, 0.5))), lwd = 2, lty=1)
curve(ETI2_d_f(x), 0, 3, ylim=c(0,2), add=TRUE, lwd = 2, n = 500)
#title("ETI([2008; 2018])", font.main=1)
title(bquote(ETI ~ "([2008; 2018]" ~ "|" ~ widetilde(Theta) ~ ")"), font.main=1)
legend("topleft", c("Median", "50%", "95%"), col = c("black", "gray65", "gray80"), 
       lwd = 2, bty="n", cex=0.9, lty=c(1,NA, NA), pch = c(NA, 15, 15), pt.cex=1.5, seg.len=1.3)
dev.off()

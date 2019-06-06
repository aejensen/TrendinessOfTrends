library(rstan)
library(MCMCglmm)

add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)
nPlot <- 2000

################################################
# Plot model 1
################################################
pdf("smoking1.pdf", height = 6, width = 8)
par(mfrow=c(2,2), bty="n", mar =  c(3, 3, 0.5, 0.5), mgp=c(1.8,0.4,0), col.lab="white", col.axis="white", bg=NA)
plot(p ~ year, smoking, pch=19, ylim=c(20, 40), xaxt="n", ylab="Proportion", xlab="Year", cex=1, yaxt="n", col="white")
lines(xPred, apply(pred1[,,1], 2, mean), lwd = 2, col = "white")
lines(xPred, apply(pred1[,,1], 2, quantile, prob = 0.025), lty = 2, lwd = 1, col = "white")
lines(xPred, apply(pred1[,,1], 2, quantile, prob = 0.975), lty = 2, lwd = 1, col = "white")
lines(xPred, apply(pred1[,,3], 2, quantile, prob = 0.025), lty = 3, lwd = 1, col = "white")
lines(xPred, apply(pred1[,,3], 2, quantile, prob = 0.975), lty = 3, lwd = 1, col = "white")
axis(1, seq(1998, 2018, 2), col = "white")
axis(2, col = "white")

matplot(xPred, t(pred1[1:nPlot,,1]), type="l", lty=1, ylim=c(20,40), xaxt="n",
        xlab="Year", ylab="f", col = add.alpha(2:6, 0.8), yaxt="n")
lines(xPred, apply(pred1[,,1], 2, mean), lwd = 4, col = "white")
axis(1, seq(1998, 2018, 2), col = "white")
axis(2, col = "white")

matplot(xPred, t(pred1[1:nPlot,,2]), type="l", lty=1, xlab="Year", xaxt="n",
        ylab="df", col = add.alpha(2:6, 0.8), ylim=c(-4,4), yaxt="n")
lines(xPred, apply(pred1[,,2], 2, mean), lwd = 4, col = "white")
axis(1, seq(1998, 2018, 2), col="white")
axis(2, col = "white")

matplot(xPred, t(pred1[1:nPlot,,4]), type="l", lty=1, xlab="Year", xaxt="n",
        ylab="TDI", col = add.alpha(2:6, 0.8), yaxt="n")
lines(xPred, apply(pred1[,,4], 2, mean), lwd = 4, col = "white")
axis(1, seq(1998, 2018, 2), col = "white")
axis(2, col = "white")
abline(h = 0.5, lty = 2, col = "white")
dev.off()

################################################
# Plot model 2
################################################
pdf("smoking2-1.pdf", height = 5, width = 8)
par(mfrow=c(1,2), bty="n", mar = c(3, 3, 2, 0), mgp=c(1.8,0.4,0), col.lab="white", col.axis="white", bg=NA, col.main="white")
plot(xPred, apply(extract(fit2, "wPred")$wPred, 2, median), ylim=c(0,1), type="l", xaxt="n", yaxt="n",
     xlab="Year", ylab="Posterior median", main = "P(G(t) = 1)", lwd = 2, font.main = 1, col="white")
axis(1, seq(1998, 2018, 2), col="white")
axis(2, col="white")
abline(h = 0.5, lty = 2, col="white")

plot(density(extract(fit2, "rho")$rho[,1]),  ylab="Posterior density estimate", ylim=c(0, 0.5), 
     main="Length-scale parameters", xlab="Length-scale", lwd = 2, font.main = 1, col="white", yaxt="n", xaxt="n")
lines(density(extract(fit2, "rho")$rho[,2]), lwd = 2, lty=2, col="white")
axis(1, col="white")
axis(2, col="white")
legend("topright", c("Class 1", "Class 2"), lty=1:2, bty="n", lwd = 2, col="white", text.col="white")
dev.off()


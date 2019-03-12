library(rstan)
library(MCMCglmm)

add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)
nPlot <- 3000

################################################
# Plot model 0
################################################
par(mfrow=c(2,2), bty="n")
plot(p ~ year, smoking, pch=19, ylim=c(20, 40), xaxt="n")
lines(xPred, apply(pred0[,,1], 2, mean), lwd = 3)
lines(xPred, apply(pred0[,,1], 2, quantile, prob = 0.025), lty = 2, lwd = 1)
lines(xPred, apply(pred0[,,1], 2, quantile, prob = 0.975), lty = 2, lwd = 1)
lines(xPred, apply(pred0[,,3], 2, quantile, prob = 0.025), lty = 3, lwd = 1)
lines(xPred, apply(pred0[,,3], 2, quantile, prob = 0.975), lty = 3, lwd = 1)
axis(1, smoking$year)

matplot(xPred, t(pred0[1:nPlot,,1]), type="l", lty=1, ylim=c(20,40), xaxt="n",
        xlab="year", ylab="Posterior f", col = add.alpha(1:6, 0.5))
lines(xPred, apply(pred0[,,1], 2, mean), lwd = 4)
axis(1, smoking$year)

matplot(xPred, t(pred0[1:nPlot,,2]), type="l", lty=1, xlab="year", xaxt="n",
        ylab="Posterior df", col = add.alpha(1:6, 0.5), ylim=c(-4,4))
lines(xPred, apply(pred0[,,2], 2, mean), lwd = 4)
axis(1, smoking$year)
abline(h = 0, lty = 2)

matplot(xPred, t(pred0[1:nPlot,,4]), type="l", lty=1, xlab="year", xaxt="n",
        ylab="TDI", col = add.alpha(1:6, 0.5))
lines(xPred, apply(pred0[,,4], 2, mean), lwd = 4)
axis(1, smoking$year)
abline(h = 0.5, lty = 2)

################################################
# Plot model 1
################################################
pdf("smoking1.pdf", height = 5, width = 8)
par(mfrow=c(2,2), bty="n", mar =  c(3, 3, 0.5, 0.5), mgp=c(1.8,0.4,0))
plot(p ~ year, smoking, pch=19, ylim=c(20, 40), xaxt="n", ylab="Proportion", xlab="Year", cex=0.6)
lines(xPred, apply(pred1[,,1], 2, mean), lwd = 2)
lines(xPred, apply(pred1[,,1], 2, quantile, prob = 0.025), lty = 2, lwd = 1)
lines(xPred, apply(pred1[,,1], 2, quantile, prob = 0.975), lty = 2, lwd = 1)
lines(xPred, apply(pred1[,,3], 2, quantile, prob = 0.025), lty = 3, lwd = 1)
lines(xPred, apply(pred1[,,3], 2, quantile, prob = 0.975), lty = 3, lwd = 1)
axis(1, seq(1998, 2018, 2))

matplot(xPred, t(pred1[1:nPlot,,1]), type="l", lty=1, ylim=c(20,40), xaxt="n",
        xlab="Year", ylab="f", col = add.alpha(2:6, 0.1))
lines(xPred, apply(pred1[,,1], 2, mean), lwd = 2)
axis(1, seq(1998, 2018, 2))

matplot(xPred, t(pred1[1:nPlot,,2]), type="l", lty=1, xlab="Year", xaxt="n",
        ylab="df", col = add.alpha(2:6, 0.1), ylim=c(-4,4))
lines(xPred, apply(pred1[,,2], 2, mean), lwd = 2)
axis(1, seq(1998, 2018, 2))
abline(h = 0, lty = 2)

matplot(xPred, t(pred1[1:nPlot,,4]), type="l", lty=1, xlab="Year", xaxt="n",
        ylab="TDI", col = add.alpha(2:6, 0.1))
lines(xPred, apply(pred1[,,4], 2, median), lwd = 2)
lines(xPred, apply(pred1[,,4], 2, quantile, prob = 0.25), lwd = 1, lty=2)
lines(xPred, apply(pred1[,,4], 2, quantile, prob = 0.75), lwd = 1, lty=2)
axis(1, seq(1998, 2018, 2))
abline(h = 0.5, lty = 2)
dev.off()

################################################
# Plot model 2
################################################
pdf("smoking2-1.pdf", height = 3.5, width = 8)
par(mfrow=c(1,2), bty="n", mar = c(3, 3, 2, 0), mgp=c(1.8,0.4,0))
plot(xPred, apply(extract(fit2, "wPred")$wPred, 2, median), ylim=c(0,1), type="l", xaxt="n",
     xlab="Year", ylab="Posterior median", main = "P(G(t) = 1)", lwd = 2, font.main = 1,)
axis(1, seq(1998, 2018, 2))
abline(h = 0.5, lty = 2)

plot(density(extract(fit2, "rho")$rho[,1]),  ylab="Posterior density estimate", ylim=c(0, 0.5), 
     main="Length-scale parameters", xlab="Length-scale", lwd = 2, font.main = 1,)
lines(density(extract(fit2, "rho")$rho[,2]), lwd = 2, lty=2)
legend("topright", c("Class 1", "Class 2"), lty=1:2, bty="n", lwd = 2)
dev.off()

posterior.mode(mcmc(extract(fit2, "rho")$rho[,1]))
posterior.mode(mcmc(extract(fit2, "rho")$rho[,2]))


#Fix this!
pdf("smoking2-2.pdf", height = 5, width = 8)
par(mfrow=c(2,2), bty="n", mar =  c(3, 3, 0.5, 0.5), mgp=c(1.8,0.4,0))
plot(p ~ year, smoking, pch=19, ylim=c(20, 40), xaxt="n", ylab="Proportion", xlab="Year", cex=0.6)
lines(xPred, apply(pred1[,,1], 2, mean), lwd = 2)
lines(xPred, apply(pred1[,,1], 2, quantile, prob = 0.025), lty = 2, lwd = 1)
lines(xPred, apply(pred1[,,1], 2, quantile, prob = 0.975), lty = 2, lwd = 1)
lines(xPred, apply(pred1[,,3], 2, quantile, prob = 0.025), lty = 3, lwd = 1)
lines(xPred, apply(pred1[,,3], 2, quantile, prob = 0.975), lty = 3, lwd = 1)
axis(1, seq(1998, 2018, 2))

matplot(xPred, t(pred1[1:nPlot,,1]), type="l", lty=1, ylim=c(20,40), xaxt="n",
        xlab="Year", ylab="f", col = add.alpha(1:6, 0.1))
lines(xPred, apply(pred1[,,1], 2, mean), lwd = 2)
axis(1, seq(1998, 2018, 2))

matplot(xPred, t(pred1[1:nPlot,,2]), type="l", lty=1, xlab="Year", xaxt="n",
        ylab="df", col = add.alpha(1:6, 0.1), ylim=c(-4,4))
lines(xPred, apply(pred1[,,2], 2, mean), lwd = 2)
axis(1, seq(1998, 2018, 2))

matplot(xPred, t(pred1[1:nPlot,,4]), type="l", lty=1, xlab="Year", xaxt="n",
        ylab="TDI", col = add.alpha(1:6, 0.1))
lines(xPred, apply(pred1[,,4], 2, mean), lwd = 2)
axis(1, seq(1998, 2018, 2))
abline(h = 0.5, lty = 2)
dev.off()

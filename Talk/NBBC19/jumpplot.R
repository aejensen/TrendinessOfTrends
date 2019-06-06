rm(list=ls())
library(numDeriv)
load("~/Dropbox/igangværende projekter/TrendinessOfTrends/data/smoking.RData")

smoking[21,] <- c(2009, 28.0)
smoking

pdf("jumpPlot.pdf", width = 8, height = 5)
par(bty="n", mar = c(2.3, 2.3, 0, 0), mgp=c(1.3,0.4,0), bg = NA, col.lab="white", col.axis="white")
plot(1999:2018, as.numeric(diff(smoking$p) > 0), type="s", xaxt="n", bty="n", 
     xlab="year", ylab=expression(1(Y(t[i+1]) - Y(t[i]) > 0)), col = "white", lwd = 2)
axis(1, 1999:2018, col = "white")
axis(2, col = "white")
dev.off()
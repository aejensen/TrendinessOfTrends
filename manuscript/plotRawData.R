rm(list=ls())

load("../Stan/smoking.RData")
dev.off()
pdf("rawDataPlot.pdf", width = 7, height = 4)
par(mar =  c(2.3, 2.3, 0.7, 0), mgp=c(1.3,0.4,0), bty="n")
plot(p ~ year, smoking, pch=19, ylim=c(20, 40), xaxt="n", bty="n", xlab="Year", ylab="Proportion [%]", 
     main = "Proportion of smokers in Denmark", font.main = 1, yaxt="n", cex.main=1, cex.lab = 0.8)
axis(1, 1998:2018, cex.axis = 0.8, cex.lab=0.2)
axis(2, cex.axis=0.8)
dev.off()
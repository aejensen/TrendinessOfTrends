rm(list=ls())

load("../applications/smoking/smoking.RData")
covid19 <- read.csv("../applications/covid19/dpc-covid19-ita-andamento-nazionale.csv")

dat <- data.frame(t = 0:(nrow(covid19)-1), 
                  y = covid19$nuovi_positivi,
                  date=as.Date(covid19$data))

pdf("../figures/rawDataPlotDouble.pdf", width = 8, height = 3)
par(mfrow=c(1,2), mar =  c(2.3, 2.3, 0.7, 0.5), mgp=c(1.3,0.4,0), bty="n")
plot(p ~ year, smoking, pch=19, ylim=c(20, 40), xaxt="n", bty="n", xlab="Year", ylab="Proportion [%]", 
     main = "Daily or occasional smokers in Denmark", font.main = 1, yaxt="n", cex.main=0.9, cex.lab = 0.8, cex=0.7)
axis(1, 1998:2018, cex.axis = 0.8, cex.lab=0.2)
axis(2, cex.axis=0.8)

plot(y ~ t, data=dat, pch=19, bty="n", xlab="Number of days since Feb 24th", ylab="New positives", 
     main = "New daily positive COVID-19 cases in Italy (2020)", xaxt="n", ylim=c(0,7000),
     font.main = 1, yaxt="n", cex.main=0.9, cex.lab = 0.8, cex=0.7, cex.axis=0.8, xlim=c(0,90))
axis(1, seq(0, 90, 10), cex.axis = 0.8, cex.lab=0.2)
axis(2, cex.axis = 0.8, cex.lab=0.2)
dev.off()

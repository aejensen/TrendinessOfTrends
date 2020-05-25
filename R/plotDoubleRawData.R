rm(list=ls())

setwd("~/ku/forskning/trendiness_of_trends/R/")

load("../data/smoking.RData")

#pdf("../figures/rawDataPlot.pdf", width = 6, height = 3)
par(mar =  c(2.3, 2.3, 0.7, 0), mgp=c(1.3,0.4,0), bty="n", mfrow=c(1, 2))
plot(p ~ year, smoking, pch=19, ylim=c(20, 40), xaxt="n", bty="n", xlab="Year", ylab="Proportion [%]", 
     main = "Daily or occasional smokers in Denmark", font.main = 1, yaxt="n", cex.main=1, cex.lab = 0.8, cex=0.7)
axis(1, 1998:2018, cex.axis = 0.8, cex.lab=0.2)
axis(2, cex.axis=0.8)



d <- read.csv("../applications/covid19/dpc-covid19-ita-andamento-nazionale.csv")

#We look at the number of new positives
dat <- data.frame(t = 0:(nrow(d)-1), 
                  y = d$nuovi_positivi,
                  date=as.Date(d$data))

plot(y ~ date, data=dat, pch=19, bty="n", xlab="Date", ylab="New cases", 
     main = "New daily positive Covid-19 cases in Italy (2020)", font.main = 1, yaxt="n", cex.main=1, cex.lab = 0.8, cex=0.7, cex.axis=0.8)
#axis(1, 1998:2018, cex.axis = 0.8, cex.lab=0.2)
axis(2, cex.axis=0.8)


#dev.off()
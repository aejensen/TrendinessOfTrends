rm(list=ls())
load("~/Dropbox/igangværende projekter/trend/rygevaner.RData")
pdf(file = "smokersData.pdf", width = 8, height = 6)
par(bg = NA, col.lab="white", col.axis="white")
#par(bg = "black", col.lab="white", col.axis="white")
plot(rygning$year, rygning$p, bty="n", pch=19, xaxt="n", xlab="Year", 
     ylab="Proportion", col="white", yaxt="n", ylim=c(20, 35))
title("Daily or occasional smokers in Denmark", col.main="white")
axis(2, col = "white")
axis(1, 1998:2018, col="white")
lines(c(2007.66, 2007.66), c(20, 35), col = "white", lty=2)
legend("topright", c("Law of smoke free environments"), lty=2, bty="n", col="white", text.col="white")
dev.off()

############################




par(bg = NA, col.lab="white", col.axis="white")
plot(fireworks$aar, fireworks$alvorlige_skader, bty="n", pch=19, xaxt="n", xlab="Year", 
     ylab="Number of serious firework injuries", col="white", yaxt="n", ylim=c(20, 100))
axis(2, seq(20, 100, 10), col = "white")
axis(1, fireworks$aar, col="white")



library(readxl)
hej <- read_excel("~/Desktop/hej.xlsx", col_names = FALSE)
hej$year <- seq(2013, 2019 - 1/12, length.out = 71)
plot(hej$year, hej$X__2, type="l", bty="n")
abline(h = 100)



plot(hej$X__2, type="l")

plot(100 * (hej$X__2 - hej$X__2[1]) / hej$X__2[1], type="l")
plot((hej$X__2 / hej$X__2[1]) * 100, type="l")

data$p
data

year <- c(1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018)
p <- c(34.6, 34.1, 33.5, 32.3, 31.0, 30.0, 27.1, 28.0, 27.7, 28.5, 28.0, 24.3, 23.4, 22.3, 22.6, 21.0, 22.5, 21.1, 21.6, 23.1)
plot(year, p, col = "white", pch = 19)

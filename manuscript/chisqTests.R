rm(list=ls())
load("../../../trend/ildipiben.RData")
data$non <- data$n - data$x
data[1:6,]

test1 <- chisq.test(matrix(c(1158, 1106, 3859, 4018), 2, 2)) #2018 vs 2017
test2 <- chisq.test(matrix(c(1158, 1070, 3859, 3992), 2, 2)) #2018 vs 2016
test3 <- chisq.test(matrix(c(1158, 1134, 3859, 3908), 2, 2)) #2018 vs 2015
test4 <- chisq.test(matrix(c(1158, 1060, 3859, 3990), 2, 2)) #2018 vs 2014
test5 <- chisq.test(matrix(c(1158, 1133, 3859, 3882), 2, 2)) #2018 vs 2013

plot(c(test1$p.value, test2$p.value, test3$p.value, test4$p.value, test5$p.value), ylim=c(0,1), pch=19)
abline(h = 0.05)

round(c(test1$p.value, test2$p.value, test3$p.value, test4$p.value, test5$p.value), 3)

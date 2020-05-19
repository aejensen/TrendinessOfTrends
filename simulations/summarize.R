rm(list=ls())
library(pracma)
library(rootSolve)

getStats <- function(sim) {
  cat("Calculating for f\n")
  f <- sapply(1:length(sim), function(r) {
    pracma::trapz(sim[[r]]$tPred, sim[[r]]$f - sim[[r]]$fMu)  
  })
  
  cat("Calculating for df\n")
  df <- sapply(1:length(sim), function(r) {
    pracma::trapz(sim[[r]]$tPred, sim[[r]]$df - sim[[r]]$dfMu)  
  })
  
  cat("Calculating for TDI\n")
  TDI <- sapply(1:length(sim), function(r) {
    pracma::trapz(sim[[r]]$tPred, as.numeric(sim[[r]]$df > 0) - sim[[r]]$TDI)
  })
  
  cat("Calculating for ETI\n")
  ETI <- sapply(1:length(sim), function(r) {
    #Roots of true df
    roots <- rootSolve::uniroot.all(approxfun(sim[[r]]$tPred, sim[[r]]$df), c(0,1))
    if(length(roots) == 0) { #special case with no crossings
      ETI_true <- function(t) 0*t 
    } else {
      #True counting process
      ETI_true <- approxfun(c(roots), 1:(length(roots)), 
                            yleft = 0, yright = length(roots), method="constant")
    }
    #Difference between true counting process and cumulative estimated dETI
    dif <- ETI_true(sim[[r]]$tPred) - pracma::cumtrapz(sim[[r]]$tPred, sim[[r]]$ETI)[,1]
    #Integrate the difference on [0;1]
    pracma::trapz(sim[[r]]$tPred, dif)
  })
  
  list(f = f, df = df, TDI = TDI, ETI = ETI)
}

getSamples <- function(sim, r = 100) {
  list(t = sim[[1]]$t, tPred = sim[[1]]$tPred,
       f = sapply(sim[1:r], function(q) q$f),
       df = sapply(sim[1:r], function(q) q$df),
       y = sapply(sim[1:r], function(q) q$y))
}

##########
# n = 25
##########
load("sim25_res.RData")
sim25_1 <- sim25_1[!sapply(sim25_1, function(q) is.null(q$fMu))]
sim25_2 <- sim25_2[!sapply(sim25_2, function(q) is.null(q$fMu))]
sim25_3 <- sim25_3[!sapply(sim25_3, function(q) is.null(q$fMu))]
sim25_4 <- sim25_4[!sapply(sim25_4, function(q) is.null(q$fMu))]

stat_25_1 <- getStats(sim25_1)
stat_25_2 <- getStats(sim25_2)
stat_25_3 <- getStats(sim25_3)
stat_25_4 <- getStats(sim25_4)

samps_25_1 <- getSamples(sim25_1)
samps_25_2 <- getSamples(sim25_2)
samps_25_3 <- getSamples(sim25_3)
samps_25_4 <- getSamples(sim25_4)

##########
# n = 50
##########
load("sim50_res.RData")
sim50_1 <- sim50_1[!sapply(sim50_1, function(q) is.null(q$fMu))]
sim50_2 <- sim50_2[!sapply(sim50_2, function(q) is.null(q$fMu))]
sim50_3 <- sim50_3[!sapply(sim50_3, function(q) is.null(q$fMu))]
sim50_4 <- sim50_4[!sapply(sim50_4, function(q) is.null(q$fMu))]

stat_50_1 <- getStats(sim50_1)
stat_50_2 <- getStats(sim50_2)
stat_50_3 <- getStats(sim50_3)
stat_50_4 <- getStats(sim50_4)

samps_50_1 <- getSamples(sim50_1)
samps_50_2 <- getSamples(sim50_2)
samps_50_3 <- getSamples(sim50_3)
samps_50_4 <- getSamples(sim50_4)

##########
# n = 100
##########
load("sim100_res.RData")
sim100_1 <- sim100_1[!sapply(sim100_1, function(q) is.null(q$fMu))]
sim100_2 <- sim100_2[!sapply(sim100_2, function(q) is.null(q$fMu))]
sim100_3 <- sim100_3[!sapply(sim100_3, function(q) is.null(q$fMu))]
sim100_4 <- sim100_4[!sapply(sim100_4, function(q) is.null(q$fMu))]

stat_100_1 <- getStats(sim100_1)
stat_100_2 <- getStats(sim100_2)
stat_100_3 <- getStats(sim100_3)
stat_100_4 <- getStats(sim100_4)

samps_100_1 <- getSamples(sim100_1)
samps_100_2 <- getSamples(sim100_2)
samps_100_3 <- getSamples(sim100_3)
samps_100_4 <- getSamples(sim100_4)

save(stat_25_1, stat_25_2, stat_25_3, stat_25_4,
     samps_25_1, samps_25_2, samps_25_3, samps_25_4,
     stat_50_1, stat_50_2, stat_50_3, stat_50_4,
     samps_50_1, samps_50_2, samps_50_3, samps_50_4,
     stat_100_1, stat_100_2, stat_100_3, stat_100_4,
     samps_100_1, samps_100_2, samps_100_3, samps_100_4,
     file = "simulationSummaries.RData")

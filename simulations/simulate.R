rm(list=ls())

source("functions.R")

set.seed(12345)

R <- 10000

########################################################
#n = 25
########################################################
set.seed(12345)
sim25_1 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(25, 0.025)
}, mc.cores = 64)

set.seed(12345)
sim25_2 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(25, 0.05)
}, mc.cores = 64)

set.seed(12345)
sim25_3 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(25, 0.1)
}, mc.cores = 64)

set.seed(12345)
sim25_4 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(25, 0.15)
}, mc.cores = 64)

save(sim25_1, sim25_2, sim25_3, sim25_4,
     file="sim25_res.RData")


########################################################
#n = 50
########################################################
set.seed(12345)
sim50_1 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(50, 0.025)
}, mc.cores = 64)

set.seed(12345)
sim50_2 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(50, 0.05)
}, mc.cores = 64)

set.seed(12345)
sim50_3 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(50, 0.1)
}, mc.cores = 64)

set.seed(12345)
sim50_4 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(50, 0.15)
}, mc.cores = 64)

save(sim50_1, sim50_2, sim50_3, sim50_4,
     file="sim50_res.RData")

########################################################
#n = 100
########################################################
set.seed(12345)
sim100_1 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(100, 0.025)
}, mc.cores = 64)

set.seed(12345)
sim100_2 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(100, 0.05)
}, mc.cores = 64)

set.seed(12345)
sim100_3 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(100, 0.1)
}, mc.cores = 64)

set.seed(12345)
sim100_4 <- mclapply(1:R, function(r) {
  cat("Simulation", r, "\n")
  doSim(100, 0.15)
}, mc.cores = 64)

save(sim100_1, sim100_2, sim100_3, sim100_4,
     file="sim100_res.RData")

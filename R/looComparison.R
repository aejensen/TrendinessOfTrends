library(DEoptim)

rm(list=ls())

load("../Data/smoking.RData")
dat <- data.frame(t = smoking$year, y = smoking$p)
dat$t <- dat$t - mean(dat$t)
tPred <- tPred - mean(dat$t)

seCov <- function(s, t, alpha, rho) {
  alpha^2 * exp(-(s-t)^2 / (2*rho^2))
}

matern32Cov <- function(s, t, alpha, rho) {
  r <- abs(s-t)
  alpha^2 * (1 + sqrt(3)*r/rho) * exp(-sqrt(3)*r/rho)
}

matern52Cov <- function(s, t, alpha, rho) {
  r <- abs(s-t)
  alpha^2 * (1 + sqrt(5)*r/rho + (5*r^2)/(3*rho^2)) * exp(-sqrt(5)*r/rho)
}

rqCov <- function(s, t, alpha, rho, nu) {
  alpha^2 * (1 + (s-t)^2 / (2 * nu * rho^2))^(-nu)
}

ctl <- DEoptim.control(itermax = 1000, trace = 100)

##################################################################
# SE covariance, constant mean
##################################################################
par.se.const <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
      mu <- rep(par[1], nrow(dTrain))
      cMat <- outer(dTrain$t, dTrain$t, seCov, par[2], par[3]) + diag(par[4]^2 , nrow(dTrain))
      -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
    }, lower = c(0,0,0,0), 
       upper = c(50,50,50,50), 
       control = ctl)
  opt$optim$bestmem
}))

msep.se.const <- mean(sapply(1:20, function(i) {
  coe <- par.se.const[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- rep(coe[1], 1)
  muY <- rep(coe[1], length(dTrain$t))
  K <- outer(dTrain$t, dTrain$t, seCov, coe[2], coe[3]) + diag(coe[4]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, seCov, coe[2], coe[3])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))

##################################################################
# RQ covariance, constant mean
##################################################################
par.rq.const <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
    mu <- rep(par[1], nrow(dTrain))
    cMat <- outer(dTrain$t, dTrain$t, rqCov, par[2], par[3], par[4]) + diag(par[5]^2 , nrow(dTrain))
    -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
  }, lower = c(0,0,0,0,0), 
  upper = c(50,10,10,5,5), 
  control = ctl)
  opt$optim$bestmem
}))

msep.rq.const <- mean(sapply(1:20, function(i) {
  coe <- par.rq.const[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- rep(coe[1], 1)
  muY <- rep(coe[1], length(dTrain$t))
  K <- outer(dTrain$t, dTrain$t, rqCov, coe[2], coe[3], coe[4]) + diag(coe[5]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, rqCov, coe[2], coe[3], coe[4])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))

##################################################################
# Matern 5/2 covariance, constant mean
##################################################################
par.matern52.const <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
    mu <- rep(par[1], nrow(dTrain))
    cMat <- outer(dTrain$t, dTrain$t, matern52Cov, par[2], par[3]) + diag(par[4]^2 , nrow(dTrain))
    -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
  }, lower = c(0,0,0,0), 
  upper = c(50,50,50,50), 
  control = ctl)
  opt$optim$bestmem
}))

msep.matern52.const <- mean(sapply(1:20, function(i) {
  coe <- par.matern52.const[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- rep(coe[1], 1)
  muY <- rep(coe[1], length(dTrain$t))
  K <- outer(dTrain$t, dTrain$t, matern52Cov, coe[2], coe[3]) + diag(coe[4]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, matern52Cov, coe[2], coe[3])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))

##################################################################
# Matern 3/2 covariance, constant mean
##################################################################
par.matern32.const <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
    mu <- rep(par[1], nrow(dTrain))
    cMat <- outer(dTrain$t, dTrain$t, matern32Cov, par[2], par[3]) + diag(par[4]^2 , nrow(dTrain))
    -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
  }, lower = c(0,0,0,0), 
  upper = c(50,50,50,50), 
  control = ctl)
  opt$optim$bestmem
}))

msep.matern32.const <- mean(sapply(1:20, function(i) {
  coe <- par.matern32.const[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- rep(coe[1], 1)
  muY <- rep(coe[1], length(dTrain$t))
  K <- outer(dTrain$t, dTrain$t, matern32Cov, coe[2], coe[3]) + diag(coe[4]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, matern32Cov, coe[2], coe[3])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))

##################################################################
# SE covariance, linear mean
##################################################################
par.se.lin <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
    mu <- par[1] + par[2] * dTrain$t
    cMat <- outer(dTrain$t, dTrain$t, seCov, par[3], par[4]) + diag(par[5]^2 , nrow(dTrain))
    -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
  }, lower = c(0,-2,0, 0, 0), 
  upper = c(50,2, 50, 50,50), 
  control = ctl)
  opt$optim$bestmem
}))

msep.se.lin <- mean(sapply(1:20, function(i) {
  coe <- par.se.lin[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- coe[1] + coe[2] * dValid$t
  muY <- coe[1] + coe[2] * dTrain$t
  K <- outer(dTrain$t, dTrain$t, seCov, coe[3], coe[4]) + diag(coe[5]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, seCov, coe[3], coe[4])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))

##################################################################
# RQ covariance, linear mean
##################################################################
par.rq.lin <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
    mu <- par[1] + par[2] * dTrain$t
    cMat <- outer(dTrain$t, dTrain$t, rqCov, par[3], par[4], par[5]) + diag(par[6]^2 , nrow(dTrain))
    -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
  }, lower = c(0,-2, 0, 0, 0, 0), 
  upper = c(50, 2, 10, 5, 500, 5), 
  control = ctl)
  opt$optim$bestmem
}))

msep.rq.lin <- mean(sapply(1:20, function(i) {
  coe <- par.rq.lin[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- coe[1] + coe[2] * dValid$t
  muY <- coe[1] + coe[2] * dTrain$t
  K <- outer(dTrain$t, dTrain$t, rqCov, coe[3], coe[4], coe[5]) + diag(coe[6]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, rqCov, coe[3], coe[4], coe[5])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))


##################################################################
# Matern 5/2 covariance, linear mean
##################################################################
par.matern52.lin <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
    mu <- par[1] + par[2] * dTrain$t
    cMat <- outer(dTrain$t, dTrain$t, matern52Cov, par[3], par[4]) + diag(par[5]^2 , nrow(dTrain))
    -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
  }, lower = c(0,-2,0, 0, 0), 
  upper = c(50,2, 50, 50,50), 
  control = ctl)
  opt$optim$bestmem
}))

msep.matern52.lin <- mean(sapply(1:20, function(i) {
  coe <- par.matern52.lin[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- coe[1] + coe[2] * dValid$t
  muY <- coe[1] + coe[2] * dTrain$t
  K <- outer(dTrain$t, dTrain$t, matern52Cov, coe[3], coe[4]) + diag(coe[5]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, matern52Cov, coe[3], coe[4])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))

##################################################################
# Matern 3/2 covariance, linear mean
##################################################################
par.matern32.lin <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
    mu <- par[1] + par[2] * dTrain$t
    cMat <- outer(dTrain$t, dTrain$t, matern32Cov, par[3], par[4]) + diag(par[5]^2 , nrow(dTrain))
    -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
  }, lower = c(0,-2,0, 0, 0), 
  upper = c(50,2, 50, 50,50), 
  control = ctl)
  opt$optim$bestmem
}))

msep.matern32.lin <- mean(sapply(1:20, function(i) {
  coe <- par.matern32.lin[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- coe[1] + coe[2] * dValid$t
  muY <- coe[1] + coe[2] * dTrain$t
  K <- outer(dTrain$t, dTrain$t, matern32Cov, coe[3], coe[4]) + diag(coe[5]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, matern32Cov, coe[3], coe[4])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))

##################################################################
# SE covariance, quadratic mean
##################################################################
par.se.quad <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
    mu <- par[1] + par[2] * dTrain$t + par[3] * dTrain$t^2
    cMat <- outer(dTrain$t, dTrain$t, seCov, par[4], par[5]) + diag(par[6]^2 , nrow(dTrain))
    -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
  }, lower = c(0,-2,-2, 0, 0, 0), 
  upper = c(50, 2, 2, 50, 50,50), 
  control = ctl)
  opt$optim$bestmem
}))

msep.se.quad <- mean(sapply(1:20, function(i) {
  coe <- par.se.quad[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- coe[1] + coe[2] * dValid$t + coe[3] * dValid$t^2
  muY <- coe[1] + coe[2] * dTrain$t + coe[3] * dTrain$t^2
  K <- outer(dTrain$t, dTrain$t, seCov, coe[4], coe[5]) + diag(coe[6]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, seCov, coe[4], coe[5])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))

##################################################################
# RQ covariance, quadratic mean
##################################################################
par.rq.quad <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
    mu <- par[1] + par[2] * dTrain$t + par[3] * dTrain$t^2
    cMat <- outer(dTrain$t, dTrain$t, rqCov, par[4], par[5], par[6]) + diag(par[7]^2 , nrow(dTrain))
    -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
  }, lower = c(0,-2,-2, 0, 0, 0, 0), 
  upper = c(50, 2, 2, 50, 50, 200, 50), 
  control = ctl)
  opt$optim$bestmem
}))

msep.rq.quad <- mean(sapply(1:20, function(i) {
  coe <- par.rq.quad[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- coe[1] + coe[2] * dValid$t + coe[3] * dValid$t^2
  muY <- coe[1] + coe[2] * dTrain$t + coe[3] * dTrain$t^2
  K <- outer(dTrain$t, dTrain$t, rqCov, coe[4], coe[5], coe[6]) + diag(coe[7]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, rqCov, coe[4], coe[5], coe[6])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))

##################################################################
# Matern 5/2 covariance, quadratic mean
##################################################################
par.matern52.quad <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
    mu <- par[1] + par[2] * dTrain$t + par[3] * dTrain$t^2
    cMat <- outer(dTrain$t, dTrain$t, matern52Cov, par[4], par[5]) + diag(par[6]^2 , nrow(dTrain))
    -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
  }, lower = c(0,-2,-2, 0, 0, 0), 
  upper = c(50, 2, 2, 50, 50,50), 
  control = ctl)
  opt$optim$bestmem
}))

msep.matern52.quad <- mean(sapply(1:20, function(i) {
  coe <- par.matern52.quad[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- coe[1] + coe[2] * dValid$t + coe[3] * dValid$t^2
  muY <- coe[1] + coe[2] * dTrain$t + coe[3] * dTrain$t^2
  K <- outer(dTrain$t, dTrain$t, matern52Cov, coe[4], coe[5]) + diag(coe[6]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, matern52Cov, coe[4], coe[5])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))

##################################################################
# Matern 3/2 covariance, quadratic mean
##################################################################
par.matern32.quad <- t(sapply(1:20, function(i) {
  cat(i, "\n")
  dTrain <- dat[-i,]
  opt <- DEoptim(function(par) {
    mu <- par[1] + par[2] * dTrain$t + par[3] * dTrain$t^2
    cMat <- outer(dTrain$t, dTrain$t, matern32Cov, par[4], par[5]) + diag(par[6]^2 , nrow(dTrain))
    -mvtnorm::dmvnorm(dTrain$y, mu, cMat, log=TRUE)
  }, lower = c(0,-2,-2, 0, 0, 0), 
  upper = c(50, 2, 2, 50, 50,50), 
  control = ctl)
  opt$optim$bestmem
}))

msep.matern32.quad <- mean(sapply(1:20, function(i) {
  coe <- par.matern32.quad[1,]
  dTrain <- dat[-i,]
  dValid <- dat[i,]
  
  muF <- coe[1] + coe[2] * dValid$t + coe[3] * dValid$t^2
  muY <- coe[1] + coe[2] * dTrain$t + coe[3] * dTrain$t^2
  K <- outer(dTrain$t, dTrain$t, matern32Cov, coe[4], coe[5]) + diag(coe[6]^2, nrow(dTrain))
  K1 <- outer(dValid$t, dTrain$t, matern32Cov, coe[4], coe[5])
  pred <- muF + K1 %*% solve(K) %*% as.matrix(dTrain$y - muY, nrow(dTrain), 1)
  (as.numeric(pred) - dValid$y)^2
}))

##################################################################
# Collect results
##################################################################

#Note: nu saturates for RQ linear mean and RQ quad mean => SE cov
msepMat <- matrix(c(msep.se.const,
                    msep.rq.const,
                    msep.matern32.const,
                    msep.matern52.const,
                    msep.se.lin,
                    msep.rq.lin,
                    msep.matern32.lin,
                    msep.matern52.lin,
                    msep.se.quad,
                    msep.rq.quad,
                    msep.matern52.quad,
                    msep.matern32.quad))
rownames(msepMat) <- c("SE const", "RQ const", "Matern 3/2 const", "Matern 5/2 const",
                       "SE lin", "RQ lin", "Matern 3/2 lin", "Matern 5/2 lin",
                       "SE quad", "RQ quad", "Matern 3/2 quad", "Matern 5/2 quad")
data.frame(msep = msepMat[order(msepMat),])


res <- rbind(c(msep.se.const, msep.se.lin, msep.se.quad),
             c(msep.rq.const, msep.rq.lin, msep.rq.quad),
             c(msep.matern32.const, msep.matern32.lin, msep.matern32.quad),
             c(msep.matern52.const, msep.matern52.lin, msep.matern52.quad))
colnames(res) <- c("Const", "Lin", "Quad")
rownames(res) <- c("SE", "RQ", "Matern 3/2", "Matern 5/2")
round(t(res), 3)

save.image(file = "looComparison.RData")

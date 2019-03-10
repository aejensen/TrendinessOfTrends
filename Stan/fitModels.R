rm(list=ls())
library(rstan)
library(invgamma)
library(DEoptim)
library(gptrendR) #Available at github.com/aejensen/gptrendR
library(splines2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

getInvGammaParam <- function(min, max, prob = 0.01, upper = 20, itermax = 10^4) {
  opt <- DEoptim(function(par) {
    cond1 <- invgamma::pinvgamma(min, par[1], par[2])
    cond2 <- invgamma::pinvgamma(max, par[1], par[2])
    (abs(cond1 - prob/2)^2 + abs(cond2 - (1 - prob/2))^2)^2
  }, lower = c(0,0), upper = c(upper, upper),
  control = DEoptim.control(itermax = itermax, trace=1000))
  opt$optim$bestmem
} 

################################################
# Setup data
################################################
load("smoking.RData")
xPred <- seq(1998, 2018, 0.2)
dat <- list(n = nrow(smoking), x = smoking$year, y = smoking$p, nPred = length(xPred), xPred = xPred)

seed <- 123456
iter <- 25000

################################################
# Fit model 0 - fixed parameters
################################################
m0 <- stan_model("gptrendFixed.stan")

# Get the empirical Bayes estimates
gp <- gptrendR::GP(dat$x, dat$y, gptrendR::seCov)
#gpP <- predict(gp, xPred)
#plot(gpP$df$t, gpP$df$mu, type="l")

d0 <- dat
d0$mu <- 27.666537
d0$alpha <- 3.997589
d0$rho <- 3.242167
d0$sigma <- 0.633646

fit0 <- sampling(m0, data = d0, iter = iter, seed = seed, algorithm = "Fixed_param")
pred0 <- extract(fit0, "pred")$pred

################################################
# Fit model 1 - Bayesian model
################################################
m1 <- stan_model("gptrend.stan")

#Get inv-gamma scale and shape
c(min(diff(smoking$year)), diff(range(smoking$year)))
invGammaOpt <- getInvGammaParam(1, 20, prob = 0.01)
invGammaOpt
c(invgamma::qinvgamma(0.01/2, invGammaOpt[1], invGammaOpt[2]), 
  invgamma::qinvgamma(1-0.01/2, invGammaOpt[1], invGammaOpt[2]))

d1 <- dat
d1$mu_mu <- mean(smoking$p)
d1$rho_shape <- 3.548762
d1$rho_scale <- 10.221723

fit1 <- sampling(m1, data = d1, iter = iter, seed = seed)
pred1 <- extract(fit1, "pred")$pred

################################################
# Fit model 2 - Latent class model
################################################
m2 <- stan_model("gptrendLatentClass.stan")

df <- 15
degree <- 3
basis <- t(bSpline(dat$x, df = df, degree = degree, intercept = FALSE))
basisPred <- t(bSpline(xPred, df = df, degree = degree, intercept = FALSE))
basisD <- t(dbs(dat$x, derivs = 1, df = df, degree = degree, incercept = FALSE))
basisPredD <- t(dbs(xPred, derivs = 1, df = df, degree = degree, incercept = FALSE))

d2 <- dat
d2$mu_mu <- mean(smoking$p)
d2$rho_shape <- 3.548762
d2$rho_scale <- 10.221723
d2$basis <- basis
d2$basisD <- basisD
d2$basisPred <- basisPred
d2$basisPredD <- basisPredD
d2$df <- nrow(basis)

fit2 <- sampling(m2, data = d2, iter = iter, seed = 1234)
pred2 <- extract(fit2, "pred")$pred

tObs <- seq(0, 1, length.out=25)
tPred <- seq(0, 1, length.out=100)

m <- stan_model("gpQuantileTrenderFixed.stan")

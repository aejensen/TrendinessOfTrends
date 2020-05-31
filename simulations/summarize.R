rm(list=ls())

library(genlasso)
library(parallel)
library(pracma)
library(rootSolve)

Seq <- function (a, b, ...) {
  if (a <= b) 
    return(seq(a, b, ...))
  else return(numeric(0))
}

doCV <- function (object, k = 5, mode = c("lambda", "df"), approx = FALSE, 
          rtol = 1e-07, btol = 1e-07, verbose = FALSE) {
  cl = match.call()
  if (all(class(object) != "trendfilter")) {
    stop("Cross-validation can only be performed for trend filtering.")
  }
  if (!is.null(object$X)) {
    stop("Cross-validation for trend filtering can only be performed when X=I, the identity matrix.")
  }
  mode = mode[[1]]
  if (!(mode %in% c("lambda", "df"))) {
    stop("Invalid mode, must be \"lambda\" or \"df\".")
  }
  y = object$y
  n = length(y)
  if (k < 2 || round(k) != k || k > n - 2) {
    stop("The number of folds must an integer between 2 and n-2.")
  }
  ord = object$ord
  pos = object$pos
  if (is.null(pos)) 
    pos = 1:n
  foldid = c(0, rep(Seq(1, k), n - 2)[Seq(1, n - 2)], 0)
  if (mode == "lambda") {
    lambda = object$lambda
    cvall = matrix(0, k, length(lambda))
    for (i in Seq(1, k)) {
      #cat(sprintf("Fold %i ... ", i))
      otr = which(foldid != i)
      ntr = length(otr)
      ytr = y[otr]
      ptr = pos[otr]
      Dtr = genlasso:::getDtfPosSparse(ntr, ord, ptr)
      out = genlasso:::dualpathWideSparse(ytr, Dtr, NULL, approx, 
                               Inf, min(lambda), rtol, btol, verbose)
      out$beta = as.matrix(ytr - t(Dtr) %*% out$u)
      out$fit = out$beta
      out$y = ytr
      out$bls = ytr
      b = coef.genlasso(out, lambda = lambda)$beta
      ote = which(foldid == i)
      yte = matrix(y[ote], length(ote), length(lambda))
      pte = pos[ote]
      ilo = which((Seq(1, n) %in% (ote - 1))[otr])
      ihi = which((Seq(1, n) %in% (ote + 1))[otr])
      a = (pte - ptr[ilo])/(ptr[ihi] - ptr[ilo])
      pred = b[ilo, ] * (1 - a) + b[ihi, ] * a
      cvall[i, ] = colMeans((yte - pred)^2)
    }
    cverr = colMeans(cvall)
    cvse = apply(cvall, 2, sd)/sqrt(k)
    names(cverr) = names(cvse) = round(lambda, 3)
    i0 = which.min(cverr)
    lam.min = lambda[i0]
    lam.1se = max(lambda[cverr <= cverr[i0] + cvse[i0]])
    i.min = which(lambda == lam.min)
    i.1se = which(lambda == lam.1se)
    out = list(err = cverr, se = cvse, mode = "lambda", lambda = lambda, 
               lambda.min = lam.min, lambda.1se = lam.1se, i.min = i.min, 
               i.1se = i.1se, call = cl)
  } else {
    stop("error")
  }
  class(out) = c("cv.trendfilter", "list")
  return(out)
}

getResults <- function(obj) {
  out <- do.call("rbind", mclapply(1:length(obj), function(i) {
    if(i %% 100 == 0) cat(i, "\n")
    q <- obj[[i]]
    
    #Do trend filtering
    m <- genlasso::trendfilter(q$y, ord=1) #linear trend filtering => derivtive is p.w. constant
    cv <- doCV(m, k = 10)
    tf <- approxfun(q$t, coef.genlasso(m, lambda=cv$lambda.min)$beta, rule = 2)
    d_tf <- function(x) pracma::fderiv(tf, x)
    
    #Calculate norms for f
    #GP model
    d_f <- pracma::trapz(q$tPred, q$f - q$fMu)
    d2_f <- pracma::trapz(q$tPred, (q$f - q$fMu)^2)
    
    #Trend filtering
    d_f_tf <- pracma::trapz(q$tPred, q$f - tf(q$tPred))
    d2_f_tf <- pracma::trapz(q$tPred, (q$f - tf(q$tPred))^2)
    
    #Calculate norms for df
    #GP model
    d_df <- pracma::trapz(q$tPred, q$df - q$dfMu)
    d2_df <- pracma::trapz(q$tPred, (q$df - q$dfMu)^2)
    
    #Trend filtering
    d_df_tf <- pracma::trapz(q$tPred, q$df - d_tf(q$tPred))
    d2_df_tf <- pracma::trapz(q$tPred, (q$df - d_tf(q$tPred))^2)
    
    #Calculate norms for TDI
    d_TDI <- pracma::trapz(q$tPred, as.numeric(q$df > 0) - q$TDI)
    d2_TDI <- pracma::trapz(q$tPred, (as.numeric(q$df > 0) - q$TDI)^2)
    
    #Calculate norms for ETI
    roots <- rootSolve::uniroot.all(approxfun(q$tPred, q$df), c(0,1))
    if(length(roots) == 0) { #special case with no crossings
      ETI_true <- function(t) 0*t 
    } else {
      #True counting process
      ETI_true <- approxfun(c(roots), 1:(length(roots)), 
                            yleft = 0, yright = length(roots), method="constant")
    }
    #Difference between true counting process and cumulative estimated dETI
    dif <- ETI_true(q$tPred) - pracma::cumtrapz(q$tPred, q$ETI)[,1]
    d_ETI <- pracma::trapz(q$tPred, dif)
    d2_ETI <- pracma::trapz(q$tPred, dif^2)
    
    ETI <- pracma::trapz(q$tPred, q$ETI)
    
    c(d_f, d_f_tf, d_df, d_df_tf,
      d2_f, d2_f_tf, d2_df, d2_df_tf,
      d_TDI, d2_TDI,
      d_ETI, d2_ETI, ETI)
  }, mc.cores=64))
  
  colnames(out) <- c("resid_f", "resid_f_tf", "resid_df", "resid_df_tf",
                     "L2_f", "L2_f_tf", "L2_df", "L2_df_tf",
                     "resid_TDI", "L2_TDI",
                     "resid_ETI", "L2_ETI", "ETI")
  out
}

load("sim25_res.RData")
sim25_1 <- sim25_1[!sapply(sim25_1, function(q) is.null(q$fMu))]
sim25_2 <- sim25_2[!sapply(sim25_2, function(q) is.null(q$fMu))]
sim25_3 <- sim25_3[!sapply(sim25_3, function(q) is.null(q$fMu))]
sim25_4 <- sim25_4[!sapply(sim25_4, function(q) is.null(q$fMu))]
sim25_5 <- sim25_5[!sapply(sim25_5, function(q) is.null(q$fMu))]

res25_1 <- getResults(sim25_1)
res25_2 <- getResults(sim25_2)
res25_3 <- getResults(sim25_3)
res25_4 <- getResults(sim25_4)
res25_5 <- getResults(sim25_5)

res25 <- list(res25_1, res25_2, res25_3, res25_4, res25_5)
save(res25, file="summary25.RData")

load("sim50_res.RData")
sim50_1 <- sim50_1[!sapply(sim50_1, function(q) is.null(q$fMu))]
sim50_2 <- sim50_2[!sapply(sim50_2, function(q) is.null(q$fMu))]
sim50_3 <- sim50_3[!sapply(sim50_3, function(q) is.null(q$fMu))]
sim50_4 <- sim50_4[!sapply(sim50_4, function(q) is.null(q$fMu))]
sim50_5 <- sim50_5[!sapply(sim50_5, function(q) is.null(q$fMu))]

res50_1 <- getResults(sim50_1)
res50_2 <- getResults(sim50_2)
res50_3 <- getResults(sim50_3)
res50_4 <- getResults(sim50_4)
res50_5 <- getResults(sim50_5)

res50 <- list(res50_1, res50_2, res50_3, res50_4, res50_5)
save(res50, file="summary50.RData")


load("sim100_res.RData")
sim100_1 <- sim100_1[!sapply(sim100_1, function(q) is.null(q$fMu))]
sim100_2 <- sim100_2[!sapply(sim100_2, function(q) is.null(q$fMu))]
sim100_3 <- sim100_3[!sapply(sim100_3, function(q) is.null(q$fMu))]
sim100_4 <- sim100_4[!sapply(sim100_4, function(q) is.null(q$fMu))]
sim100_5 <- sim100_5[!sapply(sim100_5, function(q) is.null(q$fMu))]

res100_1 <- getResults(sim100_1)
res100_2 <- getResults(sim100_2)
res100_3 <- getResults(sim100_3)
res100_4 <- getResults(sim100_4)
res100_5 <- getResults(sim100_5)

res100 <- list(res100_1, res100_2, res100_3, res100_4, res100_5)
save(res100, file="summary100.RData")

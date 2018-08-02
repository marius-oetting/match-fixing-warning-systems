# This file contains code for the bivariate Poisson distribution as a new families-object within gamboostLSS

# author: Marius Ötting
# email for correspondence: marius.oetting@uni-bielefeld

# Thanks to Andreas Mayr for his support with the implementation of bivariate distributions within gamboostLSS!


# gamboostLSS family for bivariate Poisson --------------------------------
# three parameters: lambda1, lambda2, lambda3


## sub-families 
# sub-Family for lambda1

PoissonMVlogLambda1 <- function(lambda1 = NULL, lambda2 = NULL, lambda3 = NULL, y2 = y2, scale.grad = TRUE) {
  # loss is negative log-Likelihood
  # f is the parameter to be fitted with log-link: exp(f) = lambda1 
  # y = y1
  loss <- function(lambda2 = lambda2, lambda3, y, y2, f = f) {
    res <- intm.res <- c()
    for (i in 1:length(y)) {
      kindex <- 0:min(y[i], y2[i])
      intm.res <- -((-exp(f[i]) - lambda2[i] - lambda3[i]) + (y[i] * f[i] - log(factorial(y[i]))) + 
                      (y2[i] * log(lambda2[i]) - log(factorial(y2[i]))) + log(sum((choose(y[i], kindex) * 
                                                                                     choose(y2[i], kindex) * factorial(kindex)) * (lambda3[i]/(exp(f[i]) * lambda2[i]))^kindex)))
      res <- c(res, intm.res)
    }
    res
  }
  
  # ngradient is the negative derivate w.r.t. f
  ngradient <- function(y, f, w = 1) {
    res <- intm.res <- c()
    for (i in 1:length(y)) {
      kindex <- 0:min(y[i], y2[i])
      intm.res <- y[i] - exp(f[i]) - (sum(kindex * factorial(kindex) * exp(-f[i] * kindex) * choose(y[i], kindex) * 
                                            choose(y2[i], kindex) * lambda3[i]^kindex/(lambda2[i]^kindex))/sum(exp(-f[i] * 
                                                                                                                     kindex) * factorial(kindex) * (choose(y[i], kindex) * choose(y2[i], kindex)) * (lambda3[i]/(lambda2[i]))^kindex))
      res <- c(res, intm.res)
    }
    ngr <- res
    
    if (scale.grad) {
      div <- sqrt(weighted.mean(ngr^2, w = w, na.rm = TRUE))
      div <- ifelse(div < 1e-04, 1e-04, div)
      ngr <- ngr/div
    }
    return(ngr)
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, lambda2 = lambda2, lambda3 = lambda3, y2 = y2))
  }
  
  
  offset <- function(y, w) {
    if (!is.null(lambda1)) {
      RET <- lambda1
    } else {
      RET <- log(weighted.mean(y, w = w, na.rm = TRUE))
    }
    return(RET)
  }
  
  # use the Family constructor of mboost
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f), offset = offset, 
                  name = "Bivariate Poisson distribution: lambda1 (log link)")
}


## sub-Family for lambda2

PoissonMVlogLambda2 <- function(lambda1 = NULL, lambda2 = NULL, lambda3 = NULL, y1 = y1, scale.grad = TRUE) {
  # loss is negative log-Likelihood
  # f is the parameter to be fitted with log-link: exp(f) = lambda2
  # y = y2
  loss <- function(lambda1 = lambda1, lambda3, y, y1, f = f) {
    res <- intm.res <- c()
    for (i in 1:length(y)) {
      mink <- min(y[i], y1[i])
      kindex <- 0:mink
      intm.res <- -((-lambda1[i] - exp(f[i]) - lambda3[i]) + (y1[i] * log(lambda1[i]) - log(factorial(y1[i]))) + 
                      (y[i] * f[i] - log(factorial(y[i]))) + log(sum((choose(y1[i], kindex) * choose(y[i], kindex) * 
                                                                        factorial(kindex)) * (lambda3[i]/(lambda1[i] * exp(f[i])))^kindex)))
      res <- c(res, intm.res)
    }
    res
  }
  
  # ngradient is the negative derivate of the loss w.r.t. f
  ngradient <- function(y, f, w = 1) {
    res <- intm.res <- c()
    for (i in 1:length(y)) {
      mink <- min(y[i], y1[i])
      kindex <- 0:mink
      intm.res <- y[i] - exp(f[i]) - (sum(kindex * exp(-f[i] * kindex) * factorial(kindex) * choose(y[i], kindex) * 
                                            choose(y1[i], kindex) * lambda3[i]^kindex/(lambda1[i]^kindex))/sum(exp(-f[i] * 
                                                                                                                     kindex) * factorial(kindex) * (choose(y[i], kindex) * choose(y1[i], kindex)) * (lambda3[i]/(lambda1[i]))^kindex))
      res <- c(res, intm.res)
    }
    ngr <- res
    if (scale.grad) {
      div <- sqrt(weighted.mean(ngr^2, w = w, na.rm = TRUE))
      div <- ifelse(div < 1e-04, 1e-04, div)
      ngr <- ngr/div
    }
    return(ngr)
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, lambda1 = lambda1, lambda3 = lambda3, y1 = y1))
  }
  
  
  
  offset <- function(y, w) {
    if (!is.null(lambda2)) {
      RET <- lambda2
    } else {
      RET <- log(weighted.mean(y, w = w, na.rm = TRUE))
    }
    return(RET)
  }
  # use the Family constructor of mboost
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f), offset = offset, 
                  name = "Bivariate Poisson distribution: lambda2 (log link)")
}

## sub-Family for lambda3

PoissonMVlogLambda3 <- function(lambda1 = NULL, lambda2 = NULL, lambda3 = NULL, y2 = y2, scale.grad = TRUE) {
  # loss is negative log-Likelihood
  # f is the parameter to be fitted with log-link: exp(f) = lambda3 
  # y = y1
  loss <- function(lambda1, lambda2, y, y2, f = f) {
    res <- intm.res <- c()
    for (i in 1:length(y)) {
      mink <- min(y[i], y2[i])
      kindex <- 0:mink
      intm.res <- -((-lambda1[i] - lambda2[i] - exp(f[i])) + (y[i] * log(lambda1[i]) - log(factorial(y[i]))) + 
                      (y2[i] * log(lambda2[i]) - log(factorial(y2[i]))) + log(sum((choose(y[i], kindex) * choose(y2[i], kindex) * 
                                                                                     factorial(kindex)) * (exp(f[i])/(lambda1[i] * lambda2[i]))^kindex)))
      res <- c(res, intm.res)
    }
    res
  }
  
  # ngradient is the negative derivate w.r.t. f
  ngradient <- function(y, f, w = 1) {
    res <- intm.res <- c()
    for (i in 1:length(y)) {
      mink <- min(y[i], y2[i])
      kindex <- 0:mink
      intm.res <- -exp(f[i]) + (sum(kindex * exp(f[i] * kindex) * factorial(kindex) * choose(y[i], kindex) 
                                    * choose(y2[i], kindex)/(lambda1[i]^kindex * lambda2[i]^kindex))/sum(exp(f[i] * 
                                                                                                               kindex) * factorial(kindex) * (choose(y[i], kindex) * choose(y2[i], kindex))/(lambda1[i]^kindex * 
                                                                                                                                                                                               lambda2[i]^kindex)))
      res <- c(res, intm.res)
    }
    ngr <- res
    
    if (scale.grad) {
      div <- sqrt(weighted.mean(ngr^2, w = w, na.rm = TRUE))
      div <- ifelse(div < 1e-04, 1e-04, div)
      ngr <- ngr/div
    }
    return(ngr)
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, lambda1 = lambda1, lambda2 = lambda2, y2 = y2))
  }
  
  
  
  offset <- function(y, w = 1) {
    if (!is.null(lambda3)) {
      RET <- lambda3
    } else {
      RET <- log(weighted.mean(y, w = w, na.rm = TRUE))
    }
    return(RET)
  }
  # use the Family constructor of mboost
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f), offset = offset, 
                  name = "Bivariate Poisson distribution: lambda3 (log link)")
}



# complete gamboostLSS families -------------------------------------------
PoissonMVlog <- function(lambda1 = NULL, lambda2 = NULL, lambda3 = NULL, y1 = NULL, y2 = NULL, scale.grad = TRUE) {
  if ((!is.null(lambda1) && lambda1 <= 0)) 
    stop(sQuote("lambda1"), " must be greater than zero")
  if ((!is.null(lambda2) && lambda2 <= 0)) 
    stop(sQuote("lambda2"), " must be greater than zero")
  if ((!is.null(lambda3) && lambda3 <= 0)) 
    stop(sQuote("lambda3"), " must be greater than zero")
  
  
  Families(lambda1 = PoissonMVlogLambda1(lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, y2 = y2, 
                                         scale.grad = scale.grad), lambda2 = PoissonMVlogLambda2(lambda1 = lambda1, lambda2 = lambda2, 
                                                                                                 lambda3 = lambda3, y1 = y1, scale.grad = scale.grad), lambda3 = PoissonMVlogLambda3(lambda1 = lambda1, 
                                                                                                                                                                                     lambda2 = lambda2, lambda3 = lambda3, y2 = y2, scale.grad = scale.grad), name = "PoissonMVlog")
}
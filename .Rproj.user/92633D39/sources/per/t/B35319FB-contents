library(hypergeo)
hypergeomFun <- function(x, n, lambda0, lambda1, lambda2){
  b <- n + lambda0
  a <- ifelse(x[1] < x[2], 1 - lambda2, 1 - lambda1)
  c <- ifelse(x[1] < x[2], n + lambda0 + lambda1, n + lambda0 + lambda2)
  z <- ifelse(x[1] < x[2], x[1]/x[2], x[2]/x[1])
  return(Re(hypergeo(a, b, c, z)))
}

SeriesComp1 <- function(x, n, lambda0, lambda1, lambda2){
  return(ifelse(x[1]<x[2], x[1]^(n+lambda0+lambda1-1)*x[2]^(lambda2-1),
                x[1]^(lambda1-1)*x[2]^(n+lambda0+lambda2-1)))
}

SeriesComp2 <- function(x, n, lambda0, lambda1, lambda2){
  return(ifelse(x[1]<x[2],
                gamma(n + lambda0)*gamma(lambda1)/gamma(n + lambda0 + lambda1),
                gamma(n + lambda0)*gamma(lambda2)/gamma(n + lambda0 + lambda2)))
}

SeriesComp <- function(x, n, lambda0, lambda1, lambda2){
  SeriesComp1(x, n, lambda0, lambda1, lambda2)*
    SeriesComp2(x, n, lambda0, lambda1, lambda2)*
    hypergeomFun(x, n, lambda0, lambda1, lambda2)/factorial(n)
}

SeriesSum <- function(x, n, lambda0, lambda1, lambda2){
  sum(sapply(0:n, function(i) SeriesComp(x, i, lambda0, lambda1, lambda2)))
}

Density2dim <- function(x, d, lambda0, lambda1, lambda2){
  const <- exp(-x[1]-x[2])/(gamma(lambda0)*gamma(lambda1)*gamma(lambda2))
  return(const*SeriesSum(x, d, lambda0, lambda1, lambda2))
}

MinusLogLike <- function(dataFunLike, lambda0, lambda1, lambda2){
  times <- 10
  sum(apply(dataFunLike, 2, function(x){
    -log(Density2dim(x, times, lambda0, lambda1, lambda2))
  }))
}

momentsFun <- function(x){
  m1 <- mean(x)
  m2 <- mean(x^2)
  return(m1^2/(m2-m1^2))
}

library(bbmle)
library("optimx")
estimML <- function(Y){
  Lambda1 <- momentsFun(Y[1,])
  Lambda2 <- momentsFun(Y[2,])
  l0 <- cor(Y[1,], Y[2,], method = "pearson", use = "pairwise.complete.obs")*sqrt(Lambda1*Lambda2)
  l1 <- Lambda1 - l0
  l2 <- Lambda2 - l0
  if(l1+l2 <= 0) return(NaN)
  est <- mle2(MinusLogLike,
              start = list(lambda0 = l0, lambda1 = l1, lambda2 = l2),
              lower = c(lambda0 = 0.1, lambda1 = 0.1, lambda2 = 0.1),
              upper = c(lambda0 = 50, lambda1 = 80, lambda2 = 80),
              data = list(dataFunLike = na.omit(Y)), method ="L-BFGS-B", optimizer = "optimx")
  return(est)
}
# вычисление гипергеометрической функции
library(hypergeo)
hypergeomFun <- function(x, n, lambda0, lambda1, lambda2){
  b <- n + lambda0
  a <- ifelse(x[1] < x[2], 1 - lambda2, 1 - lambda1)
  c <- ifelse(x[1] < x[2], n + lambda0 + lambda1, n + lambda0 + lambda2)
  z <- ifelse(x[1] < x[2], x[1]/x[2], x[2]/x[1])
  return(Re(hypergeo(a, b, c, z)))
}

# сомножители компоненты ряда
SeriesComp1 <- function(x, n, lambda0, lambda1, lambda2){
  return(ifelse(x[1]<x[2], x[1]^(n+lambda0+lambda1-1)*x[2]^(lambda2-1),
                x[1]^(lambda1-1)*x[2]^(n+lambda0+lambda2-1)))
}

SeriesComp2 <- function(x, n, lambda0, lambda1, lambda2){
  return(ifelse(x[1]<x[2],
                gamma(n + lambda0)*gamma(lambda1)/gamma(n + lambda0 + lambda1),
                gamma(n + lambda0)*gamma(lambda2)/gamma(n + lambda0 + lambda2)))
}

# компоненту ряда получаем перемножением сомноителей
SeriesComp <- function(x, n, lambda0, lambda1, lambda2){
  SeriesComp1(x, n, lambda0, lambda1, lambda2)*
    SeriesComp2(x, n, lambda0, lambda1, lambda2)*
    hypergeomFun(x, n, lambda0, lambda1, lambda2)/factorial(n)
}

# частичная сумма ряда
# вычисляем компоненты до тех пор, пока их значение не окажется достаточно малым
SeriesSum <- function(x, n, lambda0, lambda1, lambda2){
  sum(sapply(0:n, function(i) SeriesComp(x, i, lambda0, lambda1, lambda2)))
}

# плотность распределения получаем из суммы ряда умножением на соответсвующую константу
Density2dim <- function(x, d, lambda0, lambda1, lambda2){
  const <- exp(-x[1]-x[2])/(gamma(lambda0)*gamma(lambda1)*gamma(lambda2))
  return(const*SeriesSum(x, d, lambda0, lambda1, lambda2))
}

# для использования в функции mle2 вычислим минус логарифм функции правдоподобия (ФП)
MinusLogLike <- function(dataFunLike, lambda0, lambda1, lambda2){
  times <- 10
  sum(apply(dataFunLike, 2, function(x){
    -log(Density2dim(x, times, lambda0, lambda1, lambda2))
  }))
}

# оценка парметров формы по методу моментов
momentsFun <- function(x){
  m1 <- mean(x)
  m2 <- mean(x^2)
  return(m1^2/(m2-m1^2))
}

# определим начальное приближение из связи параметров распределения
# с параметрами формы, полученными по методу моментов
# параметры распределения должны быть положительными
# нижняя граница отделена от нуля для того, чтобы была определена ФП
library(bbmle)
library("optimx")
estimML <- function(Y){
  Lambda1 <- momentsFun(Y[1,])
  Lambda2 <- momentsFun(Y[2,])
  l0 <- cor(Y[1,], Y[2,], method = "pearson", use = "pairwise.complete.obs")*sqrt(Lambda1*Lambda2)
  l1 <- Lambda1 - l0
  l2 <- Lambda2 - l0
  # условие на сходимость ряда
  if(l1+l2 <= 0) return(NaN)
  est <- mle2(MinusLogLike,
              start = list(lambda0 = l0, lambda1 = l1, lambda2 = l2),
              lower = c(lambda0 = 0.01, lambda1 = 0.01, lambda2 = 0.01),
              upper = c(lambda0 = 50, lambda1 = 80, lambda2 = 80),
              data = list(dataFunLike = Y), method ="L-BFGS-B", optimizer = "optimx")
  return(est)
}

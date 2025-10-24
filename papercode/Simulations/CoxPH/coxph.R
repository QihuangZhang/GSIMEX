library(survival)

n = 10000

X = rnorm(n,0,1)

U = runif(n,0,1)


T = sqrt(-1*exp(-X) * log(1-U))

C = runif(n,0,max(T))

Y = (T<C) * T + (T>C) * C
delta = (T<C)*1


  fit <- coxph(Surv(Y, delta) ~ X)

fit

#  predict(fit)


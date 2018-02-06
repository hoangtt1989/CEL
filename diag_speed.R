pacman::p_load(microbenchmark, Matrix)
n <- 1000
p <- 5

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
wts <- rnorm(n)
y <- rnorm(n)

f1 <- function(X, wts) {
  R <- diag(wts)
  R %*% X
}
f2 <- function(X, wts) {
  R <- Diagonal(x = wts)
  res <- R %*% X
  as.matrix(res)
}
f3 <- function(X, wts) {
  (wts * X)
}

sum(f1(X, wts) - f3(X, wts))

R <- diag(wts)
tst1 <- t(X) %*% R
tst2 <- t(wts * X)
sum(tst1 - tst2)

sum(R %*% y - wts * y)

microbenchmark(f1(X, wts), f2(X, wts), f3(X, wts))

microbenchmark(wts %*% y, crossprod(wts, y))

microbenchmark(crossprod(tst1, tst1 %*% wts), crossprod(tst1) %*% wts)

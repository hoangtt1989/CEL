linesch_ww2 <- function (fn, gr, x0, d, fun_arg, grad_arg, c1 = 10e-4, 
                         c2 = 0.5, fvalquit = -Inf, prtlevel = 0) 
{
  stopifnot(is.numeric(x0), is.numeric(d))
  if (c1 < 0 || c1 > c2 || c2 > 1) 
    stop("Arguments 'c1','c2' must satisfy: 0 <= c1 <= c2 <= 1/n.")
  n <- length(x0)
  alpha <- 0
  xalpha <- x0
  fn0 <- do.call('fn', c(list(x0), fun_arg), quote = T)
  gr0 <- do.call('gr', c(list(x0), grad_arg), quote = T)
  falpha <- fn0
  galpha <- gr0
  beta <- Inf
  gbeta <- rep(NA, n)
  g0 <- sum(gr0 * d)
  if (g0 >= 0) 
    if (prtlevel > 0) 
      warning(paste("Linesearch: Argument 'd' is not a descent direction."))
  dnorm <- sqrt(sum(d * d))
  if (dnorm == 0) 
    stop("Linesearch: Argument 'd' must have length greater zero.")
  t <- 1
  nfeval <- 0
  nbisect <- 0
  nexpand <- 0
  nbisectmax <- max(30, round(log2(1e+05 * dnorm)))
  nexpandmax <- max(10, round(log2(1e+05/dnorm)))
  fevalrec <- c()
  fail <- 0
  done <- FALSE
  while (!done) {
    x <- x0 + t * d
    fun <- do.call('fn', c(list(x), fun_arg), quote = T)
    grd <- do.call('gr', c(list(x), grad_arg), quote = T)
    nfeval <- nfeval + 1
    fevalrec <- c(fevalrec, fun)
    if (fun < fvalquit) {
      return(list(alpha = t, xalpha = x, falpha = fun, 
                  galpha = grd, fail = fail, beta = beta, gbeta = gbeta, 
                  fevalrec = fevalrec))
    }
    gtd <- sum(grd * d)
    if (fun >= fn0 + c1 * t * g0 || is.na(fun)) {
      beta <- t
      gbeta <- grd
    }
    else if (gtd <= c2 * g0 || is.na(gtd)) {
      alpha <- t
      xalpha <- x
      falpha <- fun
      galpha <- grd
    }
    else {
      return(list(alpha = t, xalpha = x, falpha = fun, 
                  galpha = grd, fail = fail, beta = t, gbeta = grd, 
                  fevalrec = fevalrec))
    }
    if (beta < Inf) {
      if (nbisect < nbisectmax) {
        nbisect <- nbisect + 1
        t <- (alpha + beta)/2
      }
      else {
        done <- TRUE
      }
    }
    else {
      if (nexpand < nexpandmax) {
        nexpand <- nexpand + 1
        t <- 2 * alpha
      }
      else {
        done <- TRUE
      }
    }
  }
  if (is.infinite(beta)) {
    fail <- -1
    if (prtlevel > 0) 
      warning(paste("Linesearch: Function may be unbounded from below."))
  }
  else {
    fail <- 1
    if (prtlevel > 0)
      warning(paste("Linesearch: Failed to satisfy weak Wolfe conditions."))
  }
  list(alpha = t, xalpha = x, falpha = fun, galpha = grd, fail = fail, 
       beta = t, gbeta = grd, fevalrec = fevalrec)
}


lbfgs_step <- function(x_curr, hess_init, mem_rho, mem_s, mem_y, iter, fun, grad, fun_arg = list(), grad_arg = list(),
                       m = 6, ls_maxit = 100, ls_eps = 1e-12, ls_beta = .5, ls_alpha = .25) {
  # fun_arg$llog_eps <- 50
  # fun_arg$rev_sign <- T
  # grad_arg$rev_sign <- T
  fun_curr <- do.call('fun', c(list(x_curr), fun_arg), quote = T)
  grad_new <- do.call('grad', c(list(x_curr), grad_arg), quote = T)
  alpha_vals <- rep(0, m)
  # first_iter <- ifelse(iter > m, (iter - 1):(iter - m), iter)
  first_iter <- ifelse(iter > m, m:1, m:(m - iter + 1))
  for (i in first_iter) {
    alpha_vals[i] <- mem_rho[i] * (mem_s[, i] %*% grad_new)
    grad_new <- grad_new - alpha_vals[i] * mem_y[, i]
  }
  r_new <- hess_init %*% grad_new
  # second_iter <- ifelse(iter > m, (iter - m):(iter - 1), iter)
  # second_iter <- ifelse(iter > m, 1:m, (m - iter + 1):m)
  second_iter <- rev(first_iter)
  for (i in second_iter) {
    beta_new <- mem_rho[i] * mem_y[, i] %*% r_new
    r_new <- r_new + mem_s[, i] * (alpha_vals[i] - beta_new)
  }
  search_dir <- -r_new
  x_new <- x_curr + search_dir
  gradtr_pk <- crossprod(grad_new, search_dir)
  #########line search - wolfe conditions
  res <- linesch_ww2(fun, grad, x_new, search_dir, fun_arg, grad_arg, prtlevel = 1)
  x_new <- res$xalpha
  ls_converged <- res$fail == 0
  # x_new <- wolfe.linesearch(fun, as.numeric(x_new), as.numeric(search_dir), fun_arg)
  # ls_fval <- do.call('fun', c(list(x_new), fun_arg), quote = T)
  # step <- 1
  # ls_converged <- ls_fval < (fun_curr + ls_alpha * step * gradtr_pk)
  # if (!ls_converged) {
  #   for (ls_iter in 1:ls_maxit) {
  #     step <- ls_beta * step
  #     x_new <- x_curr + step * search_dir
  #     ls_fval <- do.call('fun', c(list(x_new), fun_arg), quote = T)
  #     ls_diff <- fun_curr + ls_alpha * step * gradtr_pk - ls_fval
  #     # print(ls_diff)
  #     if ((ls_diff >= 0 | abs(ls_diff) < ls_eps)) {
  #       ls_converged <- T
  #       break
  #     }
  #   }
  # }
  #########
  return(list(x_new = x_new, fval_curr = fun_curr, gval_new = grad_new, ls_converged = ls_converged))
}

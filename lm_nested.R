library(pracma) ##for orth, nullspace
library(compiler) ##for cmpfun
source('owenEL_old.R')
source('damped_newton.R')
source('lbfgs_step.R')

##################functions for hessian/gradient calculations
lm_fval_z <- function(z_reparam, X, F_reparam, s_hat_reparam, y, gamma, llog_eps = 1 / nrow(X), rev_sign = F) {
  n <- nrow(X)
  beta_new <- F_reparam %*% z_reparam + s_hat_reparam
  R_new <- as.vector(y - X %*% beta_new)
  score_new <- R_new * X
  ret <- sum(llog(1 + score_new %*% gamma, llog_eps))
  if (rev_sign) {
    ret <- -ret
  }
  return(ret)
}
lm_hess_z <- function(z_reparam, X, F_reparam, s_hat_reparam, y, gamma, llog_eps = 1 / nrow(X)) {
  ##L_z_gamma
  n <- nrow(X)
  p <- ncol(F_reparam)
  score_eq <- as.vector(y - X %*% (F_reparam %*% z_reparam + s_hat_reparam))
  score_mat <- score_eq * X
  inner_sum <- 1 + score_eq * (X %*% gamma)
  temp_plog1 <- llogp(inner_sum, llog_eps)
  temp_plog2 <- llogpp(inner_sum, llog_eps)
  temp1 <- matrix(0, nrow = p, ncol = ncol(X))
  temp2 <- temp1
  for (i in 1:n) {
    temp1 <- temp1 + temp_plog1[i] * crossprod(F_reparam, tcrossprod(X[i, ]))
    temp2 <- temp2 + temp_plog2[i] * crossprod(F_reparam, X[i, ] * as.numeric(X[i, ] %*% gamma)) %*% t(X[i, ]) * (y[i] - as.numeric(t(X[i, ]) %*% (F_reparam %*% z_reparam + s_hat_reparam)))
  }
  L_z_gamma <- -temp1 - temp2
  ##
  L_gamma_z <- t(L_z_gamma)
  L_gamma_gamma <- crossprod(score_mat, as.vector(llogpp(1 + score_mat %*% gamma, llog_eps)) * score_mat)
  ##L_z_z
  temp_res <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
    temp_res <- temp_res + temp_plog2[i] * (crossprod(F_reparam, X[i, ]) * as.numeric(gamma %*% X[i, ])^2) %*% crossprod(X[i, ], F_reparam)
  }
  L_z_z <- temp_res
  ##
  hess <- L_z_z - L_z_gamma %*% solve(L_gamma_gamma, L_gamma_z)
}
lm_part_z <- function(z_reparam, X, F_reparam, s_hat_reparam, y, gamma, llog_eps = 1 / nrow(X), rev_sign = F) {
  n <- nrow(X)
  p <- ncol(F_reparam)
  score_eq <- as.vector(y - X %*% (F_reparam %*% z_reparam + s_hat_reparam))
  inner_sum <- 1 + score_eq * (X %*% gamma)
  temp_plog1 <- llogp(inner_sum, llog_eps)
  temp_res <- as.vector(rep(0, p))
  for (i in 1:n) {
    temp_res <- temp_res + temp_plog1[i] * crossprod(F_reparam, X[i, ] * as.numeric(X[i, ] %*% gamma))
  }
  ret <- -temp_res
  if (rev_sign) {
    ret <- -ret
  }
  return(ret)
}
lm_part_z_gamma <- function(z_reparam, X, F_reparam, s_hat_reparam, y, gamma, llog_eps = 1) {
  n <- nrow(X)
  p <- ncol(F_reparam)
  score_eq <- as.vector(y - X %*% (F_reparam %*% z_reparam + s_hat_reparam))
  inner_sum <- 1 + score_eq * (X %*% gamma)
  temp_plog1 <- llogp(inner_sum, llog_eps)
  temp_plog2 <- llogpp(inner_sum, llog_eps)
  temp1 <- matrix(0, nrow = p, ncol = ncol(X))
  temp2 <- temp1
  for (i in 1:n) {
    temp1 <- temp1 + temp_plog1[i] * crossprod(F_reparam, tcrossprod(X[i, ]))
    temp2 <- temp2 + temp_plog2[i] * crossprod(F_reparam, X[i, ] * as.numeric(X[i, ] %*% gamma)) %*% t(X[i, ]) * (y[i] - as.numeric(t(X[i, ]) %*% (F_reparam %*% z_reparam + s_hat_reparam)))
  }
  -temp1 - temp2
}
lm_part_z_z <- function(z_reparam, X, F_reparam, s_hat_reparam, y, gamma, llog_eps = 1 / nrow(X)) {
  n <- nrow(X)
  p <- ncol(F_reparam)
  score_eq <- as.vector(y - X %*% (F_reparam %*% z_reparam + s_hat_reparam))
  temp_plog2 <- llogpp(1 + score_eq * (X %*% gamma), llog_eps)
  temp_res <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
    temp_res <- temp_res + temp_plog2[i] * (crossprod(F_reparam, X[i, ]) * as.numeric(gamma %*% X[i, ])^2) %*% crossprod(X[i, ], F_reparam)
  }
  temp_res
}
lm_fval_gamma <- function(gamma, score_mat, llog_eps = 1 / nrow(score_mat)) {
  -sum(llog(1 + score_mat %*% gamma, llog_eps))
}
lm_part_gamma <- function(gamma, score_mat, llog_eps = 1 / nrow(score_mat)) {
  -crossprod(score_mat, llogp(1 + score_mat %*% gamma, llog_eps))
}
##################

##################function for CEL linear regression nested - trying to write custom lbfgs code
CEL_lm_nested <- function(beta_test, X, y, F_reparam, s_hat_reparam,
                          outer_eps = 1e-8, outer_maxit = 1000, inner_opt_type = c('owen', 'LBFGS'), inner_owen_arg = list(), inner_lbfgs_arg = list(invisible = 1),
                          beta_opt_type = c('damped_newton', 'LBFGS'), beta_newton_arg = list(), beta_lbfgs_arg = list(m = 6),
                          outer_tol_type = c('fval', 'gval'), verbose = F) {
  ##checking inputs
  beta_opt_type <- beta_opt_type[1]
  outer_tol_type <- outer_tol_type[1]
  inner_opt_type <- inner_opt_type[1]
  if (!(outer_tol_type %in% c('fval', 'gval'))) {
    stop('supply a valid outer tolerance type')
  }
  ##
  ##compile owen's functions
  owen_EL <- cmpfun(elm) ##inner loop
  # plog0 <- cmpfun(llog) ##pseudo log
  # plog1 <- cmpfun(llogp) ##pseudo log first deriv
  # plog2 <- cmpfun(llogpp) ##pseudo log second deriv
  # owen_logelr <- cmpfun(owen_logelr_fun) ##logelr
  ##
  ##compile gradient/hessian functions
  lm_hess_z <- cmpfun(lm_hess_z)
  lm_part_z <- cmpfun(lm_part_z)
  lm_fval_z <- cmpfun(lm_fval_z)
  lm_fval_gamma <- cmpfun(lm_fval_gamma)
  lm_part_gamma <- cmpfun(lm_part_gamma)
  ##
  ##initial parameters
  prob_size <- dim(X)
  n <- prob_size[1]
  p <- prob_size[2]
  beta_curr <- beta_test
  z_reparam_curr <- crossprod(F_reparam, beta_test - s_hat_reparam)
  R_curr <- as.vector(y - X %*% beta_test)
  outer_converged <- F
  ##
  ##initialize holders
  outer_fval <- rep(NA, outer_maxit)
  inner_nits <- rep(NA, outer_maxit)
  ##
  ##first inner loop iteration
  score_curr <- R_curr * X
  if (inner_opt_type == 'owen') {
    inner_res <- do.call('owen_EL', c(list(score_curr), inner_owen_arg), quote = T)
    gamma_curr <- inner_res$lambda
    logelr_curr <- sum(llog(inner_res$wts, 1/n))
  } else if (inner_opt_type == 'LBFGS') {
    gamma_curr <- rep(0, p)
    inner_res <- do.call('lbfgs', c(list(lm_fval_gamma, lm_part_gamma, gamma_curr, score_mat = score_curr), inner_lbfgs_arg), quote = T)
    gamma_curr <- inner_res$par
    logelr_curr <- -sum(llog(1 + score_curr %*% gamma_curr, 1/n))
  }
  ##
  ##set up function, gradient, hessian arguments
  fun_arg <- list(X = X, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, y = y, gamma = gamma_curr)
  grad_arg <- list(X = X, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, y = y, gamma = gamma_curr)
  hess_arg <- list(X = X, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, y = y, gamma = gamma_curr)
  ##
  ##initial parameters for lbfgs
  if (beta_opt_type == 'LBFGS') {
    m <- lbfgs_arg$m
    num_par <- length(z_reparam_curr)
    mem_y <- matrix(0, nrow = num_par, ncol = m)
    mem_s <- matrix(0, nrow = num_par, ncol = m)
    mem_rho <- rep(0, m)
    # mem_y[, m] <- do.call('lm_part_z', c(list(z_reparam_curr), grad_arg), quote = T)
    # mem_s[, m] <- z_reparam_curr
    # mem_rho[m] <- 1 / as.numeric(mem_y[, m] %*% mem_s[, m])
    mem_y[, m] <- rep(0, num_par)
    mem_s[, m] <- rep(0, num_par)
    mem_rho[m] <- 0
    grad_curr <- mem_y[, 1]
    x_curr <- mem_s[, 1]
    # hess_init <- diag(num_par) * as.numeric(mem_s[, m] %*% mem_y[, m]) / as.numeric(mem_y[, m] %*% mem_y[, m])
    hess_init <- diag(num_par)
  }
  ##
  start_time <- Sys.time()
  for (j in 1:outer_maxit) {
    # print(j)
    if (verbose) {
      if (! j %% 5) {
        message(paste('Iteration', j))
      }
    }
    ##beta optimization
    if (beta_opt_type == 'damped_newton') {
      damped_res <- do.call('damped_newton_step', c(list(z_reparam_curr, lm_fval_z, lm_part_z, lm_hess_z, fun_arg, grad_arg, hess_arg), beta_newton_arg), quote = T)
      if (!damped_res$ls_converged) {
        warning(paste('Line search unsuccessful at j =', j))
      }
      x_curr <- damped_res$x_curr
      fval_curr <- damped_res$fval_curr
    } else if (beta_opt_type == 'LBFGS') {
      lbfgs_res <- do.call('lbfgs_step', c(list(z_reparam_curr, hess_init, mem_rho, mem_s, mem_y, j, lm_fval_z, lm_part_z, fun_arg, grad_arg), beta_lbfgs_arg), quote = T)
      if (!lbfgs_res$ls_converged) {
        warning(paste('Line search unsuccessful at j =', j))
      }
      mem_y[, 1:(m - 1)] <- mem_y[, 2:m]
      mem_s[, 1:(m - 1)] <- mem_s[, 2:m]
      mem_rho[1:(m - 1)] <- mem_rho[2:m]
      grad_new <- lbfgs_res$gval_new
      x_new <- lbfgs_res$x_new
      mem_y[, m] <- grad_new - grad_curr
      mem_s[, m] <- x_new - x_curr
      mem_rho[m] <- 1 / as.numeric(mem_y[, m] %*% mem_s[, m])
      grad_curr <- grad_new
      x_curr <- x_new
      fval_curr <- lbfgs_res$fval_curr
      hess_init <- diag(num_par) * as.numeric(mem_s[, m] %*% mem_y[, m]) / as.numeric(mem_y[, m] %*% mem_y[, m])
    }
    ##
    ##update values
    z_reparam_curr <- x_curr
    beta_curr <- F_reparam %*% z_reparam_curr + s_hat_reparam
    R_new <- as.vector(y - X %*% beta_curr)
    score_new <- R_new * X
    outer_fval_curr <- fval_curr
    # outer_fval_new <- damped_res$fval_new
    ##
    ####inner loop
    init_gamma <- gamma_curr
    init_gamma <- rep(0, p)
    if (inner_opt_type == 'owen') {
      inner_res <- do.call('owen_EL', c(list(score_new, lam = init_gamma), inner_owen_arg), quote = T)
      gamma_curr <- inner_res$lambda
      inner_nits[j] <- inner_res$nits
      outer_fval_new <- -sum(llog(inner_res$wts, 1/n))
    } else if (inner_opt_type == 'LBFGS') {
      inner_res <- do.call('lbfgs', c(list(lm_fval_gamma, lm_part_gamma, init_gamma, score_mat = score_new), inner_lbfgs_arg), quote = T)
      gamma_curr <- inner_res$par
      inner_nits[j] <- inner_res$convergence
      outer_fval_new <- sum(llog(1 + score_new %*% gamma_curr, 1/n))
    }
    outer_fval[j] <- outer_fval_new
    fun_arg$gamma <- gamma_curr
    grad_arg$gamma <- gamma_curr
    hess_arg$gamma <- gamma_curr
    ####
    #####stopping
    outer_tol <- switch(outer_tol_type,
                        fval = abs(outer_fval_new - outer_fval_curr)/abs(outer_fval_curr),
                        gval = sum(damped_res$gval_new^2))
    if (all(outer_tol < outer_eps) & !any(is.nan(outer_tol) || is.na(outer_tol) || is.infinite(outer_tol))) {
      outer_converged <- T
      break
    }
    #####
  }
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  if (!outer_converged) {
    warning(paste('Outer loop did not converge, current tolerance is', outer_tol))
  }

  ####cleaning up output
  outer_fval <- outer_fval[1:j]
  inner_nits <- inner_nits[1:j]
  logelr <- -outer_fval_new
  wts <- 1 / (1 + score_new %*% gamma_curr) / n
  sum_wts <- sum(wts)
  if (abs(sum_wts - 1) > 1e-6) { ##output warning if weights don't sum to 1
    warning('weights do not sum to 1 - convex hull constraint may not be satisfied')
  }
  ####

  return(list(logelr = logelr, logelr_stat = -2 * logelr, wts = wts, sum_wts = sum_wts,
              gamma = gamma_curr, beta_opt = beta_curr, outer_converged = outer_converged, time = tot_time,
              outer_fval = outer_fval, outer_tol = outer_tol, inner_nits = inner_nits, outer_nits = j))
}
##################

##function value for linear regression augmented lagrangian
# fun_lm_AL <- function(wts, gammatr_Xtr_R, Xtr_R, dual_step) {
#   -sum(log(wts)) + gammatr_Xtr_R %*% wts + .5 * dual_step * sum((Xtr_R %*% wts)^2)
# }
fun_lm_AL2 <- function(wts, gamma, X, R, dual_step) {
  Xtr_R <- t(R * X)
  -sum(log(wts)) + as.vector(gamma) %*% (Xtr_R %*% wts) + .5 * dual_step * sum((Xtr_R %*% wts)^2)
}
##
##gradient for linear regression augmented lagrangian
# grad_lm_AL <- function(wts, R_X_gamma, Xtr_R, dual_step) {
#   # -1 / wts + R_X_gamma + dual_step * t(Xtr_R) %*% (Xtr_R %*% wts)
#   -1 / wts + R_X_gamma + dual_step * crossprod(Xtr_R, Xtr_R %*% wts)
# }
grad_lm_AL2 <- function(wts, gamma, X, R, dual_step) {
  Xtr_R <- t(R * X)
  -1 / wts + crossprod(Xtr_R, gamma) + dual_step * crossprod(Xtr_R, Xtr_R %*% wts)
}
##
##bregman divergence
bregman_entropy <- function(x_new, x_curr) {
  crossprod(x_new, log(x_new)) - crossprod(x_curr, log(x_curr)) - crossprod((1 + log(x_curr)), x_new - x_curr)
}
##

mirror_descent <- function(fun, grad, wts, fun_arg = list(), grad_arg = list(), mirror_step = 1,
                           accelerate = T, maxit = 1000, vary_step = F, line_search = F, ls_eps = 1e-10, ls_maxit = 1000, ls_beta = .5,
                           tol_type = c('fval', 'logelr', 'wts'), mirror_eps = 1e-8, verbose = F) {
  ##checking inputs
  tol_type <- tol_type[1]
  if (!(tol_type %in% c('fval', 'logelr', 'wts'))) {
    stop('supply valid tolerance type for mirror descent')
  }
  ##
  ##checking function and gradient arguments - make sure wts is not already in the list
  if (!is.null(fun_arg[['wts']])) {
    fun_arg[['wts']] <- NULL
  }
  if (!is.null(grad_arg[['wts']])) {
    grad_arg[['wts']] <- NULL
  }
  ##
  ##initial parameters for acceleration
  if (accelerate) {
    x_curr <- wts
    y_curr <- wts
    nu_curr <- wts
  }
  ##
  ##initial parameters
  n <- length(wts)
  wts_curr <- wts
  converged <- F
  bad_break_flag <- F
  if (line_search) {
    mirror_step <- 1
  }
  if (vary_step && !line_search) {
    mirror_step <- sqrt(2 * log(n)) * mirror_step ## use the varying scheme from beck and teboulle
  }
  tiny <- .Machine$double.xmin ## a small number so there is no division by zero
  tiny2 <- tiny * 2
  ##
  ##initial values for tol
  if (tol_type == 'fval' || line_search) {
    fval_curr <- do.call('fun', c(list(wts = wts_curr), fun_arg), quote = T)
  }
  if (tol_type == 'logelr') {
    logelr_curr <- sum(log(wts_curr))
  }
  ##
  for (i in 1:maxit) {
    
    if (accelerate) {
      ######nesterov's 2nd acceleration
      theta_curr <- 2 / (i + 1)
      theta_new <- 2 / (i + 2)
      grad_update <- do.call('grad', c(list(wts = y_curr), grad_arg), quote = T)
      # grad_step <- mirror_step / theta_curr ##update stepsize
      if (vary_step && !line_search) {
        mirror_step <- mirror_step * 1 / sqrt(i)
      }
      nu_new <- nu_curr * exp(-mirror_step * grad_update / theta_curr)
      if (!line_search & any(is.nan(nu_new) || is.infinite(nu_new))) {
        bad_break_flag <- T
        tol <- Inf
        break
      }
      nu_new <- n * nu_new / sum(nu_new)
      x_new <- (1 - theta_curr) * x_curr + theta_curr * nu_new
      y_new <- (1 - theta_new) * x_new + theta_new * nu_new
      wts_new <- x_new
      ######line search
      if (line_search) {
        fval_new <- do.call('fun', c(list(wts = wts_new), fun_arg), quote = T)
        fval_y <- do.call('fun', c(list(wts = y_curr), fun_arg), quote = T)
        surr_curr <- fval_y + crossprod(grad_update, wts_new - y_curr) + (1/mirror_step) * bregman_entropy(wts_new, y_curr)
        ls_converged <- fval_new < surr_curr
        if (!ls_converged || any(is.nan(ls_converged) || is.infinite(ls_converged) || is.na(ls_converged))) {
          for (iter_ls in 1:ls_maxit) {
            mirror_step <- ls_beta * mirror_step
            nu_new <- nu_curr * exp(-mirror_step * grad_update / theta_curr)
            nu_new <- n * nu_new / sum(nu_new)
            x_new <- (1 - theta_curr) * x_curr + theta_curr * nu_new
            y_new <- (1 - theta_new) * x_new + theta_new * nu_new
            wts_new <- x_new
            fval_new <- do.call('fun', c(list(wts = wts_new), fun_arg), quote = T)
            surr_curr <- fval_y + crossprod(grad_update, wts_new - y_curr) + (1/mirror_step) * bregman_entropy(wts_new, y_curr)
            ls_diff <- surr_curr - fval_new
            if (!any(is.na(ls_diff) || is.nan(ls_diff) || is.infinite(ls_diff))) {
              if (ls_diff >= 0 || abs(ls_diff) < ls_eps) {
                ls_converged <- T
                break
              }
            }
          }
        }
      }
      ######
      ##update current values
      x_curr <- wts_new
      y_curr <- y_new
      nu_curr <- nu_new
      ##
    } else {
      grad_update <- do.call('grad', c(list(wts = wts_curr), grad_arg), quote = T)
      wts_new <- wts_curr * exp(-mirror_step * grad_update)
      if (any(is.nan(wts_new) || is.infinite(wts_new))) {
        bad_break_flag <- T
        tol <- Inf
        break
      }
      wts_new <- n * wts_new / sum(wts_new)
    }
    ##check tolerance
    if (tol_type == 'fval') {
      if (!line_search) {
        fval_new <- do.call('fun', c(list(wts = wts_new), fun_arg), quote = T)
      }
      fval_sc <- abs(fval_curr)
      fval_sc <- ifelse(fval_sc < tiny, tiny2, fval_sc)
      tol <- abs(fval_new - fval_curr) / fval_sc
      fval_curr <- fval_new
    } else if (tol_type == 'wts') {
      wts_curr_sc <- max(abs(wts_curr))
      wts_curr_sc <- ifelse(wts_curr_sc < tiny, tiny2, wts_curr_sc)
      tol <- max(abs(wts_new - wts_curr)) / wts_curr_sc
    } else if (tol_type == 'logelr') {
      logelr_new <- sum(log(wts_new))
      logelr_curr_sc <- abs(logelr_curr)
      logelr_curr_sc <- ifelse(logelr_curr_sc < tiny, tiny2, logelr_curr_sc)
      tol <- abs(logelr_new - logelr_curr) / logelr_curr_sc
      logelr_curr <- logelr_new
    }
    ##
    ##update weights, exit if tol
    wts_curr <- wts_new
    if (all(tol < mirror_eps) && !any(is.nan(tol) || is.na(tol) || is.infinite(tol))) {
      converged <- T
      break ##exit loop if tolerance met
    }
    ##
  }
  return(list(wts = wts_curr, iter = i, tol = tol, converged = converged, bad_break_flag = bad_break_flag))
}

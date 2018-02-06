###some prox operators
prox_affine <- function(inp, T_reparam, alpha_reparam) {
  inp - crossprod(T_reparam, solve(tcrossprod(T_reparam), T_reparam %*% inp - alpha_reparam))
}
prox_inequality <- function(inp, test_idx, test_val, dir = c('gt', 'lt')) {
  dir <- dir[1]
  if (length(test_val) != length(test_idx)) {
    test_val <- rep(test_val[1], length(test_idx))
  }
  cond_check <- inp[test_idx] >= test_val[order(test_idx)]
  if (dir == 'lt') {
    cond_check <- !cond_check
  }
  new_vals <- test_val
  new_vals[cond_check] <- inp[test_idx[cond_check]]
  ret <- inp
  ret[test_idx] <- new_vals
  return(ret)
}
prox_interval <- function(inp, test_idx, left_val, right_val) {
  if (!all(left_val <= right_val)) {
    stop('left_val must be <= right_val')
  }
  if (length(left_val) != length(right_val)) {
    stop('left_val must have same length as right_val')
  }
  if (length(left_val) != length(test_idx)) {
    left_val <- rep(left_val[1], length(test_idx))
    right_val <- rep(right_val[1], length(test_idx))
  }
  order_idx <- order(test_idx)
  left_check <- inp[test_idx] < left_val[order_idx]
  right_check <- inp[test_idx] > right_val[order_idx]
  ret <- inp
  ret[test_idx[left_check]] <- left_val[left_check]
  ret[test_idx[right_check]] <- right_val[right_check]
  return(ret)
}
###

surrogate_fun_glm_AL <- function(beta_new, beta_old, step_size, ...) {
  as.numeric(fun_glm_AL(beta_old, ...) + as.numeric(grad_glm_AL(beta_old, ...)) %*% as.numeric(beta_new - beta_old) + 1 / (2 * step_size) * sum((beta_new - beta_old)^2))
}

proximal_gradient <- function(beta, X, y, wts, gamma, dual_step, family, prox_fun, ..., 
                              tol_type = c('beta', 'fval'), step_size = 1, accelerate = T, line_search = T, 
                              ls_beta = .75, ls_eps = 1e-10, ls_max_iter = 50, maxit = 1000, prox_eps = 1e-6, verbose = F) {
  ##checking inputs
  tol_type <- tol_type[1]
  ##initializing
  beta_curr <- beta
  ##
  ## for acceleration
  if (accelerate) {
    beta_m1 <- beta_curr
    beta_m2 <- beta_curr
  }
  ##
  fun_curr <- fun_glm_AL(beta, X, y, wts, gamma, dual_step, family)
  fun_vals <- rep(0, maxit)
  tiny <- .Machine$double.eps
  tiny2 <- .Machine$double.eps^(.25)
  j <- 0
  converged <- F
  while(!converged && j < maxit) {
    j <- j + 1
    ####acceleration with previous two iterates
    if (accelerate) {
      nu <- beta_m1 + (j - 2) / (j + 1) * (beta_m1 - beta_m2)
      grad_val <- grad_glm_AL(nu, X, y, wts, gamma, dual_step, family)
      beta_new <- prox_fun(nu - step_size * grad_val, ...)
      if (line_search) {
        fun_diff <- as.numeric(fun_glm_AL(beta_new, X, y, wts, gamma, dual_step, family) - surrogate_fun_glm_AL(beta_new, nu, step_size, X = X, y = y, wts = wts, gamma = gamma, dual_step = dual_step, family = family))
        ls_iter <- 0
        while(fun_diff > 0 && abs(fun_diff) > ls_eps && ls_iter < ls_max_iter) {
          ls_iter <- ls_iter + 1
          step_size <- step_size * ls_beta
          beta_new <- prox_fun(nu - step_size * grad_val, ...)
          fun_diff <- as.numeric(fun_glm_AL(beta_new, X, y, wts, gamma, dual_step, family) - surrogate_fun_glm_AL(beta_new, nu, step_size, X = X, y = y, wts = wts, gamma = gamma, dual_step = dual_step, family = family))
        }
        if (ls_iter == ls_max_iter) {
          warning(paste('Line search unsuccessful at iteration', j))
        }
      }
    } else {
      step_size <- 1
      grad_val <- grad_glm_AL(beta_curr, X, y, wts, gamma, dual_step, family)
      beta_new <- prox_fun(beta_curr - step_size * grad_val, ...)
      if (line_search) {
        fun_diff <- as.numeric(fun_glm_AL(beta_new, X, y, wts, gamma, dual_step, family) - surrogate_fun_glm_AL(beta_new, beta_curr, step_size, X, y, wts, gamma, dual_step, family))
        ls_iter <- 0
        while(fun_diff > 0 && abs(fun_diff) > ls_eps && ls_iter < ls_max_iter) {
          ls_iter <- ls_iter + 1
          step_size <- step_size * ls_beta
          beta_new <- prox_fun(beta_curr - step_size * grad_val, ...)
          fun_diff <- as.numeric(fun_glm_AL(beta_new, X, y, wts, gamma, dual_step, family) - surrogate_fun_glm_AL(beta_new, beta_curr, step_size, X, y, wts, gamma, dual_step, family))
          # print(paste(fun_diff, j, step_size))
        }
        if (ls_iter == ls_max_iter) {
          warning(paste('Line search unsuccessful at iteration', j))
        }
      }
    }
    fun_new <- fun_glm_AL(beta_new, X, y, wts, gamma, dual_step, family)
    fun_sc <- ifelse(abs(fun_curr) < tiny, tiny2, abs(fun_curr))
    fun_tol <- max(abs(fun_new - fun_curr)) / fun_sc
    fun_curr <- fun_new
    fun_vals[j] <- fun_curr
    beta_sc <- ifelse(max(abs(beta_curr)) < tiny, tiny2, max(abs(beta_curr)))
    beta_tol <- max(abs(beta_new - beta_curr)) / beta_sc
    ## update values
    beta_curr <- beta_new
    if (accelerate) {
      beta_m2 <- beta_m1
      beta_m1 <- beta_new
    }
    ##
    tol <- switch(tol_type,
                  beta = beta_tol,
                  fval = fun_tol)
    if (tol < prox_eps && !any(is.na(tol) || is.nan(tol) || is.infinite(tol))) {
      converged <- T
    }
  }
  fun_vals <- fun_vals[1:j]
  return(list(beta_new = beta_curr, prox_step = step_size, fval = fun_vals, iter = j, converged = converged))
}

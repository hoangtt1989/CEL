quantile_thresh <- function(inp, lambda) {
  outp <- rep(0, length(inp))
  inp_dec <- order(abs(inp), decreasing = T)
  keep_idx <- inp_dec[1:lambda]
  outp[keep_idx] <- inp[keep_idx]
  return(outp)
}

fun_outlier_lm_AL <- function(delta, beta, wts, X, y, gamma, dual_step) {
  XtrW <- t(wts * X)
  loss_term <- as.vector(XtrW %*% (y - X %*% beta - delta))
  as.vector(gamma) %*% loss_term + .5 * dual_step * sum((loss_term)^2) 
}

grad_outlier_lm_AL <- function(delta, beta, wts, X, y, gamma, dual_step) {
  XtrW <- t(wts * X)
  WX <- t(XtrW)
  as.vector(dual_step * WX %*% (XtrW %*% (delta - y + X %*% beta) - gamma / dual_step))
}

surrogate_outlier_lm_AL <- function(delta_new, delta_curr, beta, wts, X, y, gamma, dual_step, step_size) {
  as.numeric(fun_outlier_lm_AL(delta_curr, beta, wts, X, y, gamma, dual_step)) +
    as.numeric(grad_outlier_lm_AL(delta_curr, beta, wts, X, y, gamma, dual_step)) %*% as.numeric(delta_new - delta_curr) +
    1 / (2 * step_size) * as.numeric(sum((delta_new - delta_curr)^2))
}

outlier_loop <- function(wts_new, X, F_reparam, y, gamma_curr, dual_step, s_hat_reparam, beta_opt_type, 
                         lbfgs_arg, fun_reparam_lm_AL, grad_reparam_lm_AL, beta_new = NULL, delta_new = NULL, q = NULL,
                         maxit = 1000, outlier_eps = 1e-5, line_search = T, ls_eps = 1e-12, ls_max_iter = 100, ls_beta = .75) {
  
  n <- length(wts_new)
  if (is.null(delta_new)) {
    delta_init <- rep(0, n)
  } else {
    delta_init <- delta_new
  }
  
  y_adj <- y - delta_init
  delta_curr <- delta_init
  
  WX <- wts_new * X
  XtrW <- t(WX)
  if (line_search) {
    delta_step <- 1
  } else {
    delta_step <- 1 / max(svd(XtrW)$d)^2
  }
  if (is.null(q)) {
    q <- floor(.1 * length(wts_new))
  }
  
  converged <- F
  
  beta_curr <- rep(0, n)
  
  for (i in 1:maxit) {
    if (line_search) {
      delta_step <- 1
    }
    ####beta update
    # beta_res <- beta_wrapper(beta_new, wts_new, X, F_reparam, y_adj, gamma_curr, dual_step, s_hat_reparam, beta_opt_type, lbfgs_arg, fun_reparam_lm_AL, grad_reparam_lm_AL)
    # beta_new <- beta_res$beta_new
    ####
    ####delta (outlier) update
    delta_grad_val <- grad_outlier_lm_AL(delta_curr, beta_new, wts_new, X, y, gamma_curr, dual_step)
    delta_new <- quantile_thresh(delta_curr - delta_step * delta_grad_val, q)
    # delta_new <- quantile_thresh(delta_curr - dual_step * delta_step * WX %*% (XtrW %*% (delta_curr - y + X %*% beta_res$beta_new) - gamma_curr / dual_step), q)
    if (line_search) {
      fun_diff <- as.numeric(fun_outlier_lm_AL(delta_new, beta_new, wts_new, X, y, gamma_curr, dual_step) - surrogate_outlier_lm_AL(delta_new, delta_curr, beta_new, wts_new, X, y, gamma_curr, dual_step, delta_step))
      ls_iter <- 0
      while(fun_diff > 0 && abs(fun_diff) > ls_eps && ls_iter < ls_max_iter) {
        ls_iter <- ls_iter + 1
        delta_step <- delta_step * ls_beta
        delta_new <- quantile_thresh(delta_curr - delta_step * delta_grad_val, q)
        fun_diff <- as.numeric(fun_outlier_lm_AL(delta_new, beta_new, wts_new, X, y, gamma_curr, dual_step) - surrogate_outlier_lm_AL(delta_new, delta_curr, beta_new, wts_new, X, y, gamma_curr, dual_step, delta_step))
      }
    }
    ####
    delta_diff <- max(abs(delta_new - delta_curr))
    # delta_diff <- sum((delta_new == 0) - (delta_curr == 0)) == 0
    # beta_diff <- as.numeric(max(abs(beta_curr - beta_new)))
    delta_curr <- delta_new
    beta_curr <- beta_new
    y_adj <- y - delta_curr
    #### or delta_diff < outlier_eps
    if (delta_diff < outlier_eps) {
      converged <- T
      break
    }
  }
  R_adj_new <- as.vector(y - X %*% beta_new - delta_curr)
  # R_adj_new <- beta_res$R_new
  return(list(beta_new = beta_new, delta_new = delta_curr, R_adj_new = R_adj_new, y_adj = y_adj, converged = converged, iter = i, tol = delta_diff))
}

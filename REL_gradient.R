quantile_thresh <- function(inp, lambda) {
  n <- length(inp)
  lambda <- min(length(inp), floor(lambda)) ##make sure lambda is an integer, less than or equal to n
  outp <- rep(0, length(inp))
  if (lambda > 0) { ##if lambda <= 0, don't do any thresholding
    inp_dec <- order(abs(inp), decreasing = T)
    keep_idx <- inp_dec[1:lambda]
    outp[keep_idx] <- inp[keep_idx]
  }
  return(outp)
}


fun_delta_lm_AL <- function(wts, delta, gamma, X, R, dual_step, delta_plus = F) {
  if (delta_plus) {
    delta <- -delta
  }
  inner_term <- as.vector(crossprod(R * X, wts - delta))
  as.vector(gamma) %*% inner_term + .5 * dual_step * sum((inner_term)^2)
}

grad_delta_lm_AL <- function(wts, delta, gamma, X, R, dual_step, delta_plus = F) {
  if (delta_plus) {
    delta <- -delta
  }
  RX <- R * X
  ret <- -RX %*% gamma - dual_step * tcrossprod(RX) %*% (wts - delta)
  if (delta_plus) {
    ret <- -ret
  }
  return(ret)
}

surrogate_outlier_lm_AL <- function(delta_new, delta_curr, wts, gamma, X, R, dual_step, step_size, delta_plus = F) {
  as.numeric(fun_delta_lm_AL(wts, delta_curr, gamma, X, R, dual_step, delta_plus)) +
    as.numeric(grad_delta_lm_AL(wts, delta_curr, gamma, X, R, dual_step, delta_plus)) %*% as.numeric(delta_new - delta_curr) +
    1 / (2 * step_size) * as.numeric(sum((delta_new - delta_curr)^2))
}


outlier_loop <- function(wts_new, X, R, gamma_curr, dual_step, 
                         delta_new = NULL, q = NULL, delta_pos = F, delta_plus = F, maxit = 1000, accelerate = T,
                         outlier_eps = 1e-5, line_search = T, ls_eps = 1e-12, 
                         ls_max_iter = 100, ls_beta = .75) {
  print('REL outlier')
  n <- length(wts_new)
  if (is.null(delta_new)) {
    delta_init <- rep(0, n)
    # delta_init <- rnorm(n)
  } else {
    delta_init <- delta_new
  }
  
  delta_curr <- delta_init
  
  RX <- R * X

  if (line_search) {
    delta_step <- 1
  } else {
    delta_step <- 1 / max(svd(dual_step * t(RX))$d)^2
  }
  if (is.null(q)) {
    q <- floor(.1 * length(wts_new))
  }
  
  ## for acceleration
  if (accelerate) {
    delta_m1 <- delta_curr
    delta_m2 <- delta_curr
  }
  ##
  
  converged <- F
  
  if (line_search) {
    delta_step <- 1
  }
  
  for (i in 1:maxit) {
    ####acceleration with previous two iterates
    if (accelerate) {
      nu <- delta_m1 + (i - 2) / (i + 1) * (delta_m1 - delta_m2)
      delta_grad_val <- grad_delta_lm_AL(wts_new, nu, gamma_curr, X, R, dual_step, delta_plus)
      update_term <- nu - delta_step * delta_grad_val
      if (delta_pos) {
        update_term <- pmax(0, update_term)
      }
      delta_new <- quantile_thresh(update_term, q)
      if (line_search) {
        fun_diff <- as.numeric(fun_delta_lm_AL(wts_new, delta_new, gamma_curr, X, R, dual_step, delta_plus) -
                                 surrogate_outlier_lm_AL(delta_new, nu, wts_new, gamma_curr, X, R, dual_step, delta_step, delta_plus))
        ls_iter <- 0
        while(fun_diff > 0 && abs(fun_diff) > ls_eps && ls_iter < ls_max_iter) {
          ls_iter <- ls_iter + 1
          delta_step <- delta_step * ls_beta
          update_term <- nu - delta_step * delta_grad_val
          if (delta_pos) {
            update_term <- pmax(0, update_term)
          }
          delta_new <- quantile_thresh(update_term, q)
          fun_diff <- as.numeric(fun_delta_lm_AL(wts_new, delta_new, gamma_curr, X, R, dual_step, delta_plus) -
                                   surrogate_outlier_lm_AL(delta_new, nu, wts_new, gamma_curr, X, R, dual_step, delta_step, delta_plus))
        }
        if (ls_iter == ls_max_iter) {
          warning(paste('Line search unsuccessful at iteration', j))
        }
      }
    } else {
      if (line_search) {
        delta_step <- 1
      }
      ####delta (outlier) update
      delta_grad_val <- grad_delta_lm_AL(wts_new, delta_curr, gamma_curr, X, R, dual_step, delta_plus)
      update_term <- delta_curr - delta_step * delta_grad_val
      if (delta_pos) {
        update_term <- pmax(0, update_term)
      }
      delta_new <- quantile_thresh(update_term, q)
      if (line_search) {
        fun_diff <- as.numeric(fun_delta_lm_AL(wts_new, delta_new, gamma_curr, X, R, dual_step, delta_plus) -
                                 surrogate_outlier_lm_AL(delta_new, delta_curr, wts_new, gamma_curr, X, R, dual_step, delta_step, delta_plus))
        ls_iter <- 0
        while(fun_diff > 0 && abs(fun_diff) > ls_eps && ls_iter < ls_max_iter) {
          ls_iter <- ls_iter + 1
          delta_step <- delta_step * ls_beta
          update_term <- delta_curr - delta_step * delta_grad_val
          if (delta_pos) {
            update_term <- pmax(0, update_term)
          }
          delta_new <- quantile_thresh(update_term, q)
          fun_diff <- as.numeric(fun_delta_lm_AL(wts_new, delta_new, gamma_curr, X, R, dual_step, delta_plus) -
                                   surrogate_outlier_lm_AL(delta_new, delta_curr, wts_new, gamma_curr, X, R, dual_step, delta_step, delta_plus))
        }
      }
      ####
    }
    delta_diff <- max(abs(delta_new - delta_curr))
    # delta_diff <- sum((delta_new == 0) - (delta_curr == 0)) == 0
    # beta_diff <- as.numeric(max(abs(beta_curr - beta_new)))
    delta_curr <- delta_new
    if (accelerate) {
      delta_m1 <- delta_curr
      delta_m2 <- delta_curr
    }
    if (delta_diff < outlier_eps) {
      converged <- T
      break
    }
  }
  return(list(delta_new = delta_curr, converged = converged, iter = i, tol = delta_diff))
}

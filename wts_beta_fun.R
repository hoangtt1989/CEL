##wrapper functions for wts/beta updates in ADMM, so their order can be swapped
# ###wts wrapper
# wts_wrapper <- function(n, fun, grad, wts_new, fun_arg, grad_arg, mirror_arg) {
#   A_max <- fun_arg$dual_step * max(abs(crossprod(fun_arg$Xtr_R))) * n ##t(Xtr_R) %*% Xtr_R
#   mirror_step <- 1 / (A_max + max(abs(fun_arg$gammatr_Xtr_R)))
#   mirror_res <- do.call('mirror_descent', c(list(fun = fun, grad = grad, wts = wts_new, fun_arg = fun_arg, grad_arg = grad_arg, mirror_step = mirror_step), mirror_arg), quote = T)
#   wts_new <- as.vector(mirror_res$wts)
#   return(list(wts_new = wts_new, nits = mirror_res$iter, converged = mirror_res$converged))
# }
# ###
###beta wrapper
beta_wrapper <- function(beta_new, wts_new, X, F_reparam, y, gamma_curr, dual_step, s_hat_reparam, beta_opt_type, lbfgs_arg, fun_reparam_lm_AL, grad_reparam_lm_AL) {
  # W <- Diagonal(x = as.numeric(wts_curr))
  # Xtr_W_X <- crossprod(X, W %*% X) ##t(X) %*% W %*% X
  if (beta_opt_type == 'closed_form') {
    z_hat_new <- closed_beta_lm_AL_cpp(wts_new, gamma_curr, y, X, F_reparam, s_hat_reparam, dual_step) ## use a cpp function for the closed form solution
    # z_hat_new <- solve(crossprod(F_reparam, crossprod(Xtr_W_X) %*% F_reparam), crossprod(F_reparam, Xtr_W_X) %*% (crossprod(X, wts_new * y) + gamma_curr / dual_step - Xtr_W_X %*% s_hat_reparam))
    converged <- T
    # z_hat_new <- solve(crossprod(F_reparam, crossprod(Xtr_W_X) %*% F_reparam), crossprod(F_reparam, Xtr_W_X) %*% (crossprod(X, W %*% y) + gamma_curr / dual_step - Xtr_W_X %*% s_hat_reparam))
  } else if (beta_opt_type == 'LBFGS') {
    Xtr_W_X <- crossprod(X, as.vector(wts_new) * X)
    Xtr_W_X_F <- Xtr_W_X %*% F_reparam
    Ftr_Xtr_W_X_gamma <- as.vector(crossprod(Xtr_W_X_F, gamma_curr))
    # gammatr_Xtr_W_X_F <- crossprod(gamma_curr, Xtr_W_X_F)
    Ftr_Xtr_W_X_Xtr_W_y <- as.vector(crossprod(Xtr_W_X_F, crossprod(X, wts_new * y)))
    Ftr_Xtr_W_X_Xtr_W_X_F <- crossprod(Xtr_W_X_F)
    Ftr_Xtr_W_X_Xtr_W_X_shat <- as.vector(crossprod(Xtr_W_X_F, Xtr_W_X %*% s_hat_reparam))
    z_hat_curr <- as.vector(crossprod(F_reparam, beta_new - s_hat_reparam))
    opt <- do.call('lbfgs', c(list(fun_reparam_lm_AL, grad_reparam_lm_AL, z_hat_curr, 
                                   Ftr_Xtr_W_X_gamma = Ftr_Xtr_W_X_gamma, Ftr_Xtr_W_X_Xtr_W_y = Ftr_Xtr_W_X_Xtr_W_y, 
                                   Ftr_Xtr_W_X_Xtr_W_X_F = Ftr_Xtr_W_X_Xtr_W_X_F, 
                                   Ftr_Xtr_W_X_Xtr_W_X_shat = Ftr_Xtr_W_X_Xtr_W_X_shat, 
                                   dual_step = dual_step), lbfgs_arg), quote = T)
    z_hat_new <- as.vector(opt$par)
    converged <- opt$convergence == 0
    # if (verbose) {
    #   if (opt$convergence < 0) {
    #     warning(paste('Errors in execution of lbfgs at iteration', j))
    #   }
    # }
  }
  beta_new <- as.vector(F_reparam %*% z_hat_new + s_hat_reparam)
  # Xtr_R <- crossprod(X, R_curr)
  R_new <- as.vector(y - X %*% beta_new)
  # Xtr_R <- t(R_new * X)
  # ##update fun and grad arguments
  # fun_arg$Xtr_R <- Xtr_R
  # grad_arg$Xtr_R <- Xtr_R
  # ##
  return(list(beta_new = beta_new, R_new = R_new, converged = converged))
}

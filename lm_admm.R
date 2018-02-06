# library(Matrix) ##for sparse Diagonal
# library(pracma) ##for orth, nullspace
library(compiler) ##for cmpfun
library(lbfgs) ##for lbfgs
library(Rcpp)
library(RcppEigen)
source('mirror_descent.R') ##for mirror descent, function and gradient calculations
source('wts_beta_fun.R')
Rcpp::sourceCpp('rcpp_fun.cpp')

###function value for reparameterized beta
fun_reparam_lm_AL <- function(z_hat, Ftr_Xtr_W_X_gamma, Ftr_Xtr_W_X_Xtr_W_y, Ftr_Xtr_W_X_Xtr_W_X_F, Ftr_Xtr_W_X_Xtr_W_X_shat, dual_step) {
  res <- -Ftr_Xtr_W_X_gamma %*% z_hat + dual_step * (-z_hat %*% Ftr_Xtr_W_X_Xtr_W_y + .5 * z_hat %*% (Ftr_Xtr_W_X_Xtr_W_X_F %*% z_hat) + z_hat %*% Ftr_Xtr_W_X_Xtr_W_X_shat)
  as.numeric(res)
}
###
###gradient for reparameterized beta
grad_reparam_lm_AL <- function(z_hat, Ftr_Xtr_W_X_gamma, Ftr_Xtr_W_X_Xtr_W_y, Ftr_Xtr_W_X_Xtr_W_X_F, Ftr_Xtr_W_X_Xtr_W_X_shat, dual_step) {
  res <- -Ftr_Xtr_W_X_gamma + dual_step * (-Ftr_Xtr_W_X_Xtr_W_y + Ftr_Xtr_W_X_Xtr_W_X_F %*% z_hat + Ftr_Xtr_W_X_Xtr_W_X_shat)
  as.numeric(res)
}
###


##################function for CEL linear regression admm (Rcpp)
CEL_lm_admm <- function(beta_test, X, y, F_reparam, s_hat_reparam, wts_init = NULL, beta_opt_type = c('closed_form', 'LBFGS'), vary_penalty = c('RB', 'RB2', 'none'), 
                        RB_mu = 10, RB_tau = 2, RB2_ksi = 2, outer_eps = 1e-8, outer_rel_eps = 1e-4, dual_step = 2, outer_maxit = 1000, dual_type = 1, wts_beta_rep = 1, random_order = F,
                        outer_tol_type = c('fval', 'primal_dual', 'gamma', 'beta', 'beta_gamma', 'primal'), 
                        mirror_arg = list(), lbfgs_arg = list(gtol = .1, invisible = 1), verbose = F) {
  ##checking inputs
  beta_opt_type <- beta_opt_type[1]
  outer_tol_type <- outer_tol_type[1]
  vary_penalty <- vary_penalty[1]
  primal_dual_flag <- ifelse(vary_penalty %in% c('RB', 'RB2') || outer_tol_type %in% c('primal', 'primal_dual'), T, F)
  beta_diff_flag <- ifelse(outer_tol_type == 'beta' || outer_tol_type == 'beta_gamma', T, F)
  gamma_diff_flag <- ifelse(outer_tol_type == 'gamma' || outer_tol_type == 'beta_gamma', T, F)
  if (!(outer_tol_type %in% c('gamma', 'beta', 'beta_gamma', 'fval', 'primal', 'primal_dual'))) {
    stop('supply a valid outer tolerance type')
  }
  if (!(vary_penalty %in% c('RB', 'RB2', 'none'))) {
    stop('supply a valid penalty varying type')
  }
  ##
  ##initializing function and gradient
  # fun <- fun_lm_AL_cpp
  # grad <- grad_lm_AL_cpp
  # mirror_descent <- cmpfun(mirror_descent)
  # wts_wrapper <- cmpfun(wts_wrapper)
  beta_wrapper <- cmpfun(beta_wrapper)
  fun_reparam_lm_AL <- cmpfun(fun_reparam_lm_AL)
  grad_reparam_lm_AL <- cmpfun(grad_reparam_lm_AL)
  ##
  ##initial parameters
  prob_size <- dim(X)
  n <- prob_size[1]
  p <- prob_size[2]
  beta_curr <- beta_test
  # z_hat_curr <- as.vector(crossprod(F_reparam, beta_test - s_hat_reparam))
  if (is.null(wts_init)) {
    wts_init <- rep(1, n)
  } else {
    if (length(wts_init) != n) {
      stop('wts_init must be a vector with n elements')
    }
  }
  wts_new <- wts_init
  wts_curr <- wts_new
  gamma_curr <- rep(0, p)
  outer_converged <- F
  ##
  ##primal_dual stopping
  if (outer_tol_type == 'primal_dual') {
    primal_eps <- outer_eps * sqrt(p)
  }
  # dual_eps <- outer_eps * sqrt(n) ## dimension of dual residual changes depending on update order
  ##
  ##initialize holders
  beta_converged_vec <- rep(NA, outer_maxit)
  outer_fval <- rep(NA, outer_maxit)
  mirror_nits <- rep(NA, outer_maxit)
  update_order <- rep(NA, outer_maxit)
  rho_vec <- rep(NA, outer_maxit)
  ##
  ####precomputing
  R_new <- as.vector(y - X %*% beta_curr)
  R_curr <- R_new
  Xtr_R <- t(R_new * X) #t(X) %*% R_curr
  ###initialize function/gradient arguments
  fun_grad_arg <- list(gamma = gamma_curr, X = X, R = R_curr, dual_step = dual_step)
  ##
  # outer_fval_curr <- do.call(fun, c(list(wts = wts_curr), fun_arg))
  outer_fval_curr <- sum(log(wts_new))
  ####
  #################main loop
  start_time <- Sys.time()
  for (j in 1:outer_maxit) {
    if (verbose) {
      if (! j %% 50) {
        message(paste('Iteration', j))
      }
    }
    ##random ordering of wts/beta update
    if (!random_order || runif(1) < .5) {
      wts_beta_order <- T ##update wts then beta
      dual_resid_elt <- sqrt(n) ##dual residual has n elements
      update_order[j] <- 'wts_beta'
    } else {
      wts_beta_order <- F ##update beta then wts
      dual_resid_elt <- sqrt(p) ##dual residual has p elements
      update_order[j] <- 'beta_wts'
    }
    if (outer_tol_type == 'primal_dual') {
      dual_eps <- outer_eps * dual_resid_elt
    }
    ##
    ##wts and beta updates can be repeated
    for (wts_beta_iter in 1:wts_beta_rep) {
      if (wts_beta_order) {
        #######mirror descent
        A_max <- dual_step * max(abs(crossprod(Xtr_R))) * n ##t(Xtr_R) %*% Xtr_R
        mirror_step <- 1 / (A_max + max(abs(crossprod(gamma_curr, Xtr_R))))
        mirror_res <- do.call('mirror_descent_lm_AL_cpp', c(list(wts = wts_new, mirror_step = mirror_step), fun_grad_arg, mirror_arg), quote = T)
        # mirror_res <- do.call('mirror_descent', c(list(fun = fun, grad = grad, wts = wts_new, fun_arg = fun_grad_arg, grad_arg = fun_grad_arg, mirror_step = mirror_step), mirror_arg), quote = T)
        wts_new <- as.vector(mirror_res$wts)
        # wts_res <- wts_wrapper(n, fun, grad, wts_new, fun_arg, grad_arg, mirror_arg)
        # wts_new <- wts_res$wts_new
        # mirror_iter <- wts_res$nits
        # mirror_converged <- wts_res$converged
        #######
        #######beta update
        beta_res <- beta_wrapper(beta_curr, wts_new, X, F_reparam, y, gamma_curr, dual_step, s_hat_reparam, beta_opt_type, lbfgs_arg, fun_reparam_lm_AL, grad_reparam_lm_AL)
        beta_new <- beta_res$beta_new
        R_new <- beta_res$R_new
        Xtr_R <- t(R_new * X)
        # Xtr_R <- beta_res$Xtr_R
        beta_converged <- beta_res$converged
        ##update fun and grad arguments
        fun_grad_arg$R <- R_new
        ##
        #######
      } else {
        #######beta update
        beta_res <- beta_wrapper(beta_curr, wts_new, X, F_reparam, y, gamma_curr, dual_step, s_hat_reparam, beta_opt_type, lbfgs_arg, fun_reparam_lm_AL, grad_reparam_lm_AL)
        beta_new <- beta_res$beta_new
        R_new <- beta_res$R_new
        Xtr_R <- t(R_new * X)
        # Xtr_R <- beta_res$Xtr_R
        beta_converged <- beta_res$converged
        ##update fun and grad arguments
        fun_grad_arg$R <- R_new
        # fun_arg$Xtr_R <- Xtr_R
        # grad_arg$Xtr_R <- Xtr_R
        ##
        #######
        #######mirror descent
        A_max <- dual_step * max(abs(crossprod(Xtr_R))) * n ##t(Xtr_R) %*% Xtr_R
        mirror_step <- 1 / (A_max + max(abs(crossprod(gamma_curr, Xtr_R))))
        mirror_res <- do.call('mirror_descent_lm_AL_cpp', c(list(wts = wts_new, mirror_step = mirror_step), fun_grad_arg, mirror_arg), quote = T)
        # mirror_res <- do.call('mirror_descent', c(list(fun = fun, grad = grad, wts = wts_new, fun_arg = fun_grad_arg, grad_arg = fun_grad_arg, mirror_step = mirror_step), mirror_arg), quote = T)
        wts_new <- as.vector(mirror_res$wts)
        # wts_res <- wts_wrapper(n, fun, grad, wts_new, fun_arg, grad_arg, mirror_arg)
        # wts_new <- wts_res$wts_new
        # mirror_iter <- wts_res$nits
        # mirror_converged <- wts_res$converged
        #######
      }
    }
    ##
    mirror_nits[j] <- mirror_res$iter
    beta_converged_vec[j] <- beta_converged
    if (verbose) {
      if (!mirror_res$converged) {
        warning(paste('Mirror descent did not converge at iteration', j))
      }
      if (!beta_res$converged) {
        warning(paste('Errors in execution of lbfgs at iteration', j))
      }
    }
    ##update beta_diff
    if (beta_diff_flag) {
      beta_sc <- max(abs(beta_curr))
      beta_sc <- ifelse(beta_sc < .Machine$double.eps, .Machine$double.eps * 2, beta_sc)
      beta_diff <- max(abs(beta_curr - beta_new)) / beta_sc
    }
    ##
    ##update R_diff and wts_diff for dual residual
    if (wts_beta_order) { ##wts were updated first
      R_diff <- as.vector(R_new - R_curr)
    } else { ##beta was updated first
      wts_diff <- as.vector(wts_new - wts_curr)
    }
    ##
    #######gamma update
    gamma_new <- gamma_curr + dual_step * Xtr_R %*% wts_new
    if (gamma_diff_flag) {
      gamma_sc <- max(abs(gamma_curr))
      gamma_sc <- ifelse(gamma_sc < .Machine$double.eps, .Machine$double.eps * 2, gamma_sc)
      gamma_diff <- max(abs(gamma_new - gamma_curr)) / gamma_sc
    }
    #######
    #######dual residual
    if (primal_dual_flag) {
      if (wts_beta_order) { ##wts were updated first
        # dual_resid <- sum((dual_step * R_curr %*% tcrossprod(X) %*% ((R_new - R_curr) %*% wts_curr))^2)
        dual_resid <- sqrt(sum((dual_step * tcrossprod(R_curr * X, R_diff * X) %*% wts_new)^2))
        dual_resid_scale <- sqrt(sum((R_curr * X %*% gamma_new)^2)) 
      } else { ##beta was updated first
        # dual_resid <- sqrt(sum(dual_step * tcrossprod(X, wts_curr * X) %*% (tcrossprod(X, wts_diff * y) - tcrossprod(X, wts_diff * X * beta_new)))^2)
        dual_resid <- sqrt(sum((dual_step * crossprod(X, wts_curr * tcrossprod(X)) %*% (wts_diff * y - wts_diff * (X %*% beta_new)))^2))
      }
    }
    #######
    #######primal residual
    if (primal_dual_flag) {
      primal_resid <- sqrt(sum((Xtr_R %*% wts_new)^2))
      wts_X <- wts_new * X
      primal_resid_scale <- max(sqrt(sum((crossprod(wts_X, y))^2)), sqrt(sum((crossprod(wts_X, X %*% beta_new))^2)))
    }
    #######
    ##update R, gamma, beta and wts
    beta_curr <- as.vector(beta_new)
    R_curr <- as.vector(R_new)
    wts_curr <- as.vector(wts_new)
    gamma_curr <- gamma_new
    ##
    ####update fun_arg and grad_arg
    fun_grad_arg$gamma <- gamma_curr
    ####
    ###function value update
    # outer_fval_new <- do.call(fun, c(list(wts = wts_curr), fun_arg))
    outer_fval_new <- sum(log(wts_curr))
    outer_fval[j] <- outer_fval_new
    if (outer_tol_type == 'fval') {
      outer_fval_sc <- max(abs(outer_fval_curr))
      outer_fval_sc <- ifelse(outer_fval_sc < .Machine$double.eps, .Machine$double.eps * 2, outer_fval_sc)
      outer_fval_diff <- abs(outer_fval_new - outer_fval_curr) / abs(outer_fval_sc)
    }
    outer_fval_curr <- outer_fval_new
    ###
    #######vary penalty parameter
    if (vary_penalty == 'RB') {
      if (primal_resid > RB_mu * dual_resid) {
        dual_step <- RB_tau * dual_step
      } else if (dual_resid > RB_mu * primal_resid) {
        dual_step <- dual_step / RB_tau
      }
      ##update fun_arg and grad_arg
      fun_grad_arg$dual_step <- dual_step
      ##
    }
    if (vary_penalty == 'RB2') {
      ###varying tau
      primal_dual_div <- primal_resid / dual_resid
      dual_primal_div <- 1 / primal_dual_div
      sqrt_primal_dual_div <- sqrt(primal_dual_div / RB2_ksi)
      sqrt_dual_primal_div <- sqrt(dual_primal_div * RB2_ksi)
      if (sqrt_primal_dual_div < RB_tau && sqrt_primal_dual_div >= 1) {
        RB2_tau <- sqrt_primal_dual_div
      } else if (sqrt_dual_primal_div < 1 && sqrt_dual_primal_div > (1 / RB_tau)) {
        RB2_tau <- sqrt_primal_dual_div
      } else {
        RB2_tau <- RB_tau
      }
      ###
      primal_resid_rel <- primal_resid / primal_resid_scale
      dual_resid_rel <- dual_resid / dual_resid_scale
      if (primal_resid_rel > RB2_ksi * RB_mu * dual_resid_rel) {
        dual_step <- RB2_tau * dual_step
      } else if (dual_resid_rel > RB_mu * primal_resid_rel / RB2_ksi) {
        dual_step <- dual_step / RB2_tau
      }
      ##update fun_arg and grad_arg
      fun_grad_arg$dual_step <- dual_step
      ##
    }
    rho_vec[j] <- dual_step
    #######
    #######stopping
    if (outer_tol_type == 'gamma') {
      outer_tol <- gamma_diff
    } else if (outer_tol_type == 'beta') {
      outer_tol <- beta_diff
    } else if (outer_tol_type == 'beta_gamma') {
      outer_tol <- min(beta_diff, gamma_diff)
    } else if (outer_tol_type == 'fval') {
      outer_tol <- outer_fval_diff
    } else if (outer_tol_type %in% c('primal', 'primal_dual')) {
      outer_tol <- primal_resid
    }
    if (outer_tol_type == 'primal_dual') {
      primal_eps_scale <- primal_eps + outer_rel_eps * primal_resid_scale / n
      dual_eps_scale <- dual_eps + outer_rel_eps * dual_resid_scale
      if (primal_resid / n < primal_eps_scale && dual_resid < dual_eps_scale) {
        outer_converged <- T
        break
      }
    } else {
      if (all(outer_tol < outer_eps) & !any(is.nan(outer_tol) || is.na(outer_tol) || is.infinite(outer_tol))) {
        outer_converged <- T
        break
      }
    }
    #######
  }
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  #################
  if (!outer_converged) {
    warning(paste('Outer loop did not converge, current tolerance is', outer_tol))
  }
  
  ####cleaning up output
  rho_vec <- rho_vec[1:j]
  outer_fval <- outer_fval[1:j]
  mirror_nits <- mirror_nits[1:j]
  beta_converged_vec <- beta_converged_vec[1:j]
  update_order <- update_order[1:j]
  # logelr <- sum(log(wts_curr))
  logelr <- outer_fval_curr
  if (!primal_dual_flag) {
    primal_resid <- 'not calculated'
    dual_resid <- 'not calculated'
  }
  ####
  
  return(list(logelr = logelr, logelr_stat = -2 * logelr, wts = wts_curr / n, sum_wts = sum(wts_curr) / n,
              gamma = gamma_curr, beta_opt = beta_curr, outer_converged = outer_converged, time = as.double(tot_time, 'secs'), 
              outer_fval = -outer_fval, primal_resid = primal_resid, dual_resid = dual_resid, rho = rho_vec,
              outer_tol = outer_tol, update_order = update_order, beta_converged = beta_converged_vec, mirror_nits = mirror_nits, outer_nits = j))
}
##################


# ##################function for CEL linear regression admm
# CEL_lm_admm <- function(beta_test, X, y, F_reparam, s_hat_reparam, beta_opt_type = c('closed_form', 'LBFGS'), vary_penalty = c('RB', 'none'), 
#                            RB_mu = 10, RB_tau = 2, outer_eps = 1e-8, dual_step = 2, outer_maxit = 1000, dual_type = 1, wts_beta_rep = 1,
#                            outer_tol_type = c('fval', 'primal_dual', 'gamma', 'beta', 'beta_gamma', 'primal'), mirror_arg = list(), lbfgs_arg = list(epsilon = 1e-8, invisible = 1), verbose = F) {
#   ##checking inputs
#   beta_opt_type <- beta_opt_type[1]
#   outer_tol_type <- outer_tol_type[1]
#   vary_penalty <- vary_penalty[1]
#   primal_dual_flag <- ifelse(vary_penalty == 'RB' || outer_tol_type %in% c('primal', 'primal_dual'), T, F)
#   beta_diff_flag <- ifelse(outer_tol_type == 'beta' || outer_tol_type == 'beta_gamma', T, F)
#   gamma_diff_flag <- ifelse(outer_tol_type == 'gamma' || outer_tol_type == 'beta_gamma', T, F)
#   if (!(outer_tol_type %in% c('gamma', 'beta', 'beta_gamma', 'fval', 'primal', 'primal_dual'))) {
#     stop('supply a valid outer tolerance type')
#   }
#   if (!(vary_penalty %in% c('RB', 'none'))) {
#     stop('supply a valid penalty varying type')
#   }
#   ##
#   ##initializing function and gradient
#   fun <- cmpfun(fun_lm_AL) ##compile function
#   grad <- cmpfun(grad_lm_AL) ##compile function
#   mirror_descent_cmp <- cmpfun(mirror_descent)
#   ##
#   ##initial parameters
#   prob_size <- dim(X)
#   n <- prob_size[1]
#   p <- prob_size[2]
#   beta_curr <- beta_test
#   z_hat_curr <- as.vector(crossprod(F_reparam, beta_test - s_hat_reparam))
#   wts_curr <- as.matrix(rep(1, n))
#   gamma_curr <- as.matrix(rep(0, p))
#   outer_converged <- F
#   ##
#   ##primal_dual stopping
#   primal_eps <- outer_eps * sqrt(p)
#   dual_eps <- outer_eps * sqrt(n)
#   ##
#   ##initialize holders
#   outer_fval <- rep(NA, outer_maxit)
#   inner_nits <- rep(NA, outer_maxit)
#   ##
#   ####precomputing
#   # R_curr <- Diagonal(x = as.numeric(y - X %*% beta_curr))
#   R_curr <- as.vector(y - X %*% beta_curr)
#   Xtr_R <- t(R_curr * X) #t(X) %*% R_curr
#   gammatr_Xtr_R <- crossprod(gamma_curr, Xtr_R) #t(gamma_curr) %*% Xtr_R
#   R_X_gamma <- crossprod(Xtr_R, gamma_curr) #t(Xtr_R) %*% gamma_curr
#   ###initialize function/gradient arguments
#   fun_arg <- list(gammatr_Xtr_R = gammatr_Xtr_R, Xtr_R = Xtr_R, dual_step = dual_step)
#   grad_arg <- list(R_X_gamma = R_X_gamma, Xtr_R = Xtr_R, dual_step = dual_step)
#   ##
#   # outer_fval_curr <- do.call(fun, c(list(wts = wts_curr), fun_arg))
#   outer_fval_curr <- sum(log(wts_curr))
#   ####
#   
#   #################main loop
#   start_time <- Sys.time()
#   for (j in 1:outer_maxit) {
#     if (verbose) {
#       if (! j %% 50) {
#         message(paste('Iteration', j))
#       }
#     }
#     ##wts and beta updates can be repeated
#     for (wts_beta_iter in 1:wts_beta_rep) {
#       #######mirror descent
#       A_max <- dual_step * max(abs(crossprod(fun_arg$Xtr_R))) * n ##t(Xtr_R) %*% Xtr_R
#       mirror_step <- 1 / (A_max + max(abs(fun_arg$gammatr_Xtr_R)))
#       mirror_res <- do.call('mirror_descent_cmp', c(list(fun = fun, grad = grad, wts = wts_curr, fun_arg = fun_arg, grad_arg = grad_arg, mirror_step = mirror_step), mirror_arg), quote = T)
#       if (verbose) {
#         if (!mirror_res$converged) {
#           warning(paste('Mirror descent did not converge at iteration', j))
#         }
#       }
#       wts_curr <- as.vector(mirror_res$wts)
#       inner_nits[j] <- mirror_res$iter
#       #######
#       #######beta update
#       # W <- Diagonal(x = as.numeric(wts_curr))
#       # Xtr_W_X <- crossprod(X, W %*% X) ##t(X) %*% W %*% X
#       Xtr_W_X <- crossprod(X, wts_curr * X)
#       if (beta_opt_type == 'closed_form') {
#         z_hat_new <- solve(crossprod(F_reparam, crossprod(Xtr_W_X) %*% F_reparam), crossprod(F_reparam, Xtr_W_X) %*% (crossprod(X, wts_curr * y) + gamma_curr / dual_step - Xtr_W_X %*% s_hat_reparam))
#         # z_hat_new <- solve(crossprod(F_reparam, crossprod(Xtr_W_X) %*% F_reparam), crossprod(F_reparam, Xtr_W_X) %*% (crossprod(X, W %*% y) + gamma_curr / dual_step - Xtr_W_X %*% s_hat_reparam))
#       } else if (beta_opt_type == 'LBFGS') {
#         Xtr_W_X_F <- Xtr_W_X %*% F_reparam
#         Ftr_Xtr_W_X_gamma <- as.vector(crossprod(Xtr_W_X_F, gamma_curr))
#         # gammatr_Xtr_W_X_F <- crossprod(gamma_curr, Xtr_W_X_F)
#         Ftr_Xtr_W_X_Xtr_W_y <- as.vector(crossprod(Xtr_W_X_F, crossprod(X, wts_curr * y)))
#         Ftr_Xtr_W_X_Xtr_W_X_F <- crossprod(Xtr_W_X_F)
#         Ftr_Xtr_W_X_Xtr_W_X_shat <- as.vector(crossprod(Xtr_W_X_F, Xtr_W_X %*% s_hat_reparam))
#         # opt <- lbfgs(fun_reparam_lm_AL, grad_reparam_lm_AL, z_hat_curr, Ftr_Xtr_W_X_gamma = Ftr_Xtr_W_X_gamma, Ftr_Xtr_W_X_Xtr_W_y = Ftr_Xtr_W_X_Xtr_W_y, Ftr_Xtr_W_X_Xtr_W_X_F = Ftr_Xtr_W_X_Xtr_W_X_F, Ftr_Xtr_W_X_Xtr_W_X_shat = Ftr_Xtr_W_X_Xtr_W_X_shat, dual_step = dual_step)
#         opt <- do.call('lbfgs', c(list(fun_reparam_lm_AL, grad_reparam_lm_AL, z_hat_curr, 
#                                        Ftr_Xtr_W_X_gamma = Ftr_Xtr_W_X_gamma, Ftr_Xtr_W_X_Xtr_W_y = Ftr_Xtr_W_X_Xtr_W_y, 
#                                        Ftr_Xtr_W_X_Xtr_W_X_F = Ftr_Xtr_W_X_Xtr_W_X_F, 
#                                        Ftr_Xtr_W_X_Xtr_W_X_shat = Ftr_Xtr_W_X_Xtr_W_X_shat, 
#                                        dual_step = dual_step), lbfgs_arg), quote = T)
#         z_hat_new <- as.vector(opt$par)
#         if (verbose) {
#           if (opt$convergence < 0) {
#             warning(paste('Errors in execution of lbfgs at iteration', j))
#           }
#         }
#       }
#       beta_new <- as.vector(F_reparam %*% z_hat_new + s_hat_reparam)
#       if (beta_diff_flag) {
#         beta_diff <- max(abs(beta_curr - beta_new)) / max(abs(beta_curr))
#       }
#       z_hat_curr <- z_hat_new
#       beta_curr <- beta_new
#       R_new <- as.vector(y - X %*% beta_curr)
#       # R_new <- Diagonal(x = as.numeric(y - X %*% beta_curr))
#       Xtr_R <- t(R_new * X)
#       # Xtr_R <- crossprod(X, R_curr)
#       ##update fun and grad arguments
#       fun_arg$Xtr_R <- Xtr_R
#       grad_arg$Xtr_R <- Xtr_R
#       ##
#       #######
#     }
#     ##update R
#     R_diff <- (R_new - R_curr)
#     R_curr <- R_new
#     ##
#     #######dual residual
#     if (primal_dual_flag) {
#       # dual_resid <- sum((dual_step * R_curr %*% tcrossprod(X) %*% ((R_new - R_curr) %*% wts_curr))^2)
#       if (dual_type == 1) {
#         dual_resid <- sqrt(sum((dual_step * tcrossprod(R_curr * X, R_diff * X) %*% wts_curr)^2))
#       } else if (dual_type == 2) {
#         dual_resid <- sqrt(sum(R_diff * X %*% (gamma_curr + dual_step * crossprod(R_diff * X, wts_curr)))^2)
#       }
#     }
#     #######
#     #######primal residual
#     if (primal_dual_flag) {
#       primal_resid <- sqrt(sum((Xtr_R %*% wts_curr)^2))
#     }
#     #######
#     #######gamma update
#     gamma_new <- gamma_curr + dual_step * Xtr_R %*% wts_curr
#     if (gamma_diff_flag) {
#       gamma_diff <- max(abs(gamma_new - gamma_curr)) / max(abs(gamma_curr))
#     }
#     gamma_curr <- gamma_new
#     #######
#     ####update fun_arg and grad_arg
#     gammatr_Xtr_R <- crossprod(gamma_curr, Xtr_R)
#     fun_arg$gammatr_Xtr_R <- gammatr_Xtr_R
#     grad_arg$R_X_gamma <- crossprod(Xtr_R, gamma_curr)
#     # grad_arg$R_X_gamma <- as.matrix(R_curr %*% (X %*% gamma_curr))
#     ####
#     ###function value update
#     # outer_fval_new <- do.call(fun, c(list(wts = wts_curr), fun_arg))
#     outer_fval_new <- sum(log(wts_curr))
#     outer_fval[j] <- outer_fval_new
#     if (outer_tol_type == 'fval') {
#       outer_fval_diff <- abs(outer_fval_new - outer_fval_curr) / abs(outer_fval_curr)
#     }
#     outer_fval_curr <- outer_fval_new
#     ###
#     #######vary penalty parameter
#     if (vary_penalty == 'RB') {
#       if (primal_resid > RB_mu * dual_resid) {
#         dual_step <- RB_tau * dual_step
#       } else if (dual_resid > RB_mu * primal_resid) {
#         dual_step <- dual_step / RB_tau
#       }
#       ##update fun_arg and grad_arg
#       fun_arg$dual_step <- dual_step
#       grad_arg$dual_step <- dual_step
#       ##
#     }
#     #######
#     #######stopping
#     if (outer_tol_type == 'gamma') {
#       outer_tol <- gamma_diff
#     } else if (outer_tol_type == 'beta') {
#       outer_tol <- beta_diff
#     } else if (outer_tol_type == 'beta_gamma') {
#       outer_tol <- min(beta_diff, gamma_diff)
#     } else if (outer_tol_type == 'fval') {
#       outer_tol <- outer_fval_diff
#     } else if (outer_tol_type %in% c('primal', 'primal_dual')) {
#       outer_tol <- primal_resid
#     }
#     if (outer_tol_type == 'primal_dual') {
#       if (primal_resid < primal_eps && dual_resid < dual_eps) {
#         outer_converged <- T
#         break
#       }
#     } else {
#       if (all(outer_tol < outer_eps) & !any(is.nan(outer_tol) || is.na(outer_tol) || is.infinite(outer_tol))) {
#         outer_converged <- T
#         break
#       }
#     }
#     #######
#   }
#   end_time <- Sys.time()
#   tot_time <- end_time - start_time
#   #################
#   if (!outer_converged) {
#     warning(paste('Outer loop did not converge, current tolerance is', outer_tol))
#   }
#   
#   ####cleaning up output
#   outer_fval <- outer_fval[1:j]
#   inner_nits <- inner_nits[1:j]
#   # logelr <- sum(log(wts_curr))
#   logelr <- outer_fval_curr
#   if (!primal_dual_flag) {
#     primal_resid <- 'not calculated'
#     dual_resid <- 'not calculated'
#   }
#   ####
#   
#   return(list(logelr = logelr, logelr_stat = -2 * logelr, wts = wts_curr / n, sum_wts = sum(wts_curr) / n,
#               gamma = gamma_curr, beta_opt = beta_curr, outer_converged = outer_converged, time = tot_time, 
#               outer_fval = -outer_fval, primal_resid = primal_resid, dual_resid = dual_resid, rho = dual_step,
#               outer_tol = outer_tol, inner_nits = inner_nits, outer_nits = j))
# }
# ##################



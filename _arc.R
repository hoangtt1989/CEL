# 
# 
# ##################function for CEL linear regression admm
# CEL_lm_admm <- function(beta_test, X, y, F_reparam, s_hat_reparam, wts_init = NULL, beta_opt_type = c('closed_form', 'LBFGS'), vary_penalty = c('RB', 'RB2', 'none'), 
#                         RB_mu = 10, RB_tau = 2, RB2_ksi = 2, outer_eps = 1e-8, outer_rel_eps = 1e-4, dual_step = 2, outer_maxit = 1000, dual_type = 1, wts_beta_rep = 1, random_order = F,
#                         outer_tol_type = c('fval', 'primal_dual', 'gamma', 'beta', 'beta_gamma', 'primal'), 
#                         mirror_arg = list(), lbfgs_arg = list(epsilon = 1e-8, invisible = 1), verbose = F) {
#   ##checking inputs
#   beta_opt_type <- beta_opt_type[1]
#   outer_tol_type <- outer_tol_type[1]
#   vary_penalty <- vary_penalty[1]
#   primal_dual_flag <- ifelse(vary_penalty %in% c('RB', 'RB2') || outer_tol_type %in% c('primal', 'primal_dual'), T, F)
#   beta_diff_flag <- ifelse(outer_tol_type == 'beta' || outer_tol_type == 'beta_gamma', T, F)
#   gamma_diff_flag <- ifelse(outer_tol_type == 'gamma' || outer_tol_type == 'beta_gamma', T, F)
#   if (!(outer_tol_type %in% c('gamma', 'beta', 'beta_gamma', 'fval', 'primal', 'primal_dual'))) {
#     stop('supply a valid outer tolerance type')
#   }
#   if (!(vary_penalty %in% c('RB', 'RB2', 'none'))) {
#     stop('supply a valid penalty varying type')
#   }
#   ##
#   ##initializing function and gradient
#   # fun <- fun_lm_AL_cpp
#   # grad <- grad_lm_AL_cpp
#   fun <- cmpfun(fun_lm_AL2) ##compile function
#   grad <- cmpfun(grad_lm_AL2) ##compile function
#   mirror_descent <- cmpfun(mirror_descent)
#   wts_wrapper <- cmpfun(wts_wrapper)
#   beta_wrapper <- cmpfun(beta_wrapper)
#   fun_reparam_lm_AL <- cmpfun(fun_reparam_lm_AL)
#   grad_reparam_lm_AL <- cmpfun(grad_reparam_lm_AL)
#   ##
#   ##initial parameters
#   prob_size <- dim(X)
#   n <- prob_size[1]
#   p <- prob_size[2]
#   beta_curr <- beta_test
#   # z_hat_curr <- as.vector(crossprod(F_reparam, beta_test - s_hat_reparam))
#   if (is.null(wts_init)) {
#     wts_init <- rep(1, n)
#   } else {
#     if (length(wts_init) != n) {
#       stop('wts_init must be a vector with n elements')
#     }
#   }
#   wts_new <- wts_init
#   wts_curr <- wts_new
#   gamma_curr <- rep(0, p)
#   outer_converged <- F
#   ##
#   ##primal_dual stopping
#   if (outer_tol_type == 'primal_dual') {
#     primal_eps <- outer_eps * sqrt(p)
#   }
#   # dual_eps <- outer_eps * sqrt(n) ## dimension of dual residual changes depending on update order
#   ##
#   ##initialize holders
#   outer_fval <- rep(NA, outer_maxit)
#   mirror_nits <- rep(NA, outer_maxit)
#   update_order <- rep(NA, outer_maxit)
#   rho_vec <- rep(NA, outer_maxit)
#   ##
#   ####precomputing
#   R_new <- as.vector(y - X %*% beta_curr)
#   R_curr <- R_new
#   Xtr_R <- t(R_new * X) #t(X) %*% R_curr
#   ###initialize function/gradient arguments
#   fun_grad_arg <- list(gamma = gamma_curr, X = X, R = R_curr, dual_step = dual_step)
#   ##
#   # outer_fval_curr <- do.call(fun, c(list(wts = wts_curr), fun_arg))
#   outer_fval_curr <- sum(log(wts_new))
#   ####
#   #################main loop
#   start_time <- Sys.time()
#   for (j in 1:outer_maxit) {
#     if (verbose) {
#       if (! j %% 50) {
#         message(paste('Iteration', j))
#       }
#     }
#     ##random ordering of wts/beta update
#     if (!random_order || runif(1) < .5) {
#       wts_beta_order <- T ##update wts then beta
#       dual_resid_elt <- sqrt(n) ##dual residual has n elements
#       update_order[j] <- 'wts_beta'
#     } else {
#       wts_beta_order <- F ##update beta then wts
#       dual_resid_elt <- sqrt(p) ##dual residual has p elements
#       update_order[j] <- 'beta_wts'
#     }
#     if (outer_tol_type == 'primal_dual') {
#       dual_eps <- outer_eps * dual_resid_elt
#     }
#     ##
#     ##wts and beta updates can be repeated
#     for (wts_beta_iter in 1:wts_beta_rep) {
#       if (wts_beta_order) {
#         #######mirror descent
#         A_max <- dual_step * max(abs(crossprod(Xtr_R))) * n ##t(Xtr_R) %*% Xtr_R
#         mirror_step <- 1 / (A_max + max(abs(crossprod(gamma_curr, Xtr_R))))
#         mirror_res <- do.call('mirror_descent', c(list(fun = fun, grad = grad, wts = wts_new, fun_arg = fun_grad_arg, grad_arg = fun_grad_arg, mirror_step = mirror_step), mirror_arg), quote = T)
#         wts_new <- as.vector(mirror_res$wts)
#         # wts_res <- wts_wrapper(n, fun, grad, wts_new, fun_arg, grad_arg, mirror_arg)
#         # wts_new <- wts_res$wts_new
#         # mirror_iter <- wts_res$nits
#         # mirror_converged <- wts_res$converged
#         #######
#         #######beta update
#         beta_res <- beta_wrapper(beta_curr, wts_new, X, F_reparam, y, gamma_curr, dual_step, s_hat_reparam, beta_opt_type, lbfgs_arg, fun_reparam_lm_AL, grad_reparam_lm_AL)
#         beta_new <- beta_res$beta_new
#         R_new <- beta_res$R_new
#         Xtr_R <- t(R_new * X)
#         # Xtr_R <- beta_res$Xtr_R
#         beta_converged <- beta_res$converged
#         ##update fun and grad arguments
#         fun_grad_arg$R <- R_new
#         ##
#         #######
#       } else {
#         #######beta update
#         beta_res <- beta_wrapper(beta_curr, wts_new, X, F_reparam, y, gamma_curr, dual_step, s_hat_reparam, beta_opt_type, lbfgs_arg, fun_reparam_lm_AL, grad_reparam_lm_AL)
#         beta_new <- beta_res$beta_new
#         R_new <- beta_res$R_new
#         Xtr_R <- t(R_new * X)
#         # Xtr_R <- beta_res$Xtr_R
#         beta_converged <- beta_res$converged
#         ##update fun and grad arguments
#         fun_grad_arg$R <- R_new
#         # fun_arg$Xtr_R <- Xtr_R
#         # grad_arg$Xtr_R <- Xtr_R
#         ##
#         #######
#         #######mirror descent
#         A_max <- dual_step * max(abs(crossprod(Xtr_R))) * n ##t(Xtr_R) %*% Xtr_R
#         mirror_step <- 1 / (A_max + max(abs(crossprod(gamma_curr, Xtr_R))))
#         mirror_res <- do.call('mirror_descent', c(list(fun = fun, grad = grad, wts = wts_new, fun_arg = fun_grad_arg, grad_arg = fun_grad_arg, mirror_step = mirror_step), mirror_arg), quote = T)
#         wts_new <- as.vector(mirror_res$wts)
#         # wts_res <- wts_wrapper(n, fun, grad, wts_new, fun_arg, grad_arg, mirror_arg)
#         # wts_new <- wts_res$wts_new
#         # mirror_iter <- wts_res$nits
#         # mirror_converged <- wts_res$converged
#         #######
#       }
#     }
#     ##
#     mirror_nits[j] <- mirror_res$iter
#     if (verbose) {
#       if (!mirror_res$converged) {
#         warning(paste('Mirror descent did not converge at iteration', j))
#       }
#       if (!beta_res$converged) {
#         warning(paste('Errors in execution of lbfgs at iteration', j))
#       }
#     }
#     ##update beta_diff
#     if (beta_diff_flag) {
#       beta_sc <- max(abs(beta_curr))
#       beta_sc <- ifelse(beta_sc < .Machine$double.eps, .Machine$double.eps * 2, beta_sc)
#       beta_diff <- max(abs(beta_curr - beta_new)) / beta_sc
#     }
#     ##
#     ##update R_diff and wts_diff for dual residual
#     if (wts_beta_order) { ##wts were updated first
#       R_diff <- as.vector(R_new - R_curr)
#     } else { ##beta was updated first
#       wts_diff <- as.vector(wts_new - wts_curr)
#     }
#     ##
#     #######gamma update
#     gamma_new <- gamma_curr + dual_step * Xtr_R %*% wts_new
#     if (gamma_diff_flag) {
#       gamma_sc <- max(abs(gamma_curr))
#       gamma_sc <- ifelse(gamma_sc < .Machine$double.eps, .Machine$double.eps * 2, gamma_sc)
#       gamma_diff <- max(abs(gamma_new - gamma_curr)) / gamma_sc
#     }
#     #######
#     #######dual residual
#     if (primal_dual_flag) {
#       if (wts_beta_order) { ##wts were updated first
#         # dual_resid <- sum((dual_step * R_curr %*% tcrossprod(X) %*% ((R_new - R_curr) %*% wts_curr))^2)
#         dual_resid <- sqrt(sum((dual_step * tcrossprod(R_curr * X, R_diff * X) %*% wts_new)^2))
#         dual_resid_scale <- sqrt(sum((R_curr * X %*% gamma_new)^2)) 
#       } else { ##beta was updated first
#         # dual_resid <- sqrt(sum(dual_step * tcrossprod(X, wts_curr * X) %*% (tcrossprod(X, wts_diff * y) - tcrossprod(X, wts_diff * X * beta_new)))^2)
#         dual_resid <- sqrt(sum((dual_step * crossprod(X, wts_curr * tcrossprod(X)) %*% (wts_diff * y - wts_diff * (X %*% beta_new)))^2))
#       }
#     }
#     #######
#     #######primal residual
#     if (primal_dual_flag) {
#       primal_resid <- sqrt(sum((Xtr_R %*% wts_new)^2))
#       wts_X <- wts_new * X
#       primal_resid_scale <- max(sqrt(sum((crossprod(wts_X, y))^2)), sqrt(sum((crossprod(wts_X, X %*% beta_curr))^2)))
#     }
#     #######
#     ##update R, gamma, beta and wts
#     beta_curr <- as.vector(beta_new)
#     R_curr <- as.vector(R_new)
#     wts_curr <- as.vector(wts_new)
#     gamma_curr <- gamma_new
#     ##
#     ####update fun_arg and grad_arg
#     fun_grad_arg$gamma <- gamma_curr
#     ####
#     ###function value update
#     # outer_fval_new <- do.call(fun, c(list(wts = wts_curr), fun_arg))
#     outer_fval_new <- sum(log(wts_curr))
#     outer_fval[j] <- outer_fval_new
#     if (outer_tol_type == 'fval') {
#       outer_fval_sc <- max(abs(outer_fval_curr))
#       outer_fval_sc <- ifelse(outer_fval_sc < .Machine$double.eps, .Machine$double.eps * 2, outer_fval_sc)
#       outer_fval_diff <- abs(outer_fval_new - outer_fval_curr) / abs(outer_fval_sc)
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
#       fun_grad_arg$dual_step <- dual_step
#       ##
#     }
#     if (vary_penalty == 'RB2') {
#       ###varying tau
#       primal_dual_div <- primal_resid / dual_resid
#       dual_primal_div <- 1 / primal_dual_div
#       sqrt_primal_dual_div <- sqrt(primal_dual_div / RB2_ksi)
#       sqrt_dual_primal_div <- sqrt(dual_primal_div * RB2_ksi)
#       if (sqrt_primal_dual_div < RB_tau && sqrt_primal_dual_div >= 1) {
#         RB2_tau <- sqrt_primal_dual_div
#       } else if (sqrt_dual_primal_div < 1 && sqrt_dual_primal_div > (1 / RB_tau)) {
#         RB2_tau <- sqrt_primal_dual_div
#       } else {
#         RB2_tau <- RB_tau
#       }
#       ###
#       primal_resid_rel <- primal_resid / primal_resid_scale
#       dual_resid_rel <- dual_resid / dual_resid_scale
#       if (primal_resid_rel > RB2_ksi * RB_mu * dual_resid_rel) {
#         dual_step <- RB2_tau * dual_step
#       } else if (dual_resid_rel > RB_mu * primal_resid_rel / RB2_ksi) {
#         dual_step <- dual_step / RB2_tau
#       }
#       ##update fun_arg and grad_arg
#       fun_grad_arg$dual_step <- dual_step
#       ##
#     }
#     rho_vec[j] <- dual_step
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
#       primal_eps_scale <- primal_eps + outer_rel_eps * primal_resid_scale / n
#       dual_eps_scale <- dual_eps + outer_rel_eps * dual_resid_scale
#       if (primal_resid / n < primal_eps_scale && dual_resid < dual_eps_scale) {
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
#   rho_vec <- rho_vec[1:j]
#   outer_fval <- outer_fval[1:j]
#   mirror_nits <- mirror_nits[1:j]
#   update_order <- update_order[1:j]
#   # logelr <- sum(log(wts_curr))
#   logelr <- outer_fval_curr
#   if (!primal_dual_flag) {
#     primal_resid <- 'not calculated'
#     dual_resid <- 'not calculated'
#   }
#   ####
#   
#   return(list(logelr = logelr, logelr_stat = -2 * logelr, wts = wts_curr / n, sum_wts = sum(wts_curr) / n,
#               gamma = gamma_curr, beta_opt = beta_curr, outer_converged = outer_converged, time = as.double(tot_time, 'secs'), 
#               outer_fval = -outer_fval, primal_resid = primal_resid, dual_resid = dual_resid, rho = rho_vec,
#               outer_tol = outer_tol, update_order = update_order, mirror_nits = mirror_nits, outer_nits = j))
# }
# ##################


# ##################function for CEL linear regression nested - use lbfgs library but doesn't work because we have to save updates in outer loop
# CEL_lm_nested <- function(beta_test, X, y, F_reparam, s_hat_reparam,
#                           outer_eps = 1e-8, outer_maxit = 1000, beta_opt_type = c('damped_newton', 'LBFGS'), newton_arg = list(), lbfgs_arg = list(epsilon = 1e-8, invisible = 1),
#                           outer_tol_type = c('fval', 'gval'), inner_arg = list(), verbose = F) {
#   ##checking inputs
#   beta_opt_type <- beta_opt_type[1]
#   outer_tol_type <- outer_tol_type[1]
#   if (!(outer_tol_type %in% c('fval', 'gval'))) {
#     stop('supply a valid outer tolerance type')
#   }
#   ##
#   ##compile owen's functions
#   owen_EL <- cmpfun(elm) ##inner loop
#   # plog0 <- cmpfun(llog) ##pseudo log
#   # plog1 <- cmpfun(llogp) ##pseudo log first deriv
#   # plog2 <- cmpfun(llogpp) ##pseudo log second deriv
#   # owen_logelr <- cmpfun(owen_logelr_fun) ##logelr
#   ##
#   ##compile gradient/hessian functions
#   lm_hess_z <- cmpfun(lm_hess_z)
#   lm_part_z <- cmpfun(lm_part_z)
#   lm_fval_z <- cmpfun(lm_fval_z)
#   ##
#   ##initial parameters
#   prob_size <- dim(X)
#   n <- prob_size[1]
#   p <- prob_size[2]
#   beta_curr <- beta_test
#   z_reparam_curr <- crossprod(F_reparam, beta_test - s_hat_reparam)
#   R_curr <- as.vector(y - X %*% beta_test)
#   outer_converged <- F
#   ##
#   ##initialize holders
#   outer_fval <- rep(NA, outer_maxit)
#   inner_nits <- rep(NA, outer_maxit)
#   ##
#   ##first inner loop iteration
#   score_curr <- R_curr * X
#   inner_res <- do.call('owen_EL', c(list(score_curr), inner_arg), quote = T)
#   gamma_curr <- inner_res$lambda
#   logler_curr <- sum(llog(inner_res$wts, 1/n))
#   ##
#   ##set up function, gradient, hessian arguments
#   fun_arg <- list(X = X, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, y = y, gamma = gamma_curr)
#   grad_arg <- list(X = X, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, y = y, gamma = gamma_curr)
#   hess_arg <- list(X = X, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, y = y, gamma = gamma_curr)
#   ##
#   ##initial parameters for lbfgs
#   # if (beta_opt_type == 'LBFGS') {
#   #   m <- lbfgs_arg$m
#   #   num_par <- length(z_reparam_curr)
#   #   mem_y <- matrix(0, nrow = num_par, ncol = m)
#   #   mem_s <- matrix(0, nrow = num_par, ncol = m)
#   #   mem_rho <- rep(0, m)
#   #   # mem_y[, m] <- do.call('lm_part_z', c(list(z_reparam_curr), grad_arg), quote = T)
#   #   # mem_s[, m] <- z_reparam_curr
#   #   # mem_rho[m] <- 1 / as.numeric(mem_y[, m] %*% mem_s[, m])
#   #   mem_y[, m] <- rep(0, num_par)
#   #   mem_s[, m] <- rep(0, num_par)
#   #   mem_rho[m] <- 0
#   #   grad_curr <- mem_y[, 1]
#   #   x_curr <- mem_s[, 1]
#   #   # hess_init <- diag(num_par) * as.numeric(mem_s[, m] %*% mem_y[, m]) / as.numeric(mem_y[, m] %*% mem_y[, m])
#   #   hess_init <- diag(num_par)
#   # }
#   ##
#   start_time <- Sys.time()
#   for (j in 1:outer_maxit) {
#     if (verbose) {
#       if (! j %% 5) {
#         message(paste('Iteration', j))
#       }
#     }
#     ##beta optimization
#     if (beta_opt_type == 'damped_newton') {
#       damped_res <- do.call('damped_newton_step', c(list(z_reparam_curr, lm_fval_z, lm_part_z, lm_hess_z, fun_arg, grad_arg, hess_arg), newton_arg), quote = T)
#       if (!damped_res$ls_converged) {
#         warning(paste('Line search unsuccessful at j =', j))
#       }
#       x_curr <- damped_res$x_curr
#       fval_curr <- damped_res$fval_curr
#     } else if (beta_opt_type == 'LBFGS') {
#       lbfgs_res <- do.call('lbfgs', c(list(lm_fval_z, lm_part_z, z_reparam_curr, X = X, F_reparam = F_reparam, 
#                                            s_hat_reparam = s_hat_reparam, y = y, 
#                                            gamma = gamma_curr, max_iterations = 2), lbfgs_arg), quote = T)
#       x_curr <- lbfgs_res$par
#       fval_curr <- do.call('lm_fval_z', c(list(x_curr), fun_arg), quote = T)
#     }
#     ##
#     ##update values
#     z_reparam_curr <- x_curr
#     beta_curr <- F_reparam %*% z_reparam_curr + s_hat_reparam
#     R_new <- as.vector(y - X %*% beta_curr)
#     score_new <- R_new * X
#     outer_fval_curr <- fval_curr
#     # outer_fval_new <- damped_res$fval_new
#     ##
#     ####inner loop
#     inner_res <- do.call('owen_EL', c(list(score_new), inner_arg), quote = T)
#     outer_fval_new <- -sum(llog(inner_res$wts, 1/n))
#     outer_fval[j] <- outer_fval_new
#     inner_nits[j] <- inner_res$nits
#     gamma_curr <- inner_res$lambda
#     fun_arg$gamma <- gamma_curr
#     grad_arg$gamma <- gamma_curr
#     hess_arg$gamma <- gamma_curr
#     ####
#     #####stopping
#     outer_tol <- switch(outer_tol_type,
#                         fval = abs(outer_fval_new - outer_fval_curr)/abs(outer_fval_curr),
#                         gval = sum(damped_res$gval_new^2))
#     if (all(outer_tol < outer_eps) & !any(is.nan(outer_tol) || is.na(outer_tol) || is.infinite(outer_tol))) {
#       outer_converged <- T
#       break
#     }
#     #####
#   }
#   end_time <- Sys.time()
#   tot_time <- end_time - start_time
#   if (!outer_converged) {
#     warning(paste('Outer loop did not converge, current tolerance is', outer_tol))
#   }
#   
#   ####cleaning up output
#   outer_fval <- outer_fval[1:j]
#   inner_nits <- inner_nits[1:j]
#   logelr <- -outer_fval_new
#   wts <- inner_res$wts/n
#   sum_wts <- sum(wts)
#   if (abs(sum_wts - 1) > 1e-4) { ##output warning if weights don't sum to 1
#     warning('weights do not sum to 1 - convex hull constraint may not be satisfied')
#   }
#   ####
#   
#   return(list(logelr = logelr, logelr_stat = -2 * logelr, wts = wts, sum_wts = sum_wts,
#               gamma = gamma_curr, beta_opt = beta_curr, outer_converged = outer_converged, time = tot_time, 
#               outer_fval = outer_fval, outer_tol = outer_tol, inner_nits = inner_nits, outer_nits = j))
# }
# ##################



# ##################function for CEL linear regression nested - damped newton written inside loop
# CEL_lm_nested <- function(beta_test, X, y, F_reparam, s_hat_reparam,
#                           outer_eps = 1e-8, outer_maxit = 1000, ls_maxit = 100, ls_eps = 1e-12, ls_beta = .7, ls_alpha = .25,
#                           outer_tol_type = c('fval', 'gval'), inner_arg = list(), verbose = F) {
#   ##checking inputs
#   outer_tol_type <- outer_tol_type[1]
#   if (!(outer_tol_type %in% c('fval', 'gval'))) {
#     stop('supply a valid outer tolerance type')
#   }
#   ##
#   ##compile owen's functions
#   owen_EL <- cmpfun(elm) ##inner loop
#   plog0 <- cmpfun(llog) ##pseudo log
#   plog1 <- cmpfun(llogp) ##pseudo log first deriv
#   plog2 <- cmpfun(llogpp) ##pseudo log second deriv
#   # owen_logelr <- cmpfun(owen_logelr_fun) ##logelr
#   ##
#   ##compile gradient/hessian functions
#   lm_part_z <- cmpfun(lm_part_z)
#   lm_part_z_gamma <- cmpfun(lm_part_z_gamma)
#   lm_part_z_z <- cmpfun(lm_part_z_z)
#   ##
#   ##initial parameters
#   prob_size <- dim(X)
#   n <- prob_size[1]
#   p <- prob_size[2]
#   beta_curr <- beta_test
#   z_reparam_curr <- crossprod(F_reparam, beta_test - s_hat_reparam)
#   R_curr <- as.vector(y - X %*% beta_test)
#   outer_converged <- F
#   ls_converged <- F
#   ls_diff <- -1
#   ##
#   ##initialize holders
#   outer_fval <- rep(NA, outer_maxit)
#   inner_nits <- rep(NA, outer_maxit)
#   ##
#   ##first inner loop iteration
#   score_curr <- R_curr * X
#   inner_res <- owen_EL(score_curr)
#   gamma_curr <- inner_res$lambda
#   ##
#   start_time <- Sys.time()
#   for (j in 1:outer_maxit) {
#     if (verbose) {
#       if (! j %% 5) {
#         message(paste('Iteration', j))
#       }
#     }
#     ##precomputing
#     R_curr <- as.vector(y - X %*% beta_curr)
#     score_curr <- R_curr * X
#     ##
#     ##current function value
#     outer_fval_curr <- sum(plog0(1 + score_curr %*% gamma_curr, 1/n)) ##current function value
#     ##
#     ##gradient
#     grad <- lm_part_z( z_reparam_curr, X, F_reparam, s_hat_reparam, y, gamma_curr, plog1)
#     ##
#     ###hessian calculations
#     L_z_gamma <- lm_part_z_gamma(z_reparam_curr, X, F_reparam, s_hat_reparam, y, gamma_curr, plog1, plog2)
#     L_gamma_z <- t(L_z_gamma)
#     L_gamma_gamma <- crossprod(score_curr, as.vector(plog2(1 + score_curr %*% gamma_curr, 1/n)) * score_curr)
#     L_z_z <- lm_part_z_z(z_reparam_curr, X, F_reparam, s_hat_reparam, y, gamma_curr, plog2)
#     ###
#     ##hessian
#     hess <- L_z_z - L_z_gamma %*% solve(L_gamma_gamma, L_gamma_z)
#     ##
#     #####update search direction
#     pk <- -solve(hess, grad)
#     gradtr_pk <- crossprod(grad, pk)
#     z_reparam_new <- z_reparam_curr + pk
#     beta_new <- F_reparam %*% z_reparam_new + s_hat_reparam
#     R_new <- as.vector(y - X %*% beta_new)
#     score_new <- R_new * X
#     #####
#     #########line search
#     ls_fval <- sum(plog0(1 + score_new %*% gamma_curr, 1/n))
#     step <- 1
#     ls_converged <- ls_fval < (outer_fval_curr + ls_alpha * step * gradtr_pk)
#     if (!ls_converged) {
#       for (ls_iter in 1:ls_maxit) {
#         step <- ls_beta * step
#         z_reparam_new <- z_reparam_curr + step * pk
#         beta_new <- F_reparam %*% z_reparam_new + s_hat_reparam
#         R_new <- as.vector(y - X %*% beta_new)
#         score_new <- R_new * X
#         ls_fval <- sum(plog0(1 + score_new %*% gamma_curr, 1/n))
#         ls_diff <- outer_fval_curr + ls_alpha * step * gradtr_pk - ls_fval
#         if (ls_diff >= 0 | abs(ls_diff) < ls_eps) {
#           ls_converged <- T
#           break
#         }
#       }
#       if (!ls_converged) {
#         warning(paste('Line search unsuccessful at j =', j))
#       }
#     }
#     #########
#     ##update values
#     beta_curr <- beta_new
#     z_reparam_curr <- z_reparam_new
#     outer_fval_new <- ls_fval
#     outer_fval[j] <- outer_fval_new
#     ##
#     ####inner loop
#     inner_res <- elm(score_new)
#     inner_nits[j] <- inner_res$nits
#     gamma_curr <- inner_res$lambda
#     ####
#     #####stopping
#     outer_tol <- switch(outer_tol_type,
#                         fval = abs(outer_fval_new - outer_fval_curr)/abs(outer_fval_curr),
#                         gval = sum(grad^2))
#     if (all(outer_tol < outer_eps) & !any(is.nan(outer_tol) || is.na(outer_tol) || is.infinite(outer_tol))) {
#       outer_converged <- T
#       break
#     }
#     #####
#   }
#   end_time <- Sys.time()
#   tot_time <- end_time - start_time
#   if (!outer_converged) {
#     warning(paste('Outer loop did not converge, current tolerance is', outer_tol))
#   }
#   
#   ####cleaning up output
#   outer_fval <- outer_fval[1:j]
#   inner_nits <- inner_nits[1:j]
#   logelr <- sum(plog0(inner_res$wts, 1/n))
#   wts <- inner_res$wts/n
#   sum_wts <- sum(wts)
#   if (abs(sum_wts - 1) > 1e-6) { ##output warning if weights don't sum to 1
#     warning('weights do not sum to 1 - convex hull constraint may not be satisfied')
#   }
#   ####
#   
#   return(list(logelr = logelr, logelr_stat = -2 * logelr, wts = wts, sum_wts = sum_wts,
#               gamma = gamma_curr, beta_opt = beta_curr, outer_converged = outer_converged, time = tot_time, 
#               outer_fval = -outer_fval, outer_tol = outer_tol, inner_nits = inner_nits, outer_nits = j))
# }
# ##################


# ##################profile a component of the beta vector for linear regression using ADMM
# emplik_lm_admm_profiler <- function(alpha_add, fix_index, X, y, lm_admm_arg = list(), mirror_arg = list(), verbose = F) {
#   ##initial parameters
#   prob_size <- dim(X)
#   p <- prob_size[2]
#   beta_OLS <- solve(crossprod(X), crossprod(X, y))
#   ##
#   ##set up the reparam
#   T_reparam <- t(as.matrix(rep(0, p)))
#   T_reparam[fix_index] <- 1
#   alpha_test <- beta_OLS[fix_index] + alpha_add
#   s_hat_reparam <- as.matrix(rep(0, p))
#   s_hat_reparam[fix_index] <- alpha_test
#   F_reparam <- orth(nullspace(T_reparam))
#   ##
#   ##set up the test
#   beta_test <- beta_OLS
#   beta_test[fix_index] <- alpha_test
#   ##
#   ##run admm
#   admm_cmp <- cmpfun(emplik_lm_admm)
#   admm_res <- do.call(admm_cmp, c(list(beta_test = beta_test, X = X, y = y, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, mirror_arg = mirror_arg, verbose = verbose), lm_admm_arg), quote = T)
#   ##
#   return(c(list(beta_OLS = beta_OLS), admm_res))
# }
# ##################

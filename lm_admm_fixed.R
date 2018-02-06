
##################function for EL (beta given) linear regression with outliers
EL_lm_AL_outlier <- function(beta_test, X, y, dual_step = .5, outer_eps = 1e-8, outer_maxit = 1000, vary_penalty = F, RB_mu = 10, RB_tau = 2,
                             outer_tol_type = c('primal', 'primal_dual', 'fval', 'gamma'), mirror_arg = list(), detect_outlier = T, q = NULL,
                             outlier_loop_arg = list(), verbose = F) {
  ##checking inputs
  outer_tol_type <- outer_tol_type[1]
  if (!(outer_tol_type %in% c('gamma', 'fval', 'primal', 'primal_dual'))) {
    stop('supply a valid outer tolerance type')
  }
  ##
  ##initializing function and gradient
  # fun <- cmpfun(fun_lm_AL) ##compile function
  # grad <- cmpfun(grad_lm_AL) ##compile function
  fun <- cmpfun(fun_lm_AL_cpp)
  grad <- cmpfun(grad_lm_AL_cpp)
  mirror_descent_cmp <- cmpfun(mirror_descent)
  ##
  ##initial parameters
  prob_size <- dim(X)
  n <- prob_size[1]
  p <- prob_size[2]
  beta_curr <- beta_test
  wts_curr <- as.matrix(rep(1, n))
  gamma_curr <- as.matrix(rep(0, p))
  outer_converged <- F
  ##
  ##primal_dual stopping
  primal_eps <- outer_eps * sqrt(p)
  dual_eps <- outer_eps * sqrt(n)
  ##
  ##initialize holders
  outer_fval <- rep(NA, outer_maxit)
  mirror_nits <- rep(NA, outer_maxit)
  primal_resids <- rep(NA, outer_maxit)
  beta_outlier_nits <- rep(NA, outer_maxit)
  beta_outlier_tol <- rep(NA, outer_maxit)
  ##
  ####precomputing
  R_new <- as.vector(y - X %*% beta_curr)
  R_curr <- R_new
  Xtr_R <- t(R_new * X) #t(X) %*% R_curr
  ###initialize function/gradient arguments
  fun_grad_arg <- list(gamma = gamma_curr, X = X, R = R_curr, dual_step = dual_step)
  ##
  # outer_fval_curr <- do.call(fun, c(list(wts = wts_curr), fun_arg))
  outer_fval_curr <- sum(log(wts_curr))
  primal_resid_curr <- sum((Xtr_R %*% wts_curr)^2)
  ####
  
  #################main loop
  start_time <- Sys.time()
  for (j in 1:outer_maxit) {
    if (verbose) {
      if (! j %% 50) {
        message(paste('Iteration', j))
      }
    }
    #######mirror descent
    A_max <- dual_step * max(abs(crossprod(Xtr_R))) * n ##t(Xtr_R) %*% Xtr_R
    mirror_step <- 1 / (A_max + max(abs(crossprod(gamma_curr, Xtr_R))))
    mirror_res <- do.call('mirror_descent_lm_AL_cpp', c(list(wts = wts_curr, mirror_step = mirror_step), fun_grad_arg, mirror_arg), quote = T)
    wts_new <- as.vector(mirror_res$wts)
    wts_curr <- wts_new
    mirror_nits[j] <- mirror_res$iter
    #######
    #######outliers
    beta_outlier_res <- outlier_loop(wts_new, X, F_reparam = diag(p), y, gamma_curr, dual_step, s_hat_reparam = rep(0, p), beta_opt_type = 'closed_form', 
                                     lbfgs_arg = NULL, fun_reparam_lm_AL, grad_reparam_lm_AL)
    delta_new <- beta_outlier_res$delta_new
    beta_new <- beta_outlier_res$beta_new
    R_new <- beta_outlier_res$R_adj_new
    Xtr_R <- t(R_new * X)
    # Xtr_R <- beta_res$Xtr_R
    beta_converged <- beta_outlier_res$converged
    ##update fun and grad arguments
    fun_grad_arg$R <- R_new
    y_adj <- beta_outlier_res$y_adj
    beta_outlier_nits[j] <- beta_outlier_res$iter
    beta_outlier_tol[j] <- beta_outlier_res$tol
    #######
    #######gamma update
    Xtr_R_wts <- Xtr_R %*% wts_curr
    gamma_new <- gamma_curr + dual_step * Xtr_R_wts
    gamma_sc <- max(abs(gamma_curr))
    gamma_sc <- ifelse(gamma_sc < .Machine$double.eps, .Machine$double.eps * 2, gamma_sc)
    gamma_diff <- max(abs(gamma_new - gamma_curr)) / gamma_sc
    gamma_curr <- gamma_new
    ####update fun_arg and grad_arg
    fun_grad_arg$gamma <- gamma_curr
    ####
    #######
    #######function value
    # outer_fval_new <- do.call(fun, c(list(wts = wts_curr), fun_arg))
    outer_fval_new <- sum(log(wts_curr))
    outer_fval_sc <- abs(outer_fval_curr)
    outer_fval_sc <- ifelse(outer_fval_sc < .Machine$double.eps, .Machine$double.eps * 2, outer_fval_sc)
    outer_fval_diff <- abs(outer_fval_new - outer_fval_curr) / outer_fval_sc
    outer_fval_curr <- outer_fval_new
    outer_fval[j] <- outer_fval_curr
    #######
    #######primal residual
    primal_resid <- sqrt(sum((Xtr_R_wts)^2))
    # dual_resid <- sqrt(sum((-1/wts_curr + R_X_gamma)^2))
    # if (vary_penalty) {
    #   if (primal_resid > RB_mu * dual_resid) {
    #     dual_step <- RB_tau * dual_step
    #   } else if (dual_resid > RB_mu * primal_resid) {
    #     dual_step <- dual_step / RB_tau
    #   }
    #   ##update fun_arg and grad_arg
    #   fun_arg$dual_step <- dual_step
    #   grad_arg$dual_step <- dual_step
    #   ##
    # }
    # # print(primal_resid_new)
    # print(sqrt(sum((-1/wts_curr + R_X_gamma)^2)))
    # primal_resids[j] <- primal_resid
    #######
    #######stopping
    outer_tol <- switch(outer_tol_type,
                        gamma = gamma_diff,
                        fval = outer_fval_diff,
                        primal = primal_resid_curr,
                        primal_dual = primal_resid_curr)
    if (outer_tol_type == 'primal_dual') {
      if (primal_resid < primal_eps && dual_resid < dual_eps) {
        outer_converged <- T
        break
      }
    }
    else {
      if (outer_tol < outer_eps) {
        outer_converged = T
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
  outer_fval <- outer_fval[1:j]
  mirror_nits <- mirror_nits[1:j]
  beta_outlier_nits <- beta_outlier_nits[1:j]
  beta_outlier_tol <- beta_outlier_tol[1:j]
  logelr <- sum(log(wts_curr))
  ####
  
  return(list(logelr = logelr, logelr_stat = -2 * logelr, wts = wts_curr / n, sum_wts = sum(wts_curr) / n,
              gamma = gamma_curr, outer_converged = outer_converged, time = tot_time, 
              outer_fval = outer_fval, primal_resid = primal_resid, rho = dual_step,
              outer_tol = outer_tol, mirror_nits = mirror_nits, beta_outlier_nits = beta_outlier_nits,
              beta_outlier_tol = beta_outlier_tol, delta_opt = delta_new, beta_opt = beta_new, outer_nits = j))
}
##################


##################function for EL (beta given) linear regression
EL_lm_AL <- function(beta_test, X, y, dual_step = .5, outer_eps = 1e-8, outer_maxit = 1000, vary_penalty = T, RB_mu = 10, RB_tau = 2,
                     outer_tol_type = c('primal', 'primal_dual', 'fval', 'gamma'), mirror_arg = list(), verbose = F) {
  ##checking inputs
  outer_tol_type <- outer_tol_type[1]
  if (!(outer_tol_type %in% c('gamma', 'fval', 'primal', 'primal_dual'))) {
    stop('supply a valid outer tolerance type')
  }
  ##
  ##initializing function and gradient
  # fun <- cmpfun(fun_lm_AL) ##compile function
  # grad <- cmpfun(grad_lm_AL) ##compile function
  fun <- cmpfun(fun_lm_AL_cpp)
  grad <- cmpfun(grad_lm_AL_cpp)
  mirror_descent_cmp <- cmpfun(mirror_descent)
  ##
  ##initial parameters
  prob_size <- dim(X)
  n <- prob_size[1]
  p <- prob_size[2]
  beta_curr <- beta_test
  wts_curr <- as.matrix(rep(1, n))
  gamma_curr <- as.matrix(rep(0, p))
  outer_converged <- F
  ##
  ##primal_dual stopping
  primal_eps <- outer_eps * sqrt(p)
  dual_eps <- outer_eps * sqrt(n)
  ##
  ##initialize holders
  outer_fval <- rep(NA, outer_maxit)
  inner_nits <- rep(NA, outer_maxit)
  primal_resids <- rep(NA, outer_maxit)
  ##
  ####precomputing
  R_curr <- as.vector(y - X %*% beta_curr)
  Xtr_R <- t(R_curr * X) #t(X) %*% R_curr
  A_max <- max(abs(crossprod(Xtr_R))) * n ##t(Xtr_R) %*% Xtr_R
  gammatr_Xtr_R <- crossprod(gamma_curr, Xtr_R) #t(gamma_curr) %*% Xtr_R
  R_X_gamma <- crossprod(Xtr_R, gamma_curr) #t(Xtr_R) %*% gamma_curr
  ###initialize function/gradient arguments
  fun_arg <- list(gamma = gamma_curr, X = X, R = R_curr, dual_step = dual_step)
  grad_arg <- fun_arg
  # fun_arg <- list(gammatr_Xtr_R = gammatr_Xtr_R, Xtr_R = Xtr_R, dual_step = dual_step)
  # grad_arg <- list(R_X_gamma = R_X_gamma, Xtr_R = Xtr_R, dual_step = dual_step)
  ##
  # outer_fval_curr <- do.call(fun, c(list(wts = wts_curr), fun_arg))
  outer_fval_curr <- sum(log(wts_curr))
  primal_resid_curr <- sum((Xtr_R %*% wts_curr)^2)
  ####
  
  #################main loop
  start_time <- Sys.time()
  for (j in 1:outer_maxit) {
    if (verbose) {
      if (! j %% 50) {
        message(paste('Iteration', j))
      }
    }
    #######mirror descent
    mirror_step <- 1 / (dual_step * A_max + max(abs(gammatr_Xtr_R)))
    mirror_res <- do.call('mirror_descent_cmp', c(list(fun = fun, grad = grad, wts = wts_curr, fun_arg = fun_arg, grad_arg = grad_arg, mirror_step = mirror_step), mirror_arg), quote = T)
    if (verbose) {
      if (!mirror_res$converged) {
        warning(paste('Mirror descent did not converge at iteration', j))
      }
    }
    wts_curr <- as.vector(mirror_res$wts)
    inner_nits[j] <- mirror_res$iter
    #######
    #######gamma update
    Xtr_R_wts <- Xtr_R %*% wts_curr
    gamma_new <- gamma_curr + dual_step * Xtr_R_wts
    gamma_sc <- max(abs(gamma_curr))
    gamma_sc <- ifelse(gamma_sc < .Machine$double.eps, .Machine$double.eps * 2, gamma_sc)
    gamma_diff <- max(abs(gamma_new - gamma_curr)) / gamma_sc
    gamma_curr <- gamma_new
    ##update computations and arguments
    gammatr_Xtr_R <- crossprod(gamma_curr, Xtr_R)
    R_X_gamma <- crossprod(Xtr_R, gamma_curr)
    fun_arg$gammatr_Xtr_R <- gammatr_Xtr_R
    grad_arg$R_X_gamma <- R_X_gamma
    ##
    #######
    #######function value
    # outer_fval_new <- do.call(fun, c(list(wts = wts_curr), fun_arg))
    outer_fval_new <- sum(log(wts_curr))
    outer_fval_sc <- abs(outer_fval_curr)
    outer_fval_sc <- ifelse(outer_fval_sc < .Machine$double.eps, .Machine$double.eps * 2, outer_fval_sc)
    outer_fval_diff <- abs(outer_fval_new - outer_fval_curr) / outer_fval_sc
    outer_fval_curr <- outer_fval_new
    outer_fval[j] <- outer_fval_curr
    #######
    #######primal residual
    primal_resid <- sqrt(sum((Xtr_R_wts)^2))
    dual_resid <- sqrt(sum((-1/wts_curr + R_X_gamma)^2))
    if (vary_penalty) {
      if (primal_resid > RB_mu * dual_resid) {
        dual_step <- RB_tau * dual_step
      } else if (dual_resid > RB_mu * primal_resid) {
        dual_step <- dual_step / RB_tau
      }
      ##update fun_arg and grad_arg
      fun_arg$dual_step <- dual_step
      grad_arg$dual_step <- dual_step
      ##
    }
    # print(primal_resid_new)
    print(sqrt(sum((-1/wts_curr + R_X_gamma)^2)))
    primal_resids[j] <- primal_resid
    #######
    #######stopping
    outer_tol <- switch(outer_tol_type,
                        gamma = gamma_diff,
                        fval = outer_fval_diff,
                        primal = primal_resid_curr,
                        primal_dual = primal_resid_curr)
    if (outer_tol_type == 'primal_dual') {
      if (primal_resid < primal_eps && dual_resid < dual_eps) {
        outer_converged <- T
        break
      }
    }
    else {
      if (outer_tol < outer_eps) {
        outer_converged = T
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
  outer_fval <- outer_fval[1:j]
  inner_nits <- inner_nits[1:j]
  logelr <- sum(log(wts_curr))
  ####
  
  return(list(logelr = logelr, logelr_stat = -2 * logelr, wts = wts_curr / n, sum_wts = sum(wts_curr) / n,
              gamma = gamma_curr, outer_converged = outer_converged, time = tot_time, 
              outer_fval = outer_fval, primal_resid = primal_resids, rho = dual_step,
              outer_tol = outer_tol, inner_nits = inner_nits, outer_nits = j))
}
##################

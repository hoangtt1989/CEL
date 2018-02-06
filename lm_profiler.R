library(pracma)
library(foreach)
library(doRNG)

##################bootstrap for coefficient
CEL_lm_beta_boot <- function(fix_index, X, y, B = 100, prob = .95, mean_init = NULL, algorithm = c('ADMM', 'nested'), parallel = T, seed = 123, ...) {
  ##checking inputs
  algorithm <- algorithm[1]
  if (!(algorithm %in% c('ADMM', 'nested'))) {
    stop('supply a valid algorithm - ADMM or nested')
  }
  ##
  ##initial parameters
  p <- ncol(X)
  if (is.null(mean_init)) {
    mean_init <- solve(crossprod(X), crossprod(X, y))
  }
  ##
  ##set up the reparam
  # T_reparam <- t(as.matrix(rep(0, p)))
  # T_reparam[fix_index] <- 1
  T_reparam <- diag(p)
  T_reparam <- T_reparam[fix_index, ]
  if (length(fix_index) == 1) {
    T_reparam <- t(as.matrix(T_reparam))
  }
  alpha_test <- mean_init[fix_index]
  s_hat_reparam <- rep(0, p)
  s_hat_reparam[fix_index] <- alpha_test
  # s_hat_reparam <- pinv(T_reparam) * alpha_test
  F_reparam <- orth(nullspace(T_reparam))
  ##
  set.seed(seed)
  if (parallel) {
    `%fun%` <- doRNG::`%dorng%`
  } else {
    `%fun%` <- `%do%`
  }
  start_time <- Sys.time()
  boot_res <- foreach::foreach(i = 1:B, .combine = cbind) %fun% {
    boot_beta_interm_fun(fix_index, alpha_test, X, y, F_reparam, s_hat_reparam, algorithm, ...)
  }
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  return(list(boot_stat = quantile(boot_res, prob), boot_vec = boot_res, time = as.double(tot_time, 'secs')))
}

boot_beta_interm_fun <- function(fix_index, alpha_test, X, y, F_reparam, s_hat_reparam, algorithm = c('ADMM', 'nested'), ...) {
  algorithm <- algorithm[1]
  ##resample data
  n <- nrow(X)
  sample_ind <- sample(1:n, n, replace = T)
  X <- X[sample_ind, ]
  y <- y[sample_ind]
  ##
  ##OLS estimate for initial value
  beta_OLS <- solve(crossprod(X), crossprod(X, y))
  beta_test <- beta_OLS
  beta_test[fix_index] <- alpha_test
  ##
  ##pass to optimization function
  res <- switch(algorithm,
                ADMM = CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, ...),
                nested = CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, ...))
  ##
  ##extract the test statistic
  return(res$logelr_stat)
  ##
}
##################

##################compute EL for a fixed prediction for linear regression using ADMM or nested algorithm
CEL_lm_pred_fixed <- function(alpha_add, fix_index, X, y, mean_init = NULL, algorithm = c('ADMM', 'nested'), ...) {
  ##checking inputs
  algorithm <- algorithm[1]
  if (!(algorithm %in% c('ADMM', 'nested'))) {
    stop('supply a valid algorithm - ADMM or nested')
  }
  ##
  ##initial parameters
  X_row <- t(as.matrix(X[fix_index, ]))
  if (is.null(mean_init)) {
    beta_init <- solve(crossprod(X), crossprod(X, y))
    pred_init <- as.numeric(X_row %*% beta_init)
  }
  ##
  ##set up the reparam
  T_reparam <- X_row
  alpha_test <- pred_init + alpha_add
  # s_hat_reparam <- as.matrix(rep(0, p))
  # s_hat_reparam[fix_index] <- alpha_test
  s_hat_reparam <- pinv(T_reparam) * alpha_test
  F_reparam <- orth(nullspace(T_reparam))
  ##
  ##set up the test
  beta_test <- beta_init
  # beta_test[fix_index] <- alpha_test
  ##
  ##run the chosen algorithm
  res <- switch(algorithm,
                ADMM = CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, ...),
                nested = CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, ...))
  ##
  return(c(list(mean_init = pred_init), res))
}
##################


##################compute EL for a fixed value of a component of the beta vector for linear regression
CEL_lm_beta_fixed <- function(beta_fix_val, fix_index, X, y, mean_init = NULL, algorithm = c('ADMM', 'nested'), ...) {
  ##checking inputs
  algorithm <- algorithm[1]
  if (!(algorithm %in% c('ADMM', 'nested'))) {
    stop('supply a valid algorithm - ADMM or nested')
  }
  ##
  ##initial parameters
  p <- ncol(X)
  if (is.null(mean_init)) {
    mean_init <- solve(crossprod(X), crossprod(X, y))
  }
  alpha_add <- beta_fix_val - mean_init[fix_index] 
  CEL_lm_beta_add(alpha_add, fix_index, X, y, mean_init, algorithm, ...)
  ##
}
##################

##################compute EL for a fixed value of a component of the beta vector (alpha_add specifies how far away from OLS) for linear regression
CEL_lm_beta_add <- function(alpha_add, fix_index, X, y, mean_init = NULL, algorithm = c('ADMM', 'nested'), ...) {
  ##checking inputs
  algorithm <- algorithm[1]
  if (!(algorithm %in% c('ADMM', 'nested'))) {
    stop('supply a valid algorithm - ADMM or nested')
  }
  ##
  ##initial parameters
  p <- ncol(X)
  if (is.null(mean_init)) {
    mean_init <- solve(crossprod(X), crossprod(X, y))
  }
  ##
  ##set up the reparam
  # T_reparam <- t(as.matrix(rep(0, p)))
  # T_reparam[fix_index] <- 1
  T_reparam <- diag(p)
  T_reparam <- T_reparam[fix_index, ]
  if (length(fix_index) == 1) {
    T_reparam <- t(as.matrix(T_reparam))
  }
  alpha_test <- mean_init[fix_index] + alpha_add
  s_hat_reparam <- rep(0, p)
  s_hat_reparam[fix_index] <- alpha_test
  # s_hat_reparam <- pinv(T_reparam) * alpha_test
  F_reparam <- orth(nullspace(T_reparam))
  ##
  ##set up the test
  beta_test <- mean_init
  beta_test[fix_index] <- alpha_test
  ##
  ##run the chosen algorithm
  res <- switch(algorithm,
                ADMM = CEL_lm_admm(beta_test, X, y, F_reparam, s_hat_reparam, ...),
                nested = CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, ...))
  ##
  return(c(list(mean_init = mean_init), res))
}
##################

##################find the EL confidence interval for a component of the beta vector for linear regression using ADMM or nested algorithm
CEL_lm_beta_profiler <- function(fix_index, X, y, conf_level = .05, test_thresh = 'chisq', upper_break = NULL, lower_break = 1e-6, upper_divide = 2, upper_increase = 1e-1, algorithm = c('ADMM', 'nested'), ...) {
  ##checking inputs
  algorithm <- algorithm[1]
  if (!(algorithm %in% c('ADMM', 'nested'))) {
    stop('supply a valid algorithm - ADMM or nested')
  }
  ##
  ##get critical value
  if (test_thresh == 'chisq') {
    test_thresh <- qchisq(1 - conf_level, df = 1)
  }
  ##
  ##get mean value
  beta_OLS <- solve(crossprod(X), crossprod(X, y))
  beta_fix <- beta_OLS[fix_index]
  ##
  ##make a function for root finding
  root_fun <- function(inp) {
    ret <- CEL_lm_beta_add(inp, fix_index, X, y, mean_init = beta_OLS, algorithm = algorithm, ...)
    ret$logelr_stat - test_thresh
  }
  ##
  #####make sure there is a sign change in the interval
  lower_init_fit <- root_fun(lower_break)
  if (lower_init_fit > 0) {
    stop('lower_break must produce a test statistic below the critical value - decrease it')
  }
  if (!is.null(upper_break)) {
    upper_init_fit <- root_fun(upper_break)
    if (upper_init_fit < 0) {
      stop('upper_break must produce a test statistic above the critical value - increase it')
    }
  } else { ## user did not supply upper_break - we will try to find one using the magnitude of beta
    upper_break <- abs(beta_fix) / upper_divide
    upper_init_fit <- root_fun(upper_break)
    while (upper_init_fit < 1e-2) {
      message('increasing upper_break to produce a test statistic above the critical value')
      upper_break <- upper_break + upper_increase
      upper_init_fit <- root_fun(upper_break)
    }
  } 
  #####
  ###use root finding to get the interval
  ##upper
  start_time <- Sys.time()
  upper_int <- uniroot(root_fun, lower = lower_break, upper = upper_break, check.conv = F)
  end_time <- Sys.time()
  upper_time <- end_time - start_time
  ##
  ##lower
  start_time <- Sys.time()
  lower_int <- uniroot(root_fun, lower = -upper_int$root, upper = -lower_break, check.conv = F, extendInt = 'downX')
  end_time <- Sys.time()
  lower_time <- end_time - start_time
  ##
  ###
  return(list(upper_int = beta_fix + upper_int$root, lower_int = beta_fix + lower_int$root, mean_val = beta_fix, time = as.double(upper_time + lower_time, unit = 'secs')))
}
##################

###UNFINISHED
##################find the EL confidence interval for a component of linear regression using ADMM or nested algorithm
CEL_lm_profiler <- function(fix_index, X, y, type = c('beta', 'prediction'), conf_level = .05, mean_val = NULL, upper_break = NULL, lower_break = 1e-6, upper_divide = 2, upper_increase = 1e-1, algorithm = c('ADMM', 'nested'), ...) {
  ##checking inputs
  type <- type[1]
  if (!(type %in% c('beta', 'prediction'))) {
    stop('supply a valid interval type - beta or prediction')
  }
  algorithm <- algorithm[1]
  if (!(algorithm %in% c('ADMM', 'nested'))) {
    stop('supply a valid algorithm - ADMM or nested')
  }
  ##
  ##get critical value
  test_thresh <- qchisq(1 - conf_level, df = 1)
  ##
  ##get mean value
  if (is.null(mean_val)) {
    beta_OLS <- solve(crossprod(X), crossprod(X, y))
    mean_val <- switch(type,
                       beta = beta_OLS[fix_index],
                       prediction = as.numeric(X[fix_index, ] %*% beta_OLS))
  }
  ##
  ##make a function for root finding
  root_fun <- function(inp) {
    ret <- switch(type,
                  beta = CEL_lm_beta_add(inp, fix_index, X, y, mean_init = beta_OLS, algorithm = algorithm, ...),
                  prediction = CEL_lm_pred_fixed())
    ret$logelr_stat - test_thresh
  }
  ##
  #####make sure there is a sign change in the interval
  lower_init_fit <- root_fun(lower_break)
  if (lower_init_fit > 0) {
    stop('lower_break must produce a test statistic below the critical value - decrease it')
  }
  if (!is.null(upper_break)) {
    upper_init_fit <- root_fun(upper_break)
    if (upper_init_fit < 0) {
      stop('upper_break must produce a test statistic above the critical value - increase it')
    }
  } else { ## user did not supply upper_break - we will try to find one using the magnitude of beta
    upper_break <- abs(beta_fix) / upper_divide
    upper_init_fit <- root_fun(upper_break)
    while (upper_init_fit < 1e-2) {
      message('increasing upper_break to produce a test statistic above the critical value')
      upper_break <- upper_break + upper_increase
      upper_init_fit <- root_fun(upper_break)
    }
  } 
  #####
  ###use root finding to get the interval
  ##upper
  start_time <- Sys.time()
  upper_int <- uniroot(root_fun, lower = lower_break, upper = upper_break, check.conv = T)
  end_time <- Sys.time()
  upper_time <- end_time - start_time
  ##
  ##lower
  start_time <- Sys.time()
  lower_int <- uniroot(root_fun, lower = -upper_break, upper = -lower_break, check.conv = T)
  end_time <- Sys.time()
  lower_time <- end_time - start_time
  ##
  ###
  return(list(upper_int = beta_fix + upper_int$root, lower_int = beta_fix + lower_int$root, mean_val = beta_fix, upper_time = upper_time, lower_time = lower_time))
}
##################

###UNFINISHED
##################get EL confidence interval for all relevant components for linear regression
CEL_lm_CI <- function(X, y, type = c('beta', 'prediction'), conf_level = .05, algorithm = c('ADMM', 'nested'), upper_break = NULL, lower_break = 1e-6, upper_divide = 2, upper_increase = 1e-1, ...) {
  ##checking inputs
  type <- type[1]
  if (!(type %in% c('beta', 'prediction'))) {
    stop('supply a valid interval type - beta or prediction')
  }
  algorithm <- algorithm[1]
  if (!(algorithm %in% c('ADMM', 'nested'))) {
    stop('supply a valid algorithm - ADMM or nested')
  }
  ##
  beta_OLS <- solve(crossprod(X), crossprod(X, y))
  mean_vals <- switch(type,
                      beta = beta_OLS,
                      prediction = X %*% beta_OLS)
  
}

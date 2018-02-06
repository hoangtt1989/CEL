grad_step <- function(x_curr, fun, grad, fun_arg = list(), grad_arg = list(),
                      ls_maxit = 100, ls_eps = 1e-12, ls_beta = .5, ls_alpha = .25) {
  fun_curr <- do.call('fun', c(list(x_curr), fun_arg), quote = T)
  grad_new <- do.call('grad', c(list(x_curr), grad_arg), quote = T)
  
  search_dir <- -grad_new
  x_new <- x_curr + search_dir
  gradtr_pk <- sum(as.numeric(grad_new) * as.numeric(search_dir))

  ls_fval <- do.call('fun', c(list(x_new), fun_arg), quote = T)
  step <- 1
  ls_converged <- ls_fval < (fun_curr + ls_alpha * step * gradtr_pk)
  if (!ls_converged) {
    for (ls_iter in 1:ls_maxit) {
      step <- ls_beta * step
      x_new <- x_curr + step * search_dir
      ls_fval <- do.call('fun', c(list(x_new), fun_arg), quote = T)
      ls_diff <- fun_curr + ls_alpha * step * gradtr_pk - ls_fval
      # print(ls_diff)
      if ((ls_diff >= 0 | abs(ls_diff) < ls_eps)) {
        ls_converged <- T
        break
      }
    }
  }
  return(list(x_curr = x_new, fval_new = ls_fval, fval_curr = fun_curr, gval_new = grad_new, ls_converged = ls_converged))
}

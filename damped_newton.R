##take a damped newton step - this is placed inside a for loop
damped_newton_step <- function(x_curr, fun, grad, hess, fun_arg = list(), grad_arg = list(), hess_arg = list(),
                          ls_maxit = 100, ls_eps = 1e-12, ls_beta = .7, ls_c1 = .25, mod_hess = F) {
  ##current function, gradient, hessian
  curr_fun <- do.call('fun', c(list(x_curr), fun_arg), quote = T)
  curr_grad <- do.call('grad', c(list(x_curr), grad_arg), quote = T)
  curr_hess <- do.call('hess', c(list(x_curr), hess_arg), quote = T)
  if (mod_hess) {
    curr_hess <- sfsmisc::posdefify(curr_hess) ###if the hessian is not positive definite - find the nearest one
  }
  ##
  #####update search direction
  pk <- -solve(curr_hess, curr_grad)
  gradtr_pk <- crossprod(curr_grad, pk)
  x_new <- x_curr + pk
  #####
  #########line search
  ls_fval <- do.call('fun', c(list(x_new), fun_arg), quote = T)
  step <- 1
  ls_converged <- ls_fval < (curr_fun + ls_c1 * step * gradtr_pk)
  if (!ls_converged) {
    for (ls_iter in 1:ls_maxit) {
      step <- ls_beta * step
      x_new <- x_curr + step * pk
      ls_fval <- do.call('fun', c(list(x_new), fun_arg), quote = T)
      ls_diff <- curr_fun + ls_c1 * step * gradtr_pk - ls_fval
      if (ls_diff >= 0 | abs(ls_diff) < ls_eps) {
        ls_converged <- T
        break
      }
    }
  }
  #########
  return(list(x_curr = x_new, fval_new = ls_fval, fval_curr = curr_fun, gval_new = curr_grad, ls_converged = ls_converged))
}
##

# ### modified cholesky decomposition
# mod_chol <- function(mat, delta = 1e-8, beta = NULL)
# {
#   nn=nrow(mat)
#   C=matrix(0,nrow=nn,ncol=nn)
#   L=C
#   d=rep(0,nn)
#   if (is.null(beta)) {
#     gamma <- max(diag(mat))
#     zi <- max(mat - diag(diag(mat)))
#     nu <- max(1, sqrt(nn^2 - 1))
#     beta <- max(gamma, zi/nu, delta)
#   }
#   #  ncol=ncol(mat)
#   for (j in 1:(nn-1))
#   {
#     
#     diag(C)=diag(mat)
#     maxi=which.max(diag(mat))
#     if (maxi!=j)
#     {
#       mat[,c(j,maxi)]=mat[,c(maxi,j)]
#       mat[c(j,maxi),]=mat[c(maxi,j),]
#     }
#     if (j==1&nn<=1)
#     {
#       C[(j+1):nn,j]=mat[(j+1):nn,j]
#     }
#     else
#     {
#       C[(j+1):nn,j]=mat[(j+1):nn,j]-C[(j+1):nn,1:(j-1)]%*%matrix(L[j,1:(j-1)],ncol=1)
#     }
#     theta=max(abs(C[(j+1):nn,j]))
#     d[j]=max(theta,(theta/beta)^2,delta)
#     L[(j+1):nn,j]=C[(j+1):nn,j]/d[j]
#     diag(C)[(j+1):nn]=diag(C)[(j+1):nn]-L[(j+1):nn,j]*C[(j+1):nn,j]
#   }
#   d[nn]=max(C[nn,nn],(C[nn,nn]/beta)^2,delta)
#   return(list(L=L+diag(1,nn),D=diag(d),C=C))
# }
# ###


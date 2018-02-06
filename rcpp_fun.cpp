#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <math.h>
using namespace Rcpp;
using namespace RcppEigen;

// [[Rcpp::depends(RcppEigen)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int isnan(double x) { return x != x; }
int isinf(double x) { return !isnan(x) && isnan(x - x); }

// [[Rcpp::export]]
double bregman_entropy_cpp(const VectorXd & x_new, const VectorXd & x_curr) {
  return x_new.dot(x_new.array().log().matrix()) - x_curr.dot(x_curr.array().log().matrix()) - (x_new - x_curr).dot((1.0 + x_curr.array().log()).matrix());
}

// [[Rcpp::export]]
double fun_lm_AL_cpp(const VectorXd & wts, const VectorXd & gamma, const MatrixXd & X, const VectorXd & R, const double & dual_step) {
  MatrixXd Xtr_R = (R.asDiagonal() * X).transpose();
  return(-wts.array().log().sum() + (Xtr_R * wts).dot(gamma) + .5 * dual_step * (Xtr_R * wts).array().pow(2).sum());
}
// double fun_lm_AL_cpp(const Map<VectorXd> & wts, const Map<VectorXd> & gamma, const Map<MatrixXd> & X, const Map<VectorXd> & R, const double & dual_step) {
//   MatrixXd Xtr_R = (R.asDiagonal() * X).transpose();
//   return(-wts.array().log().sum() + (Xtr_R * wts).dot(gamma) + .5 * dual_step * (Xtr_R * wts).array().pow(2).sum());
// }


// [[Rcpp::export]]
VectorXd grad_lm_AL_cpp(const VectorXd & wts, const VectorXd & gamma, const MatrixXd & X, const VectorXd & R, const double & dual_step) {
  MatrixXd Xtr_R = (R.asDiagonal() * X).transpose();
  return(-wts.array().inverse().matrix() + Xtr_R.adjoint() * gamma + dual_step * Xtr_R.adjoint() * (Xtr_R * wts));
}
// VectorXd grad_lm_AL_cpp(const Map<VectorXd> & wts, const Map<VectorXd> & gamma, const Map<MatrixXd> & X, const Map<VectorXd> & R, const double & dual_step) {
//   MatrixXd Xtr_R = (R.asDiagonal() * X).transpose();
//   return(-wts.array().inverse().matrix() + Xtr_R.adjoint() * gamma + dual_step * Xtr_R.adjoint() * (Xtr_R * wts));
// }


MatrixXd AtA(const MatrixXd & A) { // transpose of A times itself
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint());
}

// z_hat_new <- solve(crossprod(F_reparam, crossprod(Xtr_W_X) %*% F_reparam), crossprod(F_reparam, Xtr_W_X) %*% (crossprod(X, wts_new * y) + gamma_curr / dual_step - Xtr_W_X %*% s_hat_reparam))

// [[Rcpp::export]]
VectorXd closed_beta_lm_AL_cpp(const Map<VectorXd> & wts, const Map<VectorXd> & gamma, const Map<VectorXd> & y, const Map<MatrixXd> & X, const Map<MatrixXd> & F_reparam, const Map<VectorXd> & s_hat_reparam, const double & dual_step) {
  MatrixXd Xtr_W_X = X.adjoint() * wts.asDiagonal() * X;
  MatrixXd lhs = F_reparam.adjoint() * AtA(Xtr_W_X) * F_reparam;
  VectorXd rhs = F_reparam.adjoint() * Xtr_W_X * (X.adjoint() * (wts.asDiagonal() * y) + std::pow(dual_step, -1) * gamma - Xtr_W_X * s_hat_reparam);
  VectorXd soln = lhs.llt().solve(rhs);
  return soln;
}

// // [[Rcpp::export]]
// double tst_fun(const Map<VectorXd> & wts) {
//   return(wts.array().log().sum());
// }

// // [[Rcpp::export]]
// double tst_call(SEXP call_fval, const Map<VectorXd> & wts, SEXP env) {
//   ev = new Rcpp::EvalCompiled(call_fval, env);
// }

// [[Rcpp::export]]
List mirror_descent_lm_AL_cpp(const Map<VectorXd> & wts, const Map<VectorXd> & gamma, const Map<MatrixXd> & X, 
                              const Map<VectorXd> & R, const double & dual_step,
                              double mirror_step = 1.0, const int & maxit = 1000, const double & mirror_eps = 1e-8,
                              const std::string & tol_type = "fval", const bool & vary_step = false,
                              const bool & line_search = false, const double & ls_eps = 1e-10, const int & ls_maxit = 1000,
                              const double & ls_beta = .75) {
  //// initializing
  VectorXd x_curr = wts;
  VectorXd y_curr = wts;
  VectorXd nu_curr = wts;
  VectorXd nu_new = nu_curr;
  VectorXd x_new = x_curr;
  VectorXd y_new = y_curr;
  int n = wts.size();
  int p = gamma.size();
  VectorXd wts_curr = wts;
  VectorXd wts_new = wts_curr;
  bool converged = false;
  int i;
  double theta_curr;
  double theta_new;
  VectorXd grad_update(p);
  // line search
  double fval_y;
  double surr_curr;
  bool ls_converged = false;
  int iter_ls;
  double ls_diff;
  if (line_search) {
    mirror_step = 1.0;
  }
  //
  if (vary_step && !line_search) {
    mirror_step *= std::pow(2.0 * std::log(n), .5);
  }
  List ret;
  ////
  // initial values for tol
  double fval_new;
  double fval_curr;
  double fval_sc;
  double wts_sc;
  double logelr_new;
  double logelr_curr;
  double logelr_sc;
  double tiny = 1e-15;
  double tiny2 = std::pow(tiny, .25);
  if (!tol_type.compare("fval")) {
    fval_curr = fun_lm_AL_cpp(wts_curr, gamma, X, R, dual_step);
  } else if (!tol_type.compare("logelr")) {
    logelr_curr = wts_curr.array().log().sum();
  }
  double tol;
  //
  for (i = 1; i <= maxit; ++i) {
    theta_curr = 2.0 / (i + 1.0);
    theta_new = 2.0 / (i + 2.0);
    grad_update = grad_lm_AL_cpp(y_curr, gamma, X, R, dual_step);
    if (vary_step && !line_search) {
      mirror_step *= std::pow(i, -.5);
    }
    nu_new = (nu_curr.array() * exp((-mirror_step * grad_update / theta_curr).array())).matrix();
    nu_new = n * nu_new / nu_new.sum();
    x_new = (1 - theta_curr) * x_curr + theta_curr * nu_new;
    y_new = (1 - theta_new) * x_new + theta_new * nu_new;
    wts_new = x_new;
    //line search
    if (line_search) {
      fval_new = fun_lm_AL_cpp(wts_new, gamma, X, R, dual_step);
      fval_y = fun_lm_AL_cpp(y_curr, gamma, X, R, dual_step);
      surr_curr = fval_y + grad_update.dot(wts_new - y_curr) + std::pow(mirror_step, -1.0) * bregman_entropy_cpp(wts_new, y_curr);
      ls_converged = fval_new < surr_curr;
      if (!ls_converged || isinf(ls_converged) || isnan(ls_converged)) {
        for (iter_ls = 0; iter_ls < maxit; ++iter_ls) {
          mirror_step *= ls_beta;
          nu_new = (nu_curr.array() * exp((-mirror_step * grad_update / theta_curr).array())).matrix();
          nu_new = n * nu_new / nu_new.sum();
          x_new = (1 - theta_curr) * x_curr + theta_curr * nu_new;
          y_new = (1 - theta_new) * x_new + theta_new * nu_new;
          wts_new = x_new;
          fval_new = fun_lm_AL_cpp(wts_new, gamma, X, R, dual_step);
          surr_curr = fval_y + grad_update.dot(wts_new - y_curr) + std::pow(mirror_step, -1.0) * bregman_entropy_cpp(wts_new, y_curr);
          ls_diff = surr_curr - fval_new;
          if ((!isinf(ls_diff) && !isnan(ls_diff)) && (ls_diff >= 0 || (std::abs(ls_diff) < ls_eps))) {
            ls_converged = true;
            break;
          }
        }
      }
    }
    //
    //update current values
    x_curr = wts_new;
    y_curr = y_new;
    nu_curr = nu_new;
    //
    //check tolerance
    if (!tol_type.compare("fval")) {
      if (!line_search) {
        fval_new = fun_lm_AL_cpp(wts_new, gamma, X, R, dual_step);
      }
      fval_sc = std::abs(fval_curr);
      if (fval_sc < tiny) {
        fval_sc = tiny2;
      }
      tol = std::abs(fval_new - fval_curr) / fval_sc;
      fval_curr = fval_new;
    } else if (!tol_type.compare("wts")) {
      wts_sc = wts_curr.array().abs().maxCoeff();
      if (wts_sc < tiny) {
        wts_sc = tiny2;
      }
      tol = (wts_new - wts_curr).array().abs().maxCoeff() / wts_sc;
    } else if(!tol_type.compare("logelr")) {
      logelr_new = wts_new.array().log().sum();
      logelr_sc = std::abs(logelr_new);
      if (logelr_sc < tiny) {
        logelr_sc = tiny2;
      }
      tol = std::abs(logelr_new - logelr_curr) / logelr_sc;
      logelr_curr = logelr_new;
    }
    //
    // update weights, exit if tol
    wts_curr = wts_new;
    if (tol < mirror_eps) {
      converged = true;
      break;
    }
    //
  }
  ret["wts"] = wts_curr;
  ret["iter"] = i;
  ret["tol"] = tol;
  ret["converged"] = converged;
  ret["ls_converged"] = ls_converged;
  return ret;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//



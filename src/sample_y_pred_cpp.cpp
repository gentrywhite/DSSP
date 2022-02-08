#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::export(.sample_y_pred_cpp)]]
arma::mat sample_y_pred_cpp(const Rcpp::List& pars)
{
  int n = as<int>(pars["n"]);
  int N = as<int>(pars["N"]);
  int m = as<int>(pars["m"]);
  arma::colvec eta = Rcpp::as<arma::colvec>(pars["eta"]);
  arma::colvec delta = Rcpp::as<arma::colvec>(pars["delta"]);
  arma::colvec ev = Rcpp::as<arma::colvec>(pars["ev"]);
  arma::colvec Y = Rcpp::as<arma::colvec>(pars["Y"]);
  arma::mat v = Rcpp::as<arma::mat>(pars["v"]);
  arma::mat nu = Rcpp::as<arma::mat>(pars["nu"]);
  arma::mat VT = v.t();
  arma::mat samples = arma::randn(N, m);
  arma::mat S = arma::mat(ev.size(), ev.size(), arma::fill::ones);

  for(int i=0; i<N; ++i)
  {
    // declare all matrices outside the loop at the right dimensions (try one by one)
    // and then update - same as LAMBDA in sample_delta.cpp
    // Compute S
    // arma::mat S = v*arma::diagmat(1/(1+eta[i]*ev))*VT;
    S = v*arma::diagmat(1/(1+eta[i]*ev))*VT;
    arma::mat S_inv = v*arma::diagmat(1+eta[i]*ev)*VT;
    
    // Compute MU
    arma::mat MU = S*Y;

    // Compute Sigma
    arma::mat SIGMA = delta[i]*S ;

    // Partition MU and SIGMA
    arma::mat MU1 = MU.rows(0, n-1);
    arma::mat MU2 = MU.rows(n, n+m-1);

    arma::mat S11 = delta[i]*S_inv.submat(0,0,n-1,n-1);
    arma::mat S12 = SIGMA.submat(0,n,n-1,n+m-1);
    arma::mat S21 = SIGMA.submat(n,0,n+m-1,n-1);
    arma::mat S22 = SIGMA.submat(n,n,n+m-1,n+m-1);

    // Compute Residuals
    arma::mat RES = nu.submat(0, i, n-1, i) - MU1;

    // Compute mu_pred and sigma_pred for y_pred
    arma::mat MU_pred =  MU2+S21*S11*RES;
    arma::mat M_id = arma::mat(m,m,arma::fill::ones);
    arma::mat S_pred = delta[i]*M_id+S22+S21*S11*S12;

    samples.row(i) = MU_pred.t() + samples.row(i) * arma::chol(S_pred);
  }
  
  return samples.t();
}

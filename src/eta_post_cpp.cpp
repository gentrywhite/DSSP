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
/*
 * Need to incorporate the prior somehow, if possible to call as an external function,
 * or hard code as a Pareto prior when using Rcpp.
 *
 * Also need to include the IC bit as a parameter
 *
 * after sourcing this type: >ptr_name<-create_xptr("eta_post_cpp")
 *
 * then you can sample using
 *
 * rust::ru_rcpp(ptr_name,d = n = , etc.)
 *
 * experiments show that this is about 40-50 times faster.
 *
 */
//
// [[Rcpp::export]]
double eta_post_cpp(const Rcpp::NumericVector& x, const Rcpp::List& pars)
{
  // Read in the data
  arma::vec y = Rcpp::as<arma::vec>(x);
  double eta = arma::as_scalar(y);
  double nd = as<double>(pars["ND"]);
  arma::colvec ev = Rcpp::as<arma::colvec>(pars["EV"]);
  arma::colvec Q = Rcpp::as<arma::colvec>(pars["Q"]);
  //Compute the value
  arma::colvec arg1 = (1+eta*ev);
  arma::colvec larg1 = log(arg1);
  double a = -sum(larg1);
  arma::colvec lambda = 1.0/arg1;
  double alpha = nd*0.5;
  arma::colvec arg2 = Q%(1.0-lambda);
  double beta = 0.5*std::inner_product(Q.begin(),Q.end(),arg2.begin(),0.0);
  double ans = arma::as_scalar(0.5*a+alpha*log(eta)-(alpha)*log(beta)+lgamma(alpha));
  return ans;
}


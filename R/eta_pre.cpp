/*
 * Need to incorporate the prior somehow, if possible to call as an external function,
 * or hard code as a Pareto prior when using Rcpp.
 *
 * Also need to incldue the IC bit as a parameter
 *
 * after sourcing this type: >ptr_name<-create_xptr("preeta")
 *
 * then you can sample using
 *
 * rust::ru_rcpp(ptr_name,d = n = , etc.)
 *
 * experiments show that this is about 40-50 times faster.
 *
 */

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//[[Rcpp::export]]
double preeta(const Rcpp::NumericVector& x, const Rcpp::List& pars)
{
  arma::vec y = Rcpp::as<arma::vec>(x);
  double eta = arma::as_scalar(y);
  double nd = as<double>(pars["nd"]);
  arma::colvec ev = Rcpp::as<arma::colvec>(pars["ev"]);
  arma::colvec Q = Rcpp::as<arma::colvec>(pars["Q"]);
  arma::colvec arg1 = (1+eta*ev);
  arma::colvec larg1 = log(arg1);
  double lp = -2*log(eta+1);
  double a = -sum(larg1);
  arma::colvec lambda = 1.0/arg1;
  double alpha = nd*0.5;
  arma::colvec arg2 = Q%(1.0-lambda);
  double leta = log(eta);
  double beta = 0.5*std::inner_product(Q.begin(),Q.end(),arg2.begin(),0.0);
  double ans = arma::as_scalar(0.5*a+alpha*leta-(alpha)*log(beta)+lgamma(alpha))+lp;
  return ans;
}

// [[Rcpp::export]]
SEXP create_xptr(std::string fstr) {
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "preeta")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&preeta))) ;
  else
    return(Rcpp::XPtr<funcPtr>(R_NilValue)) ;
}

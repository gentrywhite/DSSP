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
// [[Rcpp::export(.sample_y_pred_cpp)]]
List sample_y_pred_cpp(const Rcpp::List& pars)
{
  // Read in the data
  // double N = as<double>(pars["N"]);
  // double eta = as<double>(pars["eta"]); // this will be a vector that's looped through later
  // int N = eta.size(); // for when we include the loop within this function...
  // arma::vec lambda_diag(N, arma::fill::zeros); // " " 
  
  // int i = 0;
  int i = as<int>(pars["i"]);
  int n = as<int>(pars["n"]);
  int m = as<int>(pars["m"]);
  arma::colvec eta = Rcpp::as<arma::colvec>(pars["eta"]);
  arma::colvec delta = Rcpp::as<arma::colvec>(pars["delta"]);
  arma::colvec EV = Rcpp::as<arma::colvec>(pars["EV"]);
  arma::colvec Y = Rcpp::as<arma::colvec>(pars["Y"]);
  arma::mat V = Rcpp::as<arma::mat>(pars["V"]);
  arma::mat nu = Rcpp::as<arma::mat>(pars["nu"]);
  
  arma::mat VT = V.t();
  arma::mat S = V*arma::diagmat(1/(1+eta[i]*EV))*VT;
  arma::mat S_inv = V*arma::diagmat(1+eta[i]*EV)*VT;
  arma::mat MU = S*Y;
  
  
  // ##  Compute Sigma
  arma::mat SIGMA = delta[i]*S ;
  
  // ##  Partition MU and SIGMA
  arma::mat MU1 = MU.rows(0, n-1);
  arma::mat MU2 = MU.rows(n, n+m-1);
  
  arma::mat S11 = delta[i]*S_inv.submat(0,0,n-1,n-1);
  arma::mat S12 = SIGMA.submat(0,n,n-1,n+m-1);
  arma::mat S21 = SIGMA.submat(n,0,n+m-1,n-1);
  arma::mat S22 = SIGMA.submat(n,n,n+m-1,n+m-1);

  // ##  Compute Residuals
  arma::mat RES = nu.submat(0, i, n-1, i) - MU1;

  // ##  Compute mu.pred and sigma.pred for y.pred
  arma::mat MU_pred =  MU2+S21*S11*RES;
  arma::mat M_id = arma::mat(m,m,arma::fill::ones);
  arma::mat S_pred = delta[i]*M_id+S22+S21*S11*S12;
//       MU.pred<-MU2+S21%*%S11%*%RES
//         S.pred<-delta[i]*diag(1,m)+S22+S21%*%S11%*%S12
  arma::mat SAMPLES = trans(S_pred*MU_pred);
  
  
  List ret;
  ret["S"] = S;
  ret["S_inv"] = S_inv;
  ret["MU"] = MU;
  ret["SIGMA"] = SIGMA;
  ret["MU1"] = MU1;
  ret["MU2"] = MU2;
  ret["S11"] = S11;
  ret["S12"] = S12;
  ret["S21"] = S21;
  ret["S22"] = S22;
  ret["RES"] = RES;
  ret["MU_pred"] = MU_pred;
  ret["S_pred"] = S_pred;
  ret["SAMPLES"] = SAMPLES;
  return ret;

}

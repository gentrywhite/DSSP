#include<math.h>
#include<cmath>
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
// [[Rcpp::export]]
Rcpp::List make_M(const Rcpp::NumericMatrix& x)
{
  // Read in the data
    arma::mat X = Rcpp::as<arma::mat>(x);
  //Compute the value
    int N = X.n_rows;
//    int M = X.n_cols; 
    arma::colvec Y(N,arma::fill::ones);
    arma::mat T =join_rows(Y,X);
/*    d<-ncol(Tmat)
    D<-as.matrix(dist(X))
    ind0<-D!=0
  K<-D
    K[ind0]<-tps.rbf(D[ind0],even)
    Kstar<-D
    Kstar[ind0]<-D[ind0]^2*log(D[ind0])
    TT<-tcrossprod(Tmat,Tmat)
    F.mat<-eigen(TT,symmetric = TRUE)
    F2<-F.mat$vectors[,-c(1:d)]
  KF2<-crossprod(K,F2)
    G<-cbind(Tmat,KF2)
    H<-matrix(0,n,n)
    H[-c(1:d),-c(1:d)]<-crossprod(KF2,F2)
    G.inv<-qr.solve(G)
    HG<-crossprod(H,G.inv)
    M<-crossprod(HG,G.inv)
    M<-as.matrix(Matrix::forceSymmetric(M))
    M.eigen<-eigen(M,symmetric = TRUE)
    return(list(M=M,M.eigen=M.eigen))*/
  return Rcpp::List::create(Rcpp::Named("Tmat")=T);
}


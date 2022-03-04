#' DSSP R Functions
#' @title DSSP
#' @name DSSP
#' @description does some DSSP stuff
## usethis namespace: start
#' @useDynLib DSSPcpp, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
#' TPS radial basis function
#'
#' Function to compute the thin-plate splines radial basis function for internal use by the function make.M().
#' @param x is a Euclidean distance between two points.
#' @param is.even is a logical argument indicating TRUE if the dimension of the space where the thin-plate spline smoother is being fitted is even.
#' @keywords thin-plate spline basis function
#' @return The resulting value of the thin-plate spline radial basis function.
#' @details This function computes the thin-plate spline radial basis function depending on the if d is odd or even.
#' @export
#' @examples
#' ## Use the Meuse River dataset from the package 'gstat'
#'
#' library(sp)
#' library(gstat)
#' data(meuse.all)
#' coordinates(meuse.all) <- ~ x + y
#' X <- scale(coordinates(meuse.all))
#' D <- as.matrix(dist(X))
#' K <- tps.rbf(D, TRUE)
tps.rbf <- function(x, is.even) {
  if (is.even == FALSE) {
    x^2
  } else {
    log(x) * x^2
  }
}

#' Precision Matrix Function
#'
#' This function creates the precision matrix for the spatial prior based on thin-plate splines 
#' and returns the matrix M, and its eigenvalues and eigenvectors
#' @param X a matrix of spatial coordinates. It is recommended that the coordinates be scaled and centred.
#' @keywords spatial prior, thin-plate splines
#' @return A list containing the precision matrix M and the object M.eigen containing 
#' eigenvalues and eigenvectors for the matrix M.
#' @details The M matrix is the precision matrix for the
#'  spatial effects from the direct sampling spatial prior (DSSP) model. M is based on
#'  thin plate splines basis functions, see White et. al. 2019 for more details on how the
#'  matrix M is constructed.
#' @export
#' @examples
#' ## Use the Meuse River dataset from the package 'gstat'
#'
#' library(sp)
#' library(gstat)
#' data(meuse.all)
#' coordinates(meuse.all) <- ~ x + y
#' X <- scale(coordinates(meuse.all))
#' make.M(X)
make.M <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  dimX <- ncol(X)
  even <- dimX %% 2 == 0
  deg <- trunc(dimX / 2 + 1) - 1
  Tmat <- cbind(1, stats::poly(X, degree = deg, raw = TRUE))
  d <- ncol(Tmat)
  D <- as.matrix(stats::dist(X))
  ind0 <- D != 0
  K <- D
  K[ind0] <- tps.rbf(D[ind0], even)
  Kstar <- D
  Kstar[ind0] <- D[ind0]^2 * log(D[ind0])
  TT <- tcrossprod(Tmat, Tmat)
  F.mat <- eigen(TT, symmetric = TRUE)
  F2 <- F.mat$vectors[, -c(1:d)]
  KF2 <- crossprod(K, F2)
  G <- cbind(Tmat, KF2)
  H <- matrix(0, n, n)
  H[-c(1:d), -c(1:d)] <- crossprod(KF2, F2)
  G.inv <- qr.solve(G)
  HG <- crossprod(H, G.inv)
  M <- crossprod(HG, G.inv)
  M <- as.matrix(Matrix::forceSymmetric(M))
  M.eigen <- eigen(M, symmetric = TRUE)
  list(M = M, M.eigen = M.eigen)
}

#' Data Vector Function
#'
#' This function creates the the data vector for evaluating the joint posterior distribution of 
#' the direct sampling spatial prior (DSSP) model.
#' @param y a vector of observed data.
#' @param V a matrix of eigen vectors from the precision matrix computed by \code{make.M()}.
#' @keywords spatial prior, thin-plate splines
#' @return A vector Q = y' V which is used in subsequent functions for sampling
#'  from the  joint posterior distribution.
#' @export
#' @examples
#' ## Use the Meuse River dataset from the package 'gstat'
#'
#' library(sp)
#' library(gstat)
#' data(meuse.all)
#' coordinates(meuse.all) <- ~ x + y
#' X <- scale(coordinates(meuse.all))
#' tmp <- make.M(X)
#'
#' EV <- tmp$M.eigen$values
#' V <- tmp$M.eigen$vectors
#'
#' Y <- scale(log(meuse.all$zinc))
#' Q <- make.Q(Y, V)
make.Q <- function(y, V) {
  crossprod(y, V)
}

#' Function to sample from the posterior of the smoothing parameter eta conditioned on the data y.
#'
#' This function samples from the log-posterior density of the smoothing parameter from the
#'   thin-plate splines based spatial prior using a ratio-of-uniform sampler.
#' @param N the number of samples desired.
#' @param ND the rank of the precision matrix, the default value is n-3 for spatial data.
#' @param EV eigenvalues of the precision matrix spatial prior from the function make.M().
#' @param Q the data vector from the make.Q function.
#' @param UL the upper limit for the smoothing parameter value; used for the 
#' ratio-of-uniform sampler, default is 1000.
#' @param log_prior a function of x evaluating the log of the prior density for eta
#' @keywords spatial prior, thin-plate splines
#' @return N samples drawn from the posterior of eta given the data y \eqn{\pi(eta | y)}.
#' @export
#' @examples
#' ## Use the Meuse River dataset from the package 'gstat'
#'
#' library(sp)
#' library(gstat)
#' data(meuse.all)
#' coordinates(meuse.all) <- ~ x + y
#' X <- scale(coordinates(meuse.all))
#' tmp <- make.M(X)
#'
#' EV <- tmp$M.eigen$values
#' V <- tmp$M.eigen$vectors
#'
#' M <- tmp$M
#'
#' Y <- scale(log(meuse.all$zinc))
#' Q <- make.Q(Y, V)
#'
#' ND <- nrow(X) - 3
#' f <- function(x) -x ## log-prior for exponential distribution for the smoothing parameter
#'
#' ## Draw 100 samples from the posterior of eta given the data y.
#' sample.eta(100, ND, EV, Q, UL = 1000, f)
sample.eta <- function(N, ND, EV, Q, UL = 1000, log_prior) {
  RES <- rust::ru(function(x) .eta_post_cpp(x, list(ND = ND, EV = EV, Q = Q)) + log_prior(x), n = N, d = 1, init = 1, trans = "BC", upper = UL)
  RES$sim_vals
}


#' Function to sample from the posterior of the variance parameter
#'
#' This function samples from the log-posterior density of the variance parameter from the likelihood
#' @param eta samples of the smoothing parameter from the sample.eta function.
#' @param ND the rank of the precision matrix, the default value is n-3 for spatial data.
#' @param EV eigenvalues of the precision matrix spatial prior from the function make.M().
#' @param Q the data vector from the make.Q function.
#' @param pars a vector of the prior shape and rate parameters for the 
#' inverse-gamma prior distribution of delta.
#' @keywords spatial prior, thin-plate splines
#' @return N samples drawn from the posterior of \eqn{\pi(delta | eta, y)}.
#' @export
#' @examples
#' ## Use the Meuse River dataset from the package 'gstat'
#'
#' library(sp)
#' library(gstat)
#' data(meuse.all)
#' coordinates(meuse.all) <- ~ x + y
#' X <- scale(coordinates(meuse.all))
#' tmp <- make.M(X)
#'
#' M <- tmp$M
#'
#' Y <- scale(log(meuse.all$zinc))
#'
#' ND <- nrow(X) - 3
#' M.list <- make.M(X) ##  Only Needs to return the eigenvalues and vectors
#' M <- M.list$M
#' EV <- M.list$M.eigen$values
#' V <- M.list$M.eigen$vectors
#' Q <- make.Q(Y, V)
#'
#' f <- function(x) -x ## log-prior for exponential distribution for the smoothing parameter
#' ## Draw 100 samples from the posterior of eta given the data y.
#'
#' ETA <- sample.eta(100, ND, EV, Q, f, UL = 1000)
#' DELTA <- sample.delta(ETA, ND, EV, Q, pars = c(0.001, 0.001))
#' ##  Old Slow Version of sample.nu()
#' ## sample.delta<-function(eta,nd,ev,Q,pars)
#' ## {
#' ##   N<-length(eta)
#' ##   f.beta<-function(x)
#' ##   {
#' ##     lambda<-1/(1+x*ev)
#' ##     b<-tcrossprod(Q,diag(1-lambda))
#' ##     beta<-0.5*tcrossprod(Q,b)+pars[2]
#' ##     return(beta)
#' ##   }
#' ##   alpha<-pars[1]+nd*0.5
#' ##   beta<-sapply(eta,f.beta)
#' ##   delta<-1/rgamma(N,shape=alpha,rate=beta)
#' ##   return(delta)
#' ## }
sample.delta <- function(eta, ND, EV, Q, pars) {
  .sample_delta_cpp(eta, list(ND = ND, EV = EV, Q = Q, PARS = pars))
}


#' Function to sample from the posterior of the spatial effects
#'
#' This function samples from the posterior density of the spatial effects from the direct sampling
#'  spatial prior (DSSP) model.
#' @param Y vector of observed data.
#' @param eta samples of the smoothing parameter from the \code{sample.eta} function.
#' @param delta samples of the variance parameter from the \code{sample.delta} function.
#' @param EV eigenvalues of the precision matrix spatial prior from the function \code{make.M()}.
#' @param V eigenvectors of the precision matrix spatial prior from the function \code{make.M()}.
#' @keywords spatial prior, thin-plate splines
#' @keywords spatial prior, thin-plate splines
#' @return A matrix of samples with each column a random draw from the posterior 
#' of the spatial effects from the DSSP model \eqn{\pi(nu | eta, delta, y)}.
#' @export
#' @examples
#' ## Use the Meuse River dataset from the package 'gstat'
#'
#' library(sp)
#' library(gstat)
#' data(meuse.all)
#' coordinates(meuse.all) <- ~ x + y
#' X <- scale(coordinates(meuse.all))
#' tmp <- make.M(X)
#'
#' EV <- tmp$M.eigen$values
#' V <- tmp$M.eigen$vectors
#'
#' Y <- scale(log(meuse.all$zinc))
#' Q <- make.Q(Y, V)
#'
#' ND <- nrow(X) - 3
#' f <- function(x) -x ## log-prior for exponential distribution for the smoothing parameter
#' ## Draw 100 samples from the posterior of eta given the data y.
#'
#' ETA <- sample.eta(100, ND, EV, Q, f, UL = 1000)
#' DELTA <- sample.delta(ETA, ND, EV, Q, pars = c(0.001, 0.001))
#' NU <- sample.nu(Y, ETA, DELTA, EV, V)
sample.nu <- function(Y, eta, delta, EV, V) {
  .sample_nu_cpp(Y, list(eta = eta, delta = delta, EV = EV, V = V))
}

##  Wrapper function takes X,y,num_samples and prior for eta and returns samples from joint posterior

##  DSSP (Direct Sampling Spatial Prior)

#' Wrapper function to draw samples from the Direct Sampling Spatial Prior (DSSP) model
#'
#' This function samples from the log-posterior of all parameters in the model and returns a list
#'    object containing the samples. It performs a few compatibility checks on the inputs, then
#'    calls the sample.eta(), sample.delta(), and sample.nu().
#' @param formula a two sided linear formula with the response on left and the covariates on the right.
#' @param data a \code{data.frame} or \code{sp::SpatialPointsDataFrame} containing the response variable, covariates and coordinates.
#' @param N is the number of random samples to be drawn from the joint posterior for eta, delta, and nu.
#' @param pars a vector of the prior shape and rate parameters for the inverse-gamma 
#' prior distribution of delta, the variance parameter for the Gaussian likelihood.
#' @param log_prior a function evaluating the log of the prior density of eta. Default to be \code{function(x) -x}.
#' @param fitted.values return a matrix containing samples of the fitted values at each location X, defaults to \code{FALSE}.
#' @param coords spatial coordinates passed as the \code{value} argument to \code{sp::coordinates()}.
#' @keywords spatial prior, thin-plate splines
#' @return A list containing N samples of nu, eta, delta, and the original data X and Y.
#' @details
#'  The direct sampling spatial prior model assumes that the spatial model can be written
#'  as the likelihood parameterised with mean vector nu and variance delta
#'  \deqn{(y | nu, delta) ~ N(nu, delta * I)}
#'  where I is the identity matrix.  The prior for the vector of spatial effects nu is
#'  improper but is proportional to
#'  \deqn{\pi(nu | eta) propto (det(M)/2\pi)^{1/2} * exp(-eta nu' M nu/2),}
#'  the prior for delta is assumed to be a inverse-gamma distribution
#'  \deqn{(delta) ~ IG(a,b)}
#'  and the prior for eta can be specified for the user as any valid density function for eta > 0.
#' @export
#' @examples
#' ## Use the Meuse River dataset from the package 'gstat'
#'
#' library(sp)
#' library(gstat)
#' data(meuse.all)
#' coordinates(meuse.all) <- ~ x + y
#'
#' f <- function(x) -x ## log-prior for exponential distribution for the smoothing parameter
#'
#' ## Draw 100 samples from the posterior of eta given the data y.
#' OUTPUT <- DSSP(
#'   formula = log(zinc) ~ 1, data = meuse.all, N = 100,
#'   pars = c(0.001, 0.001), log_prior = f, fitted.values = FALSE
#' )
DSSP <- function(formula, data, N, pars, log_prior=function(x) -x, fitted.values = FALSE, coords = NULL) {
  stopifnot(is.function(log_prior))
  
  if (all(class(data) != "SpatialPointsDataFrame")) {
    sp::coordinates(data) <- coords
  } else {
    if (!is.null(coords)) message("obtaining spatial coordinates from data; ignoring coords provided")
    coords <- NULL
  }

  w <- sp::coordinates(data)
  if (any(!grepl("scaled", names(attributes(w))))) {
    w <- scale(w)
    coord_scaling <- list(
      center = attr(w, "scaled:center"),
      scale = attr(w, "scaled:scale")
    )
  } else {
    coord_scaling <- list(center = NA, scale = NA)
  }

  mt <- stats::terms(formula, data = data)
  mf <- stats::lm(formula, data = data, method = "model.frame")
  y <- stats::model.extract(mf, "response")
  y <- scale(y)
  y_scaling <- list(
    center = attr(y, "scaled:center"),
    scale = attr(y, "scaled:scale")
  )
  x <- stats::model.matrix(mt, mf)

  if (attr(x, "assign") == 0 & dim(x)[2] == 1) {
    # if model is y ~ 1: x is just the coordinates and old version of DSSP can be run
    x <- w
    intercept_only <- TRUE
  } else {
    # use both covariates and coordinates for model fitting
    intercept_only <- FALSE
    return(message("Currently only works for intercept only models: formula=y~1"))
  }

  N <- as.integer(N)
  X <- as.matrix(x)
  Y <- as.numeric(y)
  n <- length(y)
  pars <- as.numeric(pars)
  delta <- numeric(length = N)
  eta <- numeric(length = N)

  ##  Declare dimensional constants
  n <- length(Y)
  d <- ncol(X) + 1
  nd <- n - d
  ND <- n - d

  ##  Compute M
  M.list <- make.M(X) ##  Only Needs to return the eigenvalues and vectors
  M <- M.list$M

  EV <- M.list$M.eigen$values
  V <- M.list$M.eigen$vectors
  Q <- make.Q(Y, V)

  ## sample eta
  eta <- sample.eta(N, ND, EV, Q, log_prior, UL = 1000)

  ##  sample delta_0
  delta <- sample.delta(eta, ND, EV, Q, pars)

  out <- list(eta = eta, delta = delta)
  if (fitted.values) {
    nu <- sample.nu(Y, eta, delta, EV, V)
    out <- append(out, list(nu = nu, y_fitted = nu * y_scaling$scale + y_scaling$center))
  }

  append(out, list(
    X = X, Y = Y, y_scaling = y_scaling, coord_scaling = coord_scaling,
    coords = coords, intercept_only = intercept_only
  ))
}

# .onUnload <- function (libpath) {
#   library.dynam.unload("DSSPcpp", libpath)
# }

##  Predict method for DSSP

#' Function to generate samples from the posterior predictive distribution of y for new locations x.new
#'
#'
#' This function samples from the log-posterior of all parameters in the model
#' @param dssp.model the model output from \code{DSSP()}.
#' @param x.pred a matrix of spatial coordinates to predict new values of y.
#' @keywords spatial prior, thin-plate splines
#' @return A matrix with N columns containing N values from the posterior predictive distribution.
#' @export
#' @examples
#' ## Use the Meuse River dataset from the package 'gstat'
#'
#' library(sp)
#' library(gstat)
#' data(meuse.all)
#' coordinates(meuse.all) <- ~ x + y
#'
#' f <- function(x) -x ## log-prior for exponential distribution for the smoothing parameter
#'
#' ## Draw 100 samples from the posterior of eta given the data y.
#' OUTPUT <- DSSP(
#'   formula = log(zinc) ~ 1, data = meuse.all[1:155, ], N = 100,
#'   pars = c(0.001, 0.001), log_prior = f
#' )
#' Y.PRED <- DSSP.predict(OUTPUT, meuse.all[156:164, ])
DSSP.predict <- function(dssp.model, x.pred) { ##  function to generate samples
  ##  Needs the conditional predictive for y_pred and the output of eta, delta, and nu
  ##  Needs to be able to sample from posterior predicitive of nu_pred as well
  ##  First draw nu_pred then use that to draw y_pred from likelihood

  ##  Fix this so that  if nu is missing it is
  ##  computed first then the predict scheme is run.

  #  Extract components from dssp.model

  if (!any(class(x.pred) %in% c("SpatialPointsDataFrame", "SpatialPoints"))) {
    sp::coordinates(x.pred) <- dssp.model$coords
  }
  w.pred <- sp::coordinates(x.pred)
  if (any(!grepl("scaled", names(attributes(w.pred))))) {
    w.pred <- scale(w.pred, center = dssp.model$coord_scaling$center, scale = dssp.model$coord_scaling$scale)
  }
  if (dssp.model$intercept_only) {
    x.pred <- w.pred
  }

  x <- dssp.model$X
  y <- dssp.model$Y
  eta <- dssp.model$eta
  delta <- dssp.model$delta

  if ("nu" %in% names(dssp.model)) {
    nu <- dssp.model$nu
  } else {
    m <- make.M(x)
    nu <- sample.nu(y, eta, delta, m$M.eigen$values, m$M.eigen$vectors)
  }
  N <- length(eta)
  n <- length(y)
  X <- rbind(x, x.pred)
  m <- nrow(X) - n
  Y <- c(y, rep(0, m))

  ##  Make augmented M matrix
  M.list <- make.M(X)

  ##  Extract Vectors
  v <- M.list$M.eigen$vectors

  ##  Extract Values
  ev <- M.list$M.eigen$values

  y.pred <- .sample_y_pred_cpp(list(N = N, eta = eta, ev = ev, v = v, Y = Y, delta = delta, n = n, m = m, nu = nu))

  y.pred * dssp.model$y_scaling$scale + dssp.model$y_scaling$center
}

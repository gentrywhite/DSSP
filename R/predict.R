predict.dsspMod <- function(object, newdata, ...) {
  if (missing(newdata)) {
    if("y_fitted" %in% names(object)) return(object$y_fitted)
    m <- make.M(object$X, intercept_only = object$intercept_only)
    nu <- sample.nu(object$Y, object$eta, object$delta, m$M.eigen$values, m$M.eigen$vectors)
    return(nu * object$y_scaling$scale + object$y_scaling$center)
  } 
  if (!any(class(newdata) %in% c("SpatialPointsDataFrame", "SpatialPoints"))) {
    sp::coordinates(newdata) <- object$coords
  }
  w.pred <- sp::coordinates(newdata)
  if (any(!grepl("scaled", names(attributes(w.pred))))) {
    w.pred <- scale(w.pred, center = object$coord_scaling$center, scale = object$coord_scaling$scale)
  }
  if (object$intercept_only) {
    newdata <- w.pred
  }
  
  x <- object$X
  y <- object$Y
  eta <- object$eta
  delta <- object$delta
  
  if ("nu" %in% names(object)) {
    nu <- object$nu
  } else {
    m <- make.M(x, intercept_only = object$intercept_only)
    nu <- sample.nu(y, eta, delta, m$M.eigen$values, m$M.eigen$vectors)
  }
  N <- length(eta)
  n <- length(y)
  X <- rbind(x, newdata)
  m <- nrow(X) - n
  Y <- c(y, rep(0, m))
  
  ##  Make augmented M matrix
  M.list <- make.M(X, intercept_only = object$intercept_only)
  
  ##  Extract Vectors
  v <- M.list$M.eigen$vectors
  
  ##  Extract Values
  ev <- M.list$M.eigen$values
  
  y.pred <- .sample_y_pred_cpp(list(N = N, eta = eta, ev = ev, v = v, Y = Y, delta = delta, n = n, m = m, nu = nu))
  
  y.pred * object$y_scaling$scale + object$y_scaling$center
}

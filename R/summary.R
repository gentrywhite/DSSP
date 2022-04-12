summary.dsspMod <- function(object, prob = 0.95, robust = FALSE, mc_se = FALSE, ...) {
  out <- list(
    formula = object$formula,
    nobs = object$nobs,
    niter = object$N
  )
  probs <- validate_ci_bounds(prob)
  
  variables <- list(eta=object$eta, delta=object$delta)
  if(!is.null(object$nu)) {
    variables$nu = object$nu
  }
  
  measures <- list()
  if (robust) {
    measures$Estimate <- stats::median
    if (mc_se) {
      measures$MCSE <- posterior::mcse_median
    }
    measures$Est.Error <- stats::mad
  } else {
    measures$Estimate <- mean
    if (mc_se) {
      measures$MCSE <- posterior::mcse_mean
    }
    measures$Est.Error <- stats::sd
  }
  full_summary_measures <- c(measures, list(
    ll = function(x) stats::quantile(x, probs=probs[1]),
    ul = function(x) stats:: quantile(x, probs=probs[2]),
    Rhat = posterior::rhat,
    ESS = ess
  ))
  
  full_summary <- lapply(full_summary_measures, function(m) sapply(variables, function(v) m(v)))
  full_summary <- as.data.frame(full_summary)
  
  prob <- probs[2] - probs[1]
  
  names(full_summary)[which(names(full_summary) %in% c("ll", "ul"))] <- paste0(c("l-", "u-"), prob * 100, "% CI")
  rownames(full_summary) <- names(variables)
  
  out <- c(out, list(full_summary=full_summary))
  if(!is.null(object$covariates_posterior)){
    cov_measures <- c(measures, list(
      min = min,
      `q0.025` = function(x) stats::quantile(x, probs=0.025),
      `q0.25` = function(x) stats::quantile(x, probs=0.25),
      `q0.50` = function(x) stats::quantile(x, probs=0.50),
      `q0.75` = function(x) stats::quantile(x, probs=0.75),
      `q0.975` = function(x) stats::quantile(x, probs=0.975),
      max = max
    ))
    
    cov_list <- stats::setNames(
      split(object$covariates_posterior, seq(nrow(object$covariates_posterior))),
      rownames(object$covariates_posterior)
    )
    
    cov_summary <- lapply(cov_measures, function(m) sapply(cov_list, function(v) m(v)))
    cov_summary <- as.data.frame(cov_summary)
    out <- c(out, list(cov_summary=cov_summary))
  }
  class(out) <- "dsspModsummary"
  out
}

print.dsspModsummary <- function(x, digits = 2, ...) {
  cat("Formula: ")
  print(x$formula)
  cat(paste0("Number of observations: ", x$nobs, "\n"))
  cat(paste0("Number of iterations: ", x$niter, "\n\n"))
  cat("Summary of model:\n")
  print_format(x$full_summary)
  if(!is.null(x$cov_summary)){
    cat("\nSummary of covariates:\n")
    print_format(x$cov_summary)
  }
}

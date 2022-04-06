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
  measures <- c(measures, list(
    ll = function(x) stats::quantile(x, probs=probs[1]),
    ul = function(x) stats:: quantile(x, probs=probs[2]),
    Rhat = posterior::rhat,
    Bulk_ESS = posterior::ess_bulk,
    Tail_ESS = posterior::ess_tail
  ))
  
  full_summary <- lapply(measures, function(m) sapply(variables, function(v) m(v)))
  full_summary <- as.data.frame(full_summary)
  
  prob <- probs[2] - probs[1]
  
  names(full_summary)[which(names(full_summary) %in% c("ll", "ul"))] <- paste0(c("l-", "u-"), prob * 100, "% CI")
  rownames(full_summary) <- names(variables)
  
  out <- c(out, list(full_summary=full_summary))
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
}
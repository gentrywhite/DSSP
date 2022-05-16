#' Diagnostic, Density and Contour Plots
#'
#' @param x an object of class \code{dsspMod}
#' @param robust_residuals whether to use robust residuals (median of predicted).
#'   Default to be \code{TRUE}.
#' @param contour_plots whether or not to return a second panel with contour plots.
#'   Defaults to \code{TRUE}
#' @param nx dimension of output grid in x direction.
#'   Used for interpolation (\code{akime::interp()}). 
#' @param ny dimension of output grid in y direction.
#'   Used for interpolation (\code{akime::interp()}). 
#' @param pal colour palette used for filled contour plot.
#' @param nlevels number of levels used in contour plots.
#' @param ... additional arguments that are passed to \code{ggplot2::scale_fill_distiller()}.
#'
#' @return NULL
#' @export 
#'
#' @examples
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
#'   pars = c(0.001, 0.001), log_prior = f
#' )
#' plot(OUTPUT)
plot.dsspMod <- function(x,
                         robust_residuals=TRUE,
                         contour_plots=TRUE,
                         nx=100, ny=100, nlevels=5,
                         ...){
  if (!inherits(x, "dsspMod"))	stop("use only with \"dsspMod\" objects")
  r <- residuals.dsspMod(x, robust=robust_residuals)
  ylim <- range(r)
  ypred <- predict.dsspMod(x)
  metric <- ifelse(robust_residuals, stats::median, mean)
  yh <- apply(ypred, 1, metric)
  y <- reverse_scaling(x$Y, x$y_scaling)
  
  df_params <- data.frame(eta=x$eta, delta=x$delta)
  df_params$n <- seq.int(nrow(df_params))
  
  resid_vs_fitted <- 
    ggplot() +
    geom_point(aes(x=yh, y=r)) +
    geom_smooth(aes(x=yh, y=r), se=FALSE, col="black") +
    labs(
      title="Residuals vs fitted values",
      x=paste("fitted", x$dep_var),
      y="residuals"
    )

  predicted_vs_actual <-
    ggplot() +
    geom_point(aes(y, yh)) +
    geom_abline(slope=1, intercept=0) +
    labs(
      title="Observe vs fitted values",
      x=paste("observed", x$dep_var),
      y=paste("fitted", x$dep_var)
    )

  eta_density <- 
    ggplot() +
    geom_density(aes(x$eta)) +
    labs(title=expression("Posterior Density of " * eta), x=expression(eta))
  
  eta_trace <-
    ggplot(data=df_params) +
    geom_line(aes(x=n, y=eta)) +
    labs(title=expression("Traceplot of " * eta), y=expression(eta))
  
  delta_density <-
    ggplot() +
    geom_density(aes(x$delta)) +
    labs(title=expression("Posterior Density of " * delta), x=expression(delta))
  
  delta_trace <-
    ggplot(data=df_params) +
    geom_line(aes(x=n, y=delta)) +
    labs(title=expression("Traceplot of " * delta), y=expression(delta))
  
  diagnostic_plots <- list(
    resid_vs_fitted=resid_vs_fitted,
    predicted_vs_actual=predicted_vs_actual,
    eta_density=eta_density,
    eta_trace=eta_trace,
    delta_density=delta_density,
    delta_trace=delta_trace
  )
  
  if(contour_plots){
    interp_y <- akima::interp(x$coords[,1], x$coords[,2], y, nx=nx, ny=ny)
    interp_df <- na.omit(as.data.frame(akima::interp2xyz(interp_y)))
    
    contour <- 
      ggplot(data=interp_df, aes(x=x, y=y, z=z)) +
      geom_contour(col="black", bins=nlevels) +
      labs(
        x=colnames(x$coords)[1],
        y=colnames(x$coords)[2],
        title="Contour plot"
      )
    
    filled_contour <- 
      ggplot(data=interp_df, aes(x=x, y=y, fill=z)) +
      geom_tile() + 
      scale_fill_distiller(...) +
      labs(
        x=colnames(x$coords)[1],
        y=colnames(x$coords)[2],
        title="Filled contour plot",
        fill=x$dep_var
      ) +
      theme(legend.position = "bottom")
    
    legend <- cowplot::get_legend(filled_contour)
    filled_contour_no_legend <- filled_contour + theme(legend.position = "none")
    
    diagnostic_plots_grid <- c(
      diagnostic_plots, 
      list(contour=contour, filled_contour=filled_contour_no_legend)
    )
    
    diagnostic_plots_return <- c(
      diagnostic_plots, list(contour=contour, filled_contour=filled_contour)
    )
    
    plots <- cowplot::plot_grid(plotlist=diagnostic_plots_grid, ncol=2)
    grid1 <- cowplot::plot_grid(plots, legend, rel_heights=c(1,0.12), ncol=1)
  } else {
    diagnostic_plots_grid <- diagnostic_plots_return <- diagnostic_plots
    grid1 <- cowplot::plot_grid(plotlist=diagnostic_plots_grid, ncol=2)
  }
  print(grid1)
  
  covariates_plots <- list()
  for(i in 1:nrow(x$covariates_posterior)) {
    cov_name <- rownames(x$covariates_posterior)[i]
    df_i <- data.frame(values=x$covariates_posterior[i,])
    df_i$n <- seq.int(nrow(df_i))
    
    trace_plot <- 
      ggplot(data=df_i, aes(x=n, y=values)) +
      geom_line() +
      labs(y=cov_name)
    
    density_plot <- 
      ggplot(data=df_i) +
      geom_density(aes(x=values)) +
      labs(x=cov_name)
    
    new_plots <- list(density_plot=density_plot, trace_plot=trace_plot)
    names(new_plots) <- paste(cov_name, names(new_plots), sep="_")
    
    covariates_plots <- c(covariates_plots, new_plots)
  }
  grid2 <- cowplot::plot_grid(plotlist=covariates_plots, ncol=2)
  print(grid2)
  
  res <- list(
    diagnostic_plots = diagnostic_plots_return,
    covariates_plots = covariates_plots,
    diagnostic_grid = grid1,
    covariates_grid = grid2
  )
  return(invisible(res))
  
  
}






















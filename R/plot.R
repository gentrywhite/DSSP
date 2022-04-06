plot.dsspMod <- function(x,
                         robust_residuals=TRUE,
                         panel=if(add.smooth) panel.smooth else points,
                         add.smooth=getOption("add.smooth"),
                         contour_plots=TRUE,
                         nx=100, ny=100, cols=hcl.colors(10, "YlOrRd"), nlevels=5,
                         ...){
  if (!inherits(x, "dsspMod"))	stop("use only with \"dsspMod\" objects")
  
  r <- residuals.dsspMod(x, robust=robust_residuals)
  ylim <- range(r)
  ypred <- predict.dsspMod(x)
  metric <- ifelse(robust_residuals, stats::median, mean)
  yh <- apply(ypred, 1, metric)
  y <- reverse_scaling(x$Y, x$y_scaling)
  
  par(mar=c(1,1,1,1))
  dev.hold()
  plot(yh, r, xlab = "Fitted values", ylab = "Residuals", ylim = ylim)
  panel(yh, r)
  dev.flush()
  
  dev.hold()
  plot(y, yh, xlab = "Actual", ylab = "Predicted", ylim=range(yh))
  abline(0,1)
  dev.flush()
  
  dev.hold()
  plot(stats::density(x$eta), 
       main = expression("Posterior Density of " * eta),
       xlab = expression(eta), ylab = "Posterior density")
  dev.flush()
  
  dev.hold()
  plot(stats::density(x$delta), 
       main = expression("Posterior Density of " * delta),
       xlab = expression(delta), ylab = "Posterior density")
  dev.flush()
  
  if(contour_plots){
    interp_y <- akima::interp(x$coords[,1], x$coords[,2], y, nx=nx, ny=ny)
    dev.hold()
    contour(interp_y)
    dev.flush()
    
    dev.hold()
    filled.contour(interp_y, nlevels=nlevels)
    dev.flush()
  }
  
  invisible()
}

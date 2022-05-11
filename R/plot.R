plot.dsspMod <- function(x,
                         robust_residuals=TRUE,
                         panel=if(add.smooth) graphics::panel.smooth else graphics::points,
                         add.smooth=getOption("add.smooth"),
                         contour_plots=TRUE,
                         nx=100, ny=100, pal=grDevices::heat.colors, nlevels=5,
                         ...){
  if (!inherits(x, "dsspMod"))	stop("use only with \"dsspMod\" objects")
  
  r <- residuals.dsspMod(x, robust=robust_residuals)
  ylim <- range(r)
  ypred <- predict.dsspMod(x)
  metric <- ifelse(robust_residuals, stats::median, mean)
  yh <- apply(ypred, 1, metric)
  y <- reverse_scaling(x$Y, x$y_scaling)
  
  graphics::par(mfrow=c(2,2))
  grDevices::dev.hold()
  plot(yh, r, xlab = "Fitted values", ylab = "Residuals", ylim = ylim)
  panel(yh, r)
  grDevices::dev.flush()
  
  grDevices::dev.hold()
  plot(y, yh, xlab = "Actual", ylab = "Predicted", ylim=range(yh))
  graphics::abline(0,1)
  grDevices::dev.flush()
  
  grDevices::dev.hold()
  plot(stats::density(x$eta), 
       main = expression("Posterior Density of " * eta),
       xlab = expression(eta), ylab = "Posterior density")
  grDevices::dev.flush()
  
  grDevices::dev.hold()
  plot(stats::density(x$delta), 
       main = expression("Posterior Density of " * delta),
       xlab = expression(delta), ylab = "Posterior density")
  grDevices::dev.flush()
  
  if(contour_plots){
    graphics::par(mfrow=c(1,2))
    interp_y <- akima::interp(x$coords[,1], x$coords[,2], y, nx=nx, ny=ny)
    
    grDevices::dev.hold()
    graphics::contour(interp_y)
    grDevices::dev.flush()
    graphics::title("Contour and filled contour plots", outer = T, line=-2)
    grDevices::dev.hold()
    filled.contour2(interp_y, color.palette=pal)
    grDevices::dev.flush()
  }
  invisible()
}


filled.contour2 <-function (x = seq(0, 1, length.out = nrow(z)),
                            y = seq(0, 1, length.out = ncol(z)), 
                            z, xlim = range(x, finite = TRUE), 
                            ylim = range(y, finite = TRUE), 
                            zlim = range(z, finite = TRUE), 
                            levels = pretty(zlim, nlevels), nlevels = 20,
                            color.palette = grDevices::heat.colors, 
                            col = color.palette(length(levels) - 1), plot.title,
                            plot.axes, key.title, key.axes, asp = NA,
                            xaxs = "i", yaxs = "i", las = 1, axes = TRUE, 
                            frame.plot = axes, mar, ...) 
  {
    # taken from here: https://www.cbr.washington.edu/qerm/index.php/R/Contour_Plots
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- graphics::par(c("mar", "las", "mfrow")))$mar
    on.exit(graphics::par(par.orig))
    w <- (3 + mar.orig[2]) * graphics::par("csi") * 2.54
    graphics::par(las = las)
    mar <- mar.orig
    graphics::plot.new()
    graphics::par(mar=mar)
    graphics::plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
      stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
      storage.mode(z) <- "double"
    graphics::.filled.contour(as.double(x), as.double(y), z, as.double(levels), col = col)
    if (missing(plot.axes)) {
      if (axes) {
        graphics::title(main = "", xlab = "", ylab = "")
        graphics::Axis(x, side = 1)
        graphics::Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) 
      graphics::box()
    if (missing(plot.title)) 
      graphics::title(...)
    else plot.title
    invisible()
  }

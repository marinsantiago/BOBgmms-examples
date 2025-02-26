# ------------------------------------------------------------------------------
# 2D Density plots
# ------------------------------------------------------------------------------

library(MASS, include.only = "kde2d")
library(RColorBrewer, include.only = "brewer.pal")

dens_plot <- function(x, y, main, xlab, ylab, col.palette = "Blues",
                      xlim = c(-10, 10), ylim = c(-10, 10),
                      x.points = NULL, y.points = NULL, z.clust = NULL) {
  
  # Estimate 2D KDE ------------------------------------------------------------
  add_margin <- c(-0.2, 0.2)
  dens2d <- MASS::kde2d(
    x = x, y = y, n = 500,
    lims = c(range(x) + add_margin, range(y) + add_margin)
  )
  x <- dens2d$x
  y <- dens2d$y
  z <- dens2d$z
  rm(dens2d, add_margin)
  
  # Set a color palette --------------------------------------------------------
  color.palette <- c("white", RColorBrewer::brewer.pal(9, col.palette)[-1])
  
  # Initialize the plot --------------------------------------------------------
  par(pty = "s")
  plot(
    1, type = "n", xlim = xlim, ylim = ylim, 
    xlab = xlab, ylab = ylab, main = main
  )
  
  # Plot the 2d density --------------------------------------------------------
  image(x, y, z, col = color.palette, xlim = xlim, ylim = ylim, add = TRUE)
  
  # Make borders better --------------------------------------------------------
  par(new = TRUE)
  plot(
    1, type = "n", xlim = xlim, ylim = ylim, 
    xlab = xlab, ylab = ylab, main = main
  )
  
  # If needed, add scatter-plot ------------------------------------------------
  chck <- c(is.null(x.points), is.null(y.points), is.null(z.clust))
  if (all(chck == FALSE)) {
    clusters  <- levels(z.clust)
    cols.scat <- c("darkorange", "magenta1", "black") # yellow1
    pchs.scat <- c(15, 17, 19)
    for (c in seq_len(length(clusters))) {
      points(
        x = x.points[z.clust == clusters[c]],
        y = y.points[z.clust == clusters[c]],
        col = cols.scat[c], pch = pchs.scat[c]
      )
    }
  }
  rm(chck, x, y, z)
  gc()
  cat("")
}


add_contours <- function(x, y, col = 2) {
  
  # Estimate 2D KDE ------------------------------------------------------------
  add_margin <- c(-0.1, 0.1)
  dens2d <- MASS::kde2d(
    x = x, y = y, n = 500,
    lims = c(range(x) + add_margin, range(y) + add_margin)
  )
  x <- dens2d$x
  y <- dens2d$y
  z <- dens2d$z
  rm(dens2d, add_margin)
  
  # Add the contour ------------------------------------------------------------
  contour(
    x, y, z, drawlabels = F, lwd = 2, col = col,
    lty = 1, method = "flattest", add = TRUE,
    levels = pretty(range(z, finite = TRUE), 10)[2] 
  )
  rm(x, y, z)
  gc()
  cat("")
}


bayes_dens <- function(x, y, col = 2, xlim, ylim, main, xlab, ylab) {
  
  # Estimate 2D KDE ------------------------------------------------------------
  add_margin <- c(-0.1, 0.1)
  dens2d <- MASS::kde2d(
    x = x, y = y, n = 500,
    lims = c(range(x) + add_margin, range(y) + add_margin)
  )
  x <- dens2d$x
  y <- dens2d$y
  z <- dens2d$z
  rm(dens2d, add_margin)
  
  # Initialize the plot --------------------------------------------------------
  par(pty = "s")
  plot(
    1, type = "n", xlim = xlim, ylim = ylim, 
    xlab = xlab, ylab = ylab, main = main
  )
  
  # Add the contour ------------------------------------------------------------
  contour(
    x, y, z, add = TRUE, lwd = 2, col = col,
    method = "flattest", drawlabels = FALSE
  )
  rm(x, y, z)
  gc()
  cat("")
}

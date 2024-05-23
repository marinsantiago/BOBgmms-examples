shhh <- suppressPackageStartupMessages
shhh(library(ggplot2, include.only = c("ggplot", "stat_density2d_filled", "aes")))
shhh(library(RColorBrewer, include.only = "brewer.pal"))
shhh(library(ggpubr, include.only = "get_legend"))
shhh(library(stringr, include.only = "str_match")) 
shhh(library(MASS, include.only = "kde2d"))

densPlot <- function(x, y, main, xlab, ylab, colPlot = "Blues", 
                     xlim = c(-2.5, 2.5),  ylim = c(-2.5, 2.5), 
                     x.points = NULL, y.points = NULL, z.clust = NULL) {
  
  get_reduced_labels <- function(plt) {
    
    legPlot <- ggpubr::get_legend(plt)
    legGrob <- 
      unlist(legPlot$grobs[[1]]$grobs[which(grepl("label",
                                            legPlot$grobs[[1]]$layout$name))])
    labels  <- 
      unname(legGrob[grepl("children.GRID.*text.*label$",
                           names(legGrob))])
    
    unname(sapply(labels, function(l) 
      stringr::str_match(l, ", \\s*(.*?)\\s*]")[,2]))
    
  }
  
  # Get the altitude labels
  Data           <- data.frame(x, y)
  colnames(Data) <- c("x", "y")
  plt.labels     <- ggplot(Data, aes(x = x, y = y)) + 
    stat_density2d_filled(contour_var = "density", bins = 9)
  
  altitude.labels <- get_reduced_labels(plt.labels)
  
  rm(plt.labels); rm(Data)
  
  # Estimate KDE
  dens2d <- MASS::kde2d(x, y, n = 500, lims = c(range(x) + c(-0.2, 0.2), 
                                                range(y) + c(-0.2, 0.2)))
  x <- dens2d$x
  y <- dens2d$y
  z <- dens2d$z
  rm(dens2d)
  
  # Set a color palette
  color.palette <- RColorBrewer::brewer.pal(9, colPlot)
  
  # Initialize the plot
  par(pty = "s")
  plot(1, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       main = main)
  
  # Add background color
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4], col = color.palette[1])
  
  # Plot the 2d density
  image(x, y, z, col = color.palette, xlim = xlim, ylim = ylim, add = T)
  
  # Make borders better
  par(new = T)
  plot(1, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       main = main)
  
  # Add altitude legend
  par(xpd = T)
  legend("topright", inset = c(-0.3, 0.15), legend = rev(altitude.labels),
         col = rev(color.palette), pch = 15, pt.cex = 2, y.intersp = 1.3,
         box.col = 0)
  
  # If needed, add scatter-plot
  if((is.null(x.points) & is.null(y.points) & is.null(z.clust)) == F) {
  
    clusters  <- levels(z.clust)
    cols.scat <- c("darkorange", "magenta1", "black") # yellow1
    pchs.scat <- c(15, 17, 19)
      
    for (c in 1:length(clusters)) {
      
      points(x.points[z.clust == clusters[c]],
             y.points[z.clust == clusters[c]],
             col = cols.scat[c], pch = pchs.scat[c])
      
    }
    
  }
  
}


addContours <- function(x, y, col = 2){
  
  dens2d <- MASS::kde2d(x, y, n = 500, lims = c(range(x) + c(-0.1, 0.1), 
                                                range(y) + c(-0.1, 0.1)))
  x <- dens2d$x
  y <- dens2d$y
  z <- dens2d$z
  rm(dens2d)
  
  contour(x, y, z, drawlabels = F, 
          levels = pretty(range(z, finite = TRUE), 10)[2], 
          lwd = 2, col = 2, lty = 1, method = "flattest", add = T)
  
}


bayesDens <- function(x, y, colour = 2, xlim, ylim, main, xlab, ylab){

  dens2d <- MASS::kde2d(x, y, n = 100, lims = c(range(x) + c(-0.1, 0.1), 
                                                range(y) + c(-0.1, 0.1)))
  xx <- dens2d$x
  yy <- dens2d$y
  zz <- dens2d$z
  rm(dens2d)
  
  plot(1, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
       main = main)
  
  contour(xx, yy, zz, add = TRUE, lwd = 2, col = colour, method = "flattest",
          drawlabels = FALSE)
}


densPlot.reduced <- function(x, y, main, xlab, ylab, colPlot = "Blues", 
                             xlim = c(-2.5, 2.5),  ylim = c(-2.5, 2.5), 
                             x.points = NULL, y.points = NULL, z.clust = NULL) {
  
  # Estimate KDE
  dens2d <- MASS::kde2d(x, y, n = 500, lims = c(range(x) + c(-0.2, 0.2), 
                                                range(y) + c(-0.2, 0.2)))
  x <- dens2d$x
  y <- dens2d$y
  z <- dens2d$z
  rm(dens2d)
  
  # Set a color palette
  color.palette <- c("#FFFFFF", RColorBrewer::brewer.pal(9, colPlot)[-1])

  # Initialize the plot
  par(pty = "s")
  plot(1, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       main = main)
  
  # Plot the 2d density
  image(x, y, z, col = color.palette, xlim = xlim, ylim = ylim, add = T)
  
  # Make borders better
  par(new = T)
  plot(1, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       main = main)
  
  # If needed, add scatter-plot
  if((is.null(x.points) & is.null(y.points) & is.null(z.clust)) == F) {
    
    clusters  <- levels(z.clust)
    cols.scat <- c("darkorange", "magenta1", "black") # yellow1
    pchs.scat <- c(15, 17, 19)
    
    for (c in 1:length(clusters)) {
      
      points(x.points[z.clust == clusters[c]],
             y.points[z.clust == clusters[c]],
             col = cols.scat[c], pch = pchs.scat[c])
      
    }
  }
}


#par(mfrow = c(3, 2))
#par(pty = "s", mar = c(5, 2, 2, 2) + 0.1)

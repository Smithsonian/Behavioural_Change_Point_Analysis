#############################################################
###############################################################
#tweaked plotting functions
bcpa.plot<-function (x, type = c("smooth", "flat")[1], threshold = 3, clusterwidth = 1, 
                     col.cp = rgb(0.5, 0, 0.5, 0.5), pt.cex = 0.5, legend = TRUE, xaxt="n",
                     rho.where = "topleft", mu.where = "nowhere", ...) 
{
  windowsweep <- x
  x <- windowsweep$x
  t.POSIX <- windowsweep$t.POSIX
  t <- windowsweep$t
  ws <- windowsweep$ws
  plot(t.POSIX, x,axes=FALSE, type = "n", ...)
  lines(t.POSIX, x, col = "grey")
  if (type == "smooth") {
    if ("pp.smooth" %in% names(windowsweep)) 
      pp <- windowsweep$pp.smooth
    else pp <- PartitionParameters(windowsweep, type = "smooth")
    GoodBreaks <- ws$Break[ws$Model > 0]
    GoodBreaks <- as.data.frame(table(GoodBreaks))
    GoodBreaks <- data.frame(GoodBreaks, t.POSIX = t.POSIX[match(GoodBreaks[, 
                                                                            1], windowsweep$t)])
    GoodBreaks[, 1] <- as.numeric(as.character(GoodBreaks[, 
                                                          1]))
    GoodBreaks <- GoodBreaks[GoodBreaks$Freq >= threshold, 
                             ]
    abline(v = GoodBreaks[, 3], lwd = GoodBreaks[, 2]/threshold * 
             2, col = col.cp)
    rho.scaled <- pp$rho.hat/max(pp$rho.hat, na.rm = 1)
    rho.int <- round(rho.scaled * 999 + 1)
    palette(topo.colors(1000))
    points(t.POSIX, x, pch = 21, col = "darkgrey", bg = rho.int, 
           cex = pt.cex, lwd = 0.5)
    lines(t.POSIX, pp$mu.hat, lwd = 1.5)
    lines(t.POSIX, pp$mu.hat + pp$s.hat, col = "red", lwd = 1.5)
    lines(t.POSIX, pp$mu.hat - pp$s.hat, col = "red", lwd = 1.5)
    rho.hat <- pp$rho.hat
  }
  if (type == "flat") {
    cpsummary <- ChangePointSummary(windowsweep, clusterwidth = clusterwidth)
    phases <- cpsummary$phases
    breaks <- cpsummary$breaks
    whichphase <- findInterval(t, phases$t0)
    rho.hat <- phases$rho.hat[whichphase]
    rho.int <- round(rho.hat/max(rho.hat, na.rm = TRUE) * 
                       999 + 1)
    palette(topo.colors(1000))
    points(t.POSIX, x, pch = 21, col = "darkgrey", bg = rho.int, 
           cex = pt.cex, lwd = 0.5)
    closematch <- rep(NA, length = nrow(phases))
    for (i in 1:nrow(phases)) closematch[i] <- which(abs(t - 
                                                           phases$t0[i]) == min(abs(t - phases$t0[i])))[1]
    phases$t0.POSIX <- t.POSIX[closematch]
    phases$t1.POSIX <- t.POSIX[c(closematch[-1], length(t))]
    t.mid <- (windowsweep$t[-1] + windowsweep$t[-length(windowsweep$t)])/2
    segments(phases$t0.POSIX, phases$mu.hat, phases$t1.POSIX, 
             phases$mu.hat, lwd = 1.5)
    segments(phases$t0.POSIX, phases$mu.hat - phases$s.hat, 
             phases$t1.POSIX, phases$mu.hat - phases$s.hat, col = "red", 
             lwd = 1.5)
    segments(phases$t0.POSIX, phases$mu.hat + phases$s.hat, 
             phases$t1.POSIX, phases$mu.hat + phases$s.hat, col = "red", 
             lwd = 1.5)
    abline(v = phases$t0.POSIX[-1], lwd = breaks$size/max(breaks$size) * 
             4, col = col.cp)
  }
  if (legend) {
    legend.cols <- topo.colors(1000)[seq(1, 1000, length = 5)]
    legend.rhos <- signif(seq(0, max(rho.hat, na.rm = TRUE), 
                              length = 5), 2)
    if (rho.where != "nowhere") 
      legend(rho.where, bg = "white", legend = c(expression(hat(rho)), 
                                                 legend.rhos),cex=0.8, bty="o",pch = 19, ncol = 3, col = c(0, 
                                                                                                           legend.cols), xjust = 0.2, yjust = 0.2)
    if (mu.where != "nowhere") 
      legend(mu.where, bg = "white", legend = c(expression(hat(mu)), 
                                                expression(hat(mu) %+-% hat(sigma))), lty = 1, 
             lwd = 2:1, col = c("black", "red"), xjust = 0.5, 
             yjust = 0.5, cex=0.8, bty="o")
  }
  palette("default")
  x<-round(x,3)
  
  axis.POSIXct(side=1, t.POSIX, format="%H:%M")
  axis(side=2,at=seq(min(x),max(x), by=0.001 ))
 
  
  box()
  
}


###############################################################

#tweaked PathPlot function to display colour by Persistence velocity instead of autocorrelation
PathPlot.Persistence<-function (Data, windowsweep, type = c("smooth", "flat")[1], clusterwidth = 1, 
                                plotlegend = TRUE, tauwhere = "topright", n.legend = 5, ncol.legend = 1, 
                                bty.legend = "n", ...) 
{
  if (!("Z" %in% names(Data))) 
    z <- Data$X + (0+1i) * Data$Y
  else z <- Data$Z
  if (type == "flat") 
    pp <- PartitionParameters(windowsweep, type = type, clusterwidth = clusterwidth)
  if (type == "smooth") 
    if ("pp.smooth" %in% names(windowsweep)) 
      pp <- windowsweep$pp.smooth
    else pp <- PartitionParameters(windowsweep, type = type)
    Segments <- function(z, col = col, lwd = lwd) {
      n <- length(z)
      segments(Re(z[-n]), Im(z[-n]), Re(z[-1]), Im(z[-1]), 
               col = col, lwd = lwd)
    }
    mu.hat <- pp$mu.hat
    rho.hat <- pp$rho.hat
    rho.max <- max(rho.hat, na.rm = 1)
    rho.int <- round(rho.hat/rho.max * 999 + 1)
    mu.max <- max(mu.hat, na.rm = 1)
    mu.int <- abs(round(mu.hat/mu.max * 999 + 1))
    
    
    palette(rev(heat.colors(1000)))
    plot(z, asp = 1, pch = 19, cex = 0.5, col = "grey", ylab="", ...)
    title (ylab="Latitude",xlab="Longitude")
    points(z[c(1, length(z))], pch = c(24, 23), bg = c("green", "blue"), cex = 1, lwd = 1.5, col = "darkgrey")
    
    
    
    Segments(z, col = mu.int, lwd = abs(mu.hat/max(mu.hat, na.rm = TRUE)) *4)
    if (plotlegend) 
      legend(tauwhere, lty = 1, title = ("Persistence Velocity"), 
             ncol = ncol.legend, bty = bty.legend, lwd = 2, col = seq(0, 
                                                                      999, length = n.legend) + 1, legend = signif(seq(0, 
                                                                                                                       max(mu.hat), length = n.legend), 2)) 
    palette("default")
}



#########################################################


#tweaked PathPlot function to display colour by variance in persistence velocity instead of autocorrelation
PathPlot.Variance<-function (Data, windowsweep, type = c("smooth", "flat")[1], clusterwidth = 1, 
                             plotlegend = TRUE, tauwhere = "topright", n.legend = 5, ncol.legend = 1, 
                             bty.legend = "n", ...) 
{
  if (!("Z" %in% names(Data))) 
    z <- Data$X + (0+1i) * Data$Y
  else z <- Data$Z
  if (type == "flat") 
    pp <- PartitionParameters(windowsweep, type = type, clusterwidth = clusterwidth)
  if (type == "smooth") 
    if ("pp.smooth" %in% names(windowsweep)) 
      pp <- windowsweep$pp.smooth
    else pp <- PartitionParameters(windowsweep, type = type)
    Segments <- function(z, col = col, lwd = lwd) {
      n <- length(z)
      segments(Re(z[-n]), Im(z[-n]), Re(z[-1]), Im(z[-1]), 
               col = col, lwd = lwd)
    }
    s.hat <- pp$s.hat
    rho.hat <- pp$rho.hat
    rho.max <- max(rho.hat, na.rm = 1)
    rho.int <- round(rho.hat/rho.max * 999 + 1)
    s.max <- max(s.hat, na.rm = 1)
    s.int <- abs(round(s.hat/s.max * 999 + 1))
    
    
    palette(rev(rainbow(16)))
    plot(z, asp = 1, pch = 19, cex = 0.5, col = "grey", ...)
    points(z[c(1, length(z))], pch = c(24, 23), bg = c("green", "blue"), cex = 1, lwd = 1.5, col = "darkgrey")
    
    
    
    Segments(z, col = s.int, lwd = 2)
    if (plotlegend) 
      legend(tauwhere, lty = 1, title = ("Variance in PV"), 
             ncol = ncol.legend, bty = bty.legend, lwd = 2, col = seq(0, 
                                                                      999, length = n.legend) + 1, legend = signif(seq(0, 
                                                                                                                       max(s.hat), length = n.legend), 2)) 
    palette("default")
}



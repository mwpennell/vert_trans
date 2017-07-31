library(scales)
library(ggplot2)

## colours 
cols <- c("#053061", "#78B7C5", "firebrick3", "#EE6363", "#0B775E", "#D8B70A",
          "#999999", "#252A30")
names(cols) <- c("mhet", "mhom", "fhet", "fhom", "esd", "herm", "lg", "dg")

add.profile.shading <- diversitree:::add.profile.shading

profile.block <- function(y, col.line, col.fill, col.l, col.g,
                          xlim=NULL, ymax=NULL,
                          a.fill=0.5, a.bg=0.4, lwd=2,
                          xlab="Parameter estimate",
                          ylab="Posterior density",...){
  if (missing(xlim))
    xlim <- range(y)
  dx <- diff(xlim)/(49)
  f <- function(yi){
    n <- (range(yi) - xlim)/dx
    n <- c(floor(n[1]), ceiling(n[2]))
    ri <- xlim + n * dx
    hist(yi, seq(ri[1], ri[2], by = dx), plot = FALSE)
  }
  hh <- f(y)
  ci <- hdr(y)
  if (is.null(ymax))
    ymax <- max(hh$density)
  
  par(mar=c(5,5,4,2))
  ylim <- c(0, 1.05 * ymax)
  plot(NA, xlim=xlim, ylim=ylim, type="n", yaxt="n",
       xlab=xlab, ylab=ylab, ...)
  
  lim <- par("usr")
  rect(lim[1]-1, lim[3]-1, 0, lim[4]+1, border = alpha(col.l, a.bg),
       col = alpha(col.l,a.bg))
  rect(0, lim[3]-1, lim[2]+1, lim[4]+1, border = alpha(col.g,a.bg),
       col = alpha(col.g, a.bg))
  hist.outline(hh, col=col.line, lwd=lwd)
  add.profile.shading(hh, ci, alpha(col.fill, a.fill))
}

## Draw the outline of a histogram
hist.outline <- function(h, ..., density=TRUE) {
  xy <- hist.xy(h, density)
  lines(xy, ...)
}
hist.fill <- function(h, ..., density=TRUE) {
  xy <- hist.xy(h, density)
  polygon(xy, ...)
}

hist.xy <- function(h, density=TRUE) {
  dx <- diff(h$mids[1:2])
  xx <- rep(with(h, c(mids - dx/2, mids[length(mids)] + 
                        dx/2)), each = 2)
  yy <- c(0, rep(if (density) h$density else h$counts, each = 2), 0)
  list(x=xx, y=yy)
}





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

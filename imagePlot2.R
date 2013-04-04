if(!require("fields")){
    stop("fields package must be installed")
}


image.plot2 <- function (..., add = FALSE, nlevel = 64, horizontal = FALSE, 
                         legend.shrink = 0.9, legend.width = 1.2,
                         legend.mar = ifelse(horizontal,
                         3.1, 5.1), legend.lab = NULL,
                         graphics.reset = FALSE, bigplot = NULL,
                         smallplot = NULL, legend.only = FALSE,
                         col = tim.colors(nlevel), 
                         lab.breaks = NULL, axis.args = NULL,
                         legend.args = NULL, 
                         midpoint = FALSE,
                         cexLabAxis=1.5) {
    old.par <- par(no.readonly = TRUE)
    info <- image.plot.info(...)
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    temp <- image.plot.plt(add = add, legend.shrink = legend.shrink, 
                           legend.width = legend.width,
                           legend.mar = legend.mar, 
                           horizontal = horizontal,
                           bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        if (!info$poly.grid) {
            image(..., add = add, col = col)
        }
        else {
            poly.image(..., add = add, col = col, midpoint = midpoint)
        }
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    ix <- 1
    minz <- info$zlim[1]
    maxz <- info$zlim[2]
    binwidth <- (maxz - minz)/nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    breaks <- list(...)$breaks
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    if (!is.null(breaks) & !is.null(lab.breaks)) {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                            mgp = c(3, 1, 0),
                            las = ifelse(horizontal, 0, 2), 
                            at = breaks, labels = lab.breaks),
                       axis.args)
    }
    else {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)),
                       cex=cexLabAxis, cex.axis=cexLabAxis,
                       axis.args)
    }
    ##axisCex <- match.call(expand.dots = )$cex.axis
    ##if(is.null(axisCex)) {
    ##    axisCex <- 1.0
    ##}
    
    if (!horizontal) {
        if (is.null(breaks)) {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                  ylab = "", col = col)##, cex.axis=1.5)
        }
        else {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col, breaks = breaks)##, cex.axis=axisCex)
        }
    }
    else {
        if (is.null(breaks)) {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
                  ylab = "", col = col)##, cex.axis=axisCex)
        }
        else {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col, breaks = breaks)##, cex.axis=axisCex )
        }
    }
    do.call("axis", axis.args)
    box()
    if (!is.null(legend.lab)) {
        legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
                                               1, 4), ##line = legend.mar - 2,
                            ## change mar as cex changes
                            line=legend.mar -.5,
                            cex=cexLabAxis)
    }
    if (!is.null(legend.args)) {
        do.call(mtext, legend.args)
    }
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}


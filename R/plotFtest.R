
## allows user ot over ride y axis label.
## this is an unresolved issue in multitaper package
# the issue happens as .plotFtest is called by a parent function.
plotFtest <- function(x, 
                       ftbase=1.01, 
                       siglines=NULL, 
                       xlab="Frequency",
                       ylab="Harmonic F-test Statistic",
                       ...) {

    if(is.null(x$mtm$Ftest) || !("Ftest" %in% class(x))) {
      stop(paste("Ftest not computed for given mtm object!"))
    }
   
    log <- match.call(expand.dots = )$log
    ylab <- match.call(expand.dots = )$ylab

    ##if(is.null(ylab)) ylab <- "Harmonic F-test Statistic"
    
    ylog = "n"
    if(is.null(log) || log == "yes") {
        ylog = "y"
    }

    ftestVals = x$mtm$Ftest
    ftestVals[ftestVals < ftbase] <- ftbase
    ftmax <- max(ftestVals)

    plot(x$freq, ftestVals, log=ylog, ylab=ylab, xlab=xlab, 
         ylim=c(ftbase,ftmax), type="l",...)
    
    ## add siglines if defined
    if(!is.null(siglines)) {
        for(j in 1:length(siglines)) {
            if(is.numeric(siglines[j]) && 0.80 <= siglines[j] && 1.000000 >= siglines[j]) {
                ## degree of freedom correction P&W page 499 changed to 2, 2*k-2 date Sept 30 2012
                sig0 <- qf(siglines[j],2,2*x$mtm$k-2) 
                abline(h=sig0, col="red", lty=2, lwd=1) 
                mtext(paste(siglines[j]*100,"%",sep=""), side=4, line=0, at=ceiling(sig0), col="red")
            }
        } ## end for
    } ## end logical
}

## some tools based on the SAPA lisp code
## Spectral analysis for physical application
## Percival, D.B. and Walden, A.T. (1993)
## this is based on the multitapaer.lisp file.
## http://faculty.washington.edu/dbp/sapabook.html


mtmCI <- function(mtmObj, ciWidth=.95) {
    stopifnot(mtmObj$mtm != NULL)
    stopifnot(ciWidth < 1 && ciWidth > 0)

    lowerSpec <- mtmObj$spec
    upperSpec <- mtmObj$spec
    dofs <- mtmObj$mtm$dofs
    plower <- (1 - ciWidth)/2
    pupper <- 1 - plower

    lowerSpec <- lowerSpec * dofs/qchisq(plower, dofs, lower.tail=FALSE)
    upperSpec <- upperSpec * dofs/qchisq(pupper, dofs, lower.tail=FALSE)

    list(upperSpec=upperSpec, lowerSpec=lowerSpec)
}

## example of use
## data set not provided
## mtm1 <-  spec.mtm(someTS, kVal, nw=nwVal, nFFT=nFFT,
##                         Ftest=TRUE,returnZeroFreq=TRUE, plot=FALSE)
## plot(mt1, log="y")
## mtmCI1 <- mtmCI(mtm1)
## lines(mtm1$freq, mtmCI1$upperSpec, lty=2)
## lines(mtm1$freq, mtmCI1$lowerSpec, lty=2)



## to use this, select a power value on the graph and draw a line, you may also draw a cross
## indicating ci values and bandwidth (2w).
## this cross may look appropriate on the top right, bottom left, or under a pronounced peak.
##segments(xVal, powerVal*lowerScale, y1=powerVal*upperScale, col="red")

## valid values for fun1, are mean, median, min, and max.
## this will use the median, mean, min or max degrees of freedom
## the min degrees of freedom will give the most conservative CI
##
mtmCIOneValScale <- function(mtmObj, ciWidth=.95, fun1=median) {
    stopifnot(mtmObj$mtm != NULL)
    stopifnot(ciWidth < 1 && ciWidth > 0)
    
    
    dofs <- mtmObj$mtm$dofs

    dof1 <- fun1(dofs)
   
    plower <- (1 - ciWidth)/2
    pupper <- 1 - plower

    lowerScale <- dof1/qchisq(plower, dof1, lower.tail=FALSE)
    upperScale <- dof1/qchisq(pupper, dof1, lower.tail=FALSE)

    list(upperScale=upperScale, lowerScale=lowerScale)
}


drawMTMCI <- function(f0, s0, mtmObj, ciWidth=.95, fun1=median, ...) {
    res1 <- mtmCIOneValScale(mtmObj, ciWidth=.95, fun1=median)
    s1 <- s0*res1$lowerScale
    s2 <- s0*res1$upperScale

    segments(f0, s1, y1=s2, ...)
}
    

draw2WlineAt <- function(f0, s0, mtmObj, subtractOneFromLeft=TRUE, ...) {
    w <- mtmObj$mtm$nw/(mtmObj$origin.n*mtmObj$mtm$deltaT)

    ## hight is not exactly centred on s0, but this is the correct width
    ## based on the location of the frequency.
    ## the upper ci will appear longer
    ## this is consistent with the theory.
    
    ## we can subtract one bin from the left bin in drawing 2w, so the
    ## line is not quite centred at f0, it is not completly balanced.
    
    ## get freqBin, note zeroth freq may not be returned.
    freqBin <- 0 
    if(subtractOneFromLeft) {
        freqBin <- mtmObj$freq[2] - mtmObj$freq[1]
    }
    
    f1 <- f0-w
    f2 <- f0+w - freqBin

    segments(f1, s0, x1=f2, ...)
}


drawMTMCross <- function(f0, s0, mtmObj, ciWidth=.95, fun1=median,
                         subtractOneFromLeft=TRUE,
                         ...) {
    drawMTMCI(f0, s0, mtmObj, ciWidth=ciWidth, fun1=fun1, ...)
    draw2WlineAt(f0, s0, mtmObj,
                 subtractOneFromLeft=subtractOneFromLeft,
                 ...)
}

## bandwidth on bartlett and spectrograms
bandwidthSegOnSpectrogram <- function(w, f0, idx) {
    flo <- f0 - w
    fhi <- f0 + w
    segments(idx, flo, idx, fhi)
}

bandwidthSegOnBartlet <- function(w, f0, yVal) {
    flo <- f0 - w
    fhi <- f0 + w
    segments(flo, yVal, fhi, yVal)
}

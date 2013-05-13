require(multitaper)

## the coherences in the R package example are way outside what I expected
##

## helper function required
plotJkcoh1 <- function(freqs,TRmsc, NTvar, k, nfreqs, percentGreater=NULL,
                       nehlim=10, nehc=4,
                       cdfQuantiles=c(0.001, 0.002, 0.005, 0.010, 0.020,
                       0.050, 0.100, 0.200, 0.500, 0.800, 0.900, 0.950,
                       0.980, 0.990, 0.995, 0.998, 0.999),
                       mscTicks=seq(0, .9, .1),
                       drawPercentLines=TRUE,
                       percentG=c(.1,.2,.5,.8,.9)) {
    
    ##nehlim and nehc are for smoothing 
    ## currently we plot the smoothed transformed coherence
    ## and lower CI after smoothing the variance
    plotTRmsc <- multitaper:::.lftr3p(TRmsc, NTvar, nfreqs,
                                      nehlim,nehc, "even", "ext")
    trnrm_ <- multitaper:::.trnrm(k)
    par(oma=c(2,4,0,2))
    plot.new()
    plot.window(range(freqs), range(plotTRmsc[,2]))
    xy <- xy.coords(freqs,plotTRmsc[,2])
    plot.xy(xy, type="l", lwd=1)
    lines(freqs, plotTRmsc[,1], lty=3, lwd=1)
    box()
    axis(1)
    mtext("Frequency", side=1, line=3, cex=par()$cex)
    TRmscTicks <- seq(0, max(plotTRmsc[,2]), .5)
    axis(2, at=TRmscTicks)
    mtext("Inverse Transform of Magnitude Squared Coherence",
          side=2, line=2, cex=par()$cex)
    
    ##  outer MSC axis on the left
    msc <- multitaper:::.FtoMSC(plotTRmsc[,2], trnrm_)
    maxMSC <- max(msc)
    ## ########################################################################
    ## # this is also causing problems in the example provided in the package
    ##mscTicks <- seq(0,maxMSC, .1)
    
    TRmscTicks <- multitaper:::.C2toF(mscTicks, trnrm_)
    axis(2, at=TRmscTicks, labels=mscTicks, outer=T)
    mtext("Magnitude Squared Coherence", side=2, line=6, cex=par()$cex)
    
    ## right cdf axis
    ## this function still as the name DJT gave it
    ## multitaper:::.paxpt7()$out
    ## [1] 0.001 0.002 0.005 0.010 0.020 0.050 0.100 0.200
    ## 0.500 0.800 0.900 0.950
    ## [13] 0.980 0.990 0.995 0.998 0.999
    ## CDFT <- multitaper:::.paxpt7()$out
    ## CDFT <- 1-10^(-(6:12))
    CDFT <- cdfQuantiles

    Qlvl <- multitaper:::.cdfToMSqCoh(CDFT, k)
    Qlvl <- Qlvl[Qlvl <= maxMSC]
    ## ####################################################################################3
    ## The following is the problem with the mtm.coh example in version 1.04
    ## msc <- multitaper:::.FtoMSC(plotTRmsc[,2], trnrm_)
    ## mscLB <-  multitaper:::.FtoMSC(plotTRmsc[,1], trnrm_)
    ## plot(freqs, msc, type="l")
    ## lines(freqs, mscLB, lty=2)
    
    ##      Qlvl >= min(msc)
    ##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [13] FALSE FALSE FALSE FALSE FALSE
    
    TRQlvl <- multitaper:::.C2toF(Qlvl, trnrm_)
    lenLessThanMax <-  length(Qlvl)
    CDFT <- CDFT[1:lenLessThanMax]
    axis(4, at=TRQlvl, labels=CDFT)
    mtext("Cumulative Distribution Function for Independent Data",
          side=4, line=2, cex=par()$cex) 
    
    if(drawPercentLines == T) {
        percentG <- multitaper:::.C2toF(multitaper:::.cdfToMSqCoh(percentG, k),  trnrm_)
        lenPercentG <- length(percentG)
        for(i in 1:lenPercentG) {
            lines(freqs, array(percentG[i], nfreqs), lty=2)
        }
    }
    
    if(!is.null(percentGreater)) {
        mtext(paste("CDF for C=   10.0% 20.0% 50.0% 80.0% 90.0%"),
              side=1, line=4, adj=-1, cex=.8)
        mtext(paste("% of data > Q     ",
                    100*round( percentGreater[1], digits=3),
                    "% ",
                    100*round( percentGreater[2], digits=3),
                    "% ",
                    100*round( percentGreater[3], digits=3),
                    "% ",
                    100*round( percentGreater[4], digits=3),
                    "% ",
                    100*round( percentGreater[5], digits=3),
                    "%", sep=""),
              side=1, line=5, adj=-1, cex=0.8)
    }
}



## these are the modified versions of the functions, and one helper that I plan to place in the multitaper package


## attempt to make a nicer function
mscToCDFquantiles <- function(msc, k) {
    1 - (1-msc)^(k-1)
}

plotJkcoh2 <- function(freqs, TRmsc, NTvar, k, nfreqs, percentGreater=NULL,
                       nehlim=10, nehc=4,
                       cdfQuantilesTicks=NULL,
                       drawPercentLines=TRUE,
                       percentG=c(.1,.2,.5,.8,.9)) {

    ##freqs is the frequencies required
    ## TRmsc is the inverse transformed msc
    ## NTvar is one standard deviation jackkife variance.

    ## smoothing
    ##nehlim and nehc are for smoothing
    plotTRmsc <- multitaper:::.lftr3p(TRmsc, NTvar, nfreqs,
                       nehlim,nehc, "even", "ext")

    ## trnrm_ is 2k - 2
    trnrm_ <- multitaper:::.trnrm(k)

    ## set up fancy ploting dimensions
    par(oma=c(2,4,0,2))
    plot.new()
    plot.window(range(freqs), range(plotTRmsc[,2]))
    ## plot smoothed msc
    xy <- xy.coords(freqs,plotTRmsc[,2])
    plot.xy(xy, type="l", lwd=1)
    ## plot one sd dev lower jackknife variance
    lines(freqs, plotTRmsc[,1], lty=3, lwd=1)
    box()
    axis(1)
    mtext("Frequency", side=1, line=3, cex=par()$cex)

    ## basic left axis
    axis(2) 
    mtext("Inverse Transform of Magnitude Squared Coherence",
          side=2, line=2)

    ##  outer MSC axis on the left
    ## get msc and ticks
    msc <- multitaper:::.FtoMSC(plotTRmsc[,2], trnrm_)
    mscTicks <- pretty(msc)

    ## transform ticks for at
    ##C2toF is coherence to inverse transform
    TRmscTicks <- multitaper:::.C2toF(mscTicks, trnrm_)
    axis(2, at=TRmscTicks, labels=mscTicks, outer=T)
    mtext("Magnitude Squared Coherence", side=2, line=6, cex=par()$cex)

    ##mscToCDF values may have issues for highly coherent values
    ## values over .9 will cause issues
    if(is.null(cdfQuantilesTicks)) {
        cdfQuantiles <- mscToCDFquantiles(msc, k)
        cdfQuantilesTicks <- pretty(cdfQuantiles)
    }

    ## put right axis
    Qlvl <- multitaper:::.cdfToMSqCoh(cdfQuantilesTicks, k)
    TRQlvl <- multitaper:::.C2toF(Qlvl, trnrm_)
    
    cumulativeDistVals <- multitaper:::.C2toF(msc, trnrm_)
    axis(4, at=TRQlvl, labels=cdfQuantilesTicks)
    mtext("Cumulative Distribution Function for Independent Data",
          side=4, line=2, cex=par()$cex) 
    
    if(drawPercentLines == TRUE) {
        percentG <- multitaper:::.C2toF(multitaper:::.cdfToMSqCoh(percentG, k),  trnrm_)
        lenPercentG <- length(percentG)
        for(i in 1:lenPercentG) {
            lines(freqs, array(percentG[i], nfreqs), lty=2)
        }
    }
    
    if(!is.null(percentGreater)) {
        mtext(paste("CDF for C=   10.0% 20.0% 50.0% 80.0% 90.0%"),
              side=1, line=4, adj=-1, cex=.8)
        mtext(paste("% of data > Q     ",
                    100*round( percentGreater[1], digits=3),
                    "% ",
                    100*round( percentGreater[2], digits=3),
                    "% ",
                    100*round( percentGreater[3], digits=3),
                    "% ",
                    100*round( percentGreater[4], digits=3),
                    "% ",
                    100*round( percentGreater[5], digits=3),
                    "%", sep=""),
              side=1, line=5, adj=-1, cex=0.8)
    }
}


###########################################################################################
## The following represents rewritten mtm.coh example
###

## Data(HadCRUTnh)
## data(mlco2)
## spec1 <- spec.mtm(HadCRUTnh, nw=5.0, k=8, plot=FALSE, returnInternals=TRUE, dtUnits="month", dT=1.0)
## spec2 <- spec.mtm(mlco2, nw=5.0, k=8, plot=FALSE, returnInternals=TRUE, dtUnits="month", dT=1.0)
## cohRes <- mtm.coh(spec1, spec2, plot=FALSE)


## ## example with simple modification
## plotJkcoh1(spec1$freq, cohRes$NTmsc, cohRes$NTvar, cohRes$k, cohRes$nfreqs,
##            cdfQuantiles= 1-10^(-(6:12)), mscTicks=c(.95, .96, .97, .98))

## ## example using R's function pretty
## plotJkcoh2(spec1$freq, cohRes$NTmsc, cohRes$NTvar, cohRes$k, cohRes$nfreqs)
## ## example with R's function pretty allowing cdfQuantiles to be specified.
## plotJkcoh2(spec1$freq, cohRes$NTmsc, cohRes$NTvar, cohRes$k, cohRes$nfreqs,
##           cdfQuantilesTicks=  1-10^(-(6:12)))


################################################################################################
### Below represents some tests and some examples that show how to access the internals
## of the msc function and plot basic msc




## set some defaults
## percentGreater=NULL
## nehlim=10
## nehc=4
## drawPercentLines=TRUE
## percentG=c(.1,.2,.5,.8,.9)

## ## simplify variables

#### set some values for easy access
## freqs <- spec1$freq
## TRmsc <-  cohRes$NTmsc
## NTvar <- cohRes$NTvar
## k <- cohRes$k
## nfreqs <- cohRes$nfreqs


## ## smooth
## plotTRmsc <- multitaper:::.lftr3p(TRmsc, NTvar, nfreqs,
##                                   nehlim,nehc, "even", "ext")

## ## scaling factor
## trnrm_ <- multitaper:::.trnrm(k)

##########################################################################################
## plot one, tansformed msc
## ## plot normal transformed msc
## plot(freqs, TRmsc, type="l")

##############################################################################################
## plot two plot normal transformed msc
## ## plot smoothed normal transformed msc
## plot(freqs, plotTRmsc[,2], type="l")
## ## jackkinfed lower bound value
## lines(freqs, plotTRmsc[,1], lty=2)


############################################################################################
## plot 3 plot msc and jackknife lower bound
## msc <- multitaper:::.FtoMSC(plotTRmsc[,2], trnrm_)
## mscLB <-  multitaper:::.FtoMSC(plotTRmsc[,1], trnrm_)
## plot(freqs, msc, type="l")
## lines(freqs, mscLB, lty=2)




## these options may not work on all systems
## I use them in Linux.

setGoldenRatioPlot <- function() {
    ## the constants come from rounding
    ## and using GR constant phi
    ## phi <- (1 + sqrt(5))/2
    ##
    ## 5*phi
    ## [1] 8.09017
    ## 7/phi
    ## [1] 4.326238

    ## it appears pdf dimensions control the font size
    ## larger pdf dimensions == smaller font size
    ## font size can also be set with cex, cex.axis, cex.lab
    ## and cex.main
    ## these units are in inches.
    X11.options(width=8.09, height=5)
    pdf.options(width=7, height=4.33)
}

resetPlotRatios <- function() {
    X11.options(reset=TRUE)
    pdf.options(reset=TRUE)
}

setGoldenRatioPlot2Col <- function() {
    X11.options(width=8.09, height=5)
    ##   3.5/phi
    ## [1] 2.163119

    pdf.options(width=3.5, height=2.16)
}

set9inchWidthGoldenRatio <- function() {
    ## phi <- (1 + sqrt(5))/2
    ## 9/phi
    ## [1] 5.562306

    X11.options(width=8.09, height=5)
    pdf.options(width=9, height=5.56)
}

##figure in grape harvest section of thesis.
## figure settings
## this makes the images big but the font sizes shrink
## there is a tradeoff between image size and font size
##pdf.options(width=11.34, height=7)
##X11.options(width=11.34, height=7)

## this resests things
##X11.options(reset=TRUE)
##pdf.option(reset=TRUE)

## this seems decent for papers with a width of about 7 inches
## you will need a smaller plot for a half size on a two column paper
##X11.options(width=7, height=4.33)
##pdf.options(width=7, height=4.33)

## untested two column defaults
##X11.options(width=3.5, height=2.16)
##pdf.options(width=3.5, height=2.16)

## margins need not be identical,
## and we can mix and match 
##X11.options(width=11.34, height=7)
##X11.options(width=8.09, height=5)
##pdf.options(width=7, height=4.33)

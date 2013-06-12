examplePercivalWaldenNearPeaksCh10s8e1 <- function(t=(1:256), varError=10^(-10),
                                                phi1=0, phi2=0, phi3=0) {
    ## this is from percival and walden page 486
    ## note that the dpss with nw=2 is similar to the Hanning data taper in this case
    ## also not that this example was designed to favour the nw=4 data taper
    res <-  0.0316 * cos(0.2943 * pi * t + phi1) +
        cos(0.3333 * pi * t + phi2) +
            0.0001 * cos( 0.3971 * pi * t + phi3)
    wn <- rnorm(length(t), sd=sqrt(varError))
    res <- res + wn
    res
}


examplePercivalWaldenOnePeakCh10s8e2 <- function(t=(0:100), f1=7.25,
                                               deltaT=.01,
                                               varError=0.0001,
                                               phi1=pi/4) {
    ## this is from percival and walden page 487--488
    ## in this case the periodgram biases the actual frequency higher than it should
    res <- cos( 2 * pi * f1 * t * deltaT + phi1)
    wn <- rnorm(length(t), sd=sqrt(varError))
    res <- res <-  wn
    res
}

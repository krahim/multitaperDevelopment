if(!require("multitaper")){
    stop("multitaper package must be installed")
}


## based on German Prieto's et. al  2005 IEEE paper
## MULTITAPER WIGNER-VILLE SPECTRUM FOR DETECTING DISPERSIVE SIGNALS
## FROM EARTHQUAKE RECORDS
##  1­4244­0132­1/05/

loeveSpectrum <- function (data, nw, nord, nFFT=length(data), deltat=1.0)  {
    
    res <- spec.mtm(data, ##willamette,
                    nw=nw, k=nord, nFFT=nFFT,
                    dT=deltat,
                    returnInternals=TRUE, plot=FALSE)

    ## with the frequencies we do not get the complex
    ##nord <- 8
    nFreqs <- res$mtm$nfreqs
    eigenCoefs <- res$mtm$eigenCoefs
    dk2 <- res$mtm$eigenCoefWt
    dk <- sqrt(dk2)
    lambdak <- res$mtm$dpss$eigen
    spec <- res$spec
    sqrtSumdk2 <- sqrt(apply(dk2, 1, sum))
    
    gamma <- array(NA, dim=c(nFreqs, nFreqs))
    gamma2 <- array(NA, dim=c(nFreqs, nFreqs))
    ## very bloody slow
    for( f1 in 1:nFreqs) {
        for( f2 in 1:nFreqs) {
            for( k in 1:nord) {
                gamma[f1, f2] <-lambdak[k] * dk[f1, k] *
                Conj(eigenCoefs[f1, k]) * dk[f2, k] *
                    eigenCoefs[f2, k]
                gamma[f1, f2] <-  gamma[f1, f2] /
                    (sqrtSumdk2[f1] * sqrtSumdk2[f2])
            }
    }
    }
    
    
    for( f1 in 1:nFreqs) {
        for( f2 in 1:nFreqs) {
            gamma2[f1,f2] <- abs(gamma[f1,f2])**2/(spec[f1]* spec[f2])
        }
    }

    ## plot this in dB for example image.plot(10*log10(gamma2)
    ## where gamma2 is from fields package
    gamma2
}
    





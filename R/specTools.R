### from rpaper package
##Function to generate theoretical spectrum included in paper plot
## based on pw code. 
## moved form specTools.R in MultitaperDevelopmentPrivate June 20, 2013

## ## this is in parametric.R which should be in the PWLISPR package, moved

## arCoefsToSdf <- function(arCoeffs, innovationsVariance, deltaT, nNonZeroFreqs) {
##     ## from rpaper package
##     ##Function to generate theoretical spectrum included in paper plot
##     ## based on pw code. 

##     prewhitenFilter <- c(1, -arCoeffs)
##     scratch <- c(prewhitenFilter, rep(0.0, 2 *
##                                       (nNonZeroFreqs - 1) -
##                                       length(prewhitenFilter)))
##     scratch <- fft(scratch)
##     numerator <- innovationsVariance * deltaT
##     resultSdf <- numerator / abs(scratch[1:(1+nNonZeroFreqs)])**2
##     return(resultSdf)
## }

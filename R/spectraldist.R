##     The multitaper R package (development code)
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2013 Karim Rahim 
##
##     Written by Karim Rahim.
##
##     This file is part of the multitaper package for R.
##     http://cran.r-project.org/web/packages/multitaper/index.html
## 
##     The multitaper package is free software: you can redistribute it and 
##     or modify it under the terms of the GNU General Public License as 
##     published by the Free Software Foundation, either version 2 of the 
##     License, or any later version.
##
##     The multitaper package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
##
##     If you wish to report bugs please contact the author:
## 
##     Karim Rahim
##     karim.rahim@gmail.com
##     Jeffery Hall, Queen's University, Kingston Ontario
##     Canada, K7L 3N6

library("multitaper")

## Utilities for spectral distance code used with multitaper spectrograms.

spectraldist <- function(spectra1, spectra2) {
    dist1 <- (log(spectra1) - log(spectra2))^2
    ##do we need a better integration technique?
    sumDist <- sum(dist1)
    return(list(dist=dist1, sumDist=sumDist))
}

getMaxSpectralDistAllBlocks <- function(blockSdfs) {
    nBlocks <- dim(blockSdfs)[2]
    res <- array(NA, nBlocks -1)

    for(i in 1:(nBlocks-1) ) {
        res[i] <- spectraldist(blockSdfs[,i], blockSdfs[,i+1])$sumDist
    }
    return(res)

}


findLocalMax <- function(dists, cutoff=0) {

    n <- length(dists)
    distvals <- array(NA, n-1)
    for(i in 1:(n -1) ){
        distvals[i] <- as.numeric((dists[i] < dists[i+1]))
    }

    idxs <- NULL
    maxs <- NULL
    
    for(i in 2:(n-1) ) {
        if( (distvals[i] == 0) && (distvals[i-1] == 1) &&
           (dists[i] > cutoff)) {
            idxs <- c(idxs, i)
            maxs <- c(maxs, dists[i])
        }
    }

    ## endpoints
    if((dists[1] > dists[2]) && (dists[1] > cutoff)) {
        idxs <- c(1, idxs)
        maxs <- c(dists[1], maxs)
    }

    if((dists[n] > dists[n-1]) && (dists[n] > cutoff)) {
        idxs <- c(idxs, n)
        maxs <- c(maxs, dists[n])
    }
    return(list(idxs=idxs, maxs=maxs))
}

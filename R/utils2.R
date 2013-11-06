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

##library("multitaper")

## Utilities for and qi code

sizeplotForRightAxis <- function() {
    par(mar=c(5.1,4.1,4.1,4.1))
}

scaleFnToFit <- function(fnToFit, fnPlotted) {
    const1 <- min(fnPlotted) - min(fnToFit)
    scale1 <- diff(range(fnPlotted)) / diff(range(fnToFit))
    const1+scale1*fnToFit
}

##error checks?
## we don't need no stinking error checks

getFittedAR <- function(dat1, ar1, addMean=FALSE) {

    n1 <- length(dat1)
    order1 <- ar1$order
    dat1 <- dat1 - ar1$x.mean
    fitted1 <- convolve(dat1, rev(ar1$ar), type="filter")[1:(n1-order1)]
    c(rep(NA, order1), fitted1) + if(addMean) ar1$x.mean else 0
}
    

## helper function may be able to use R
## seq(along=x)[x == max(x)]
## may have to condsider duplicates

getIndex <- function(items, searchArray) {
   nItems <- length(items);
   if(nItems==0) {
      return();
   }
   N <- length(searchArray);
   result <- array(NA, nItems);
   for( i in 1:nItems) {
      for(j in 1:N) {
         if(items[i] == searchArray[j]) {
		      result[i] <- j;
		      break;
          }
      }
   }
   return(result);
}

descriptiveStats <- function(X) {
    cm4pp(X)
}


## stats matching those used in djt's fortran libraries.

kurtosis <- function(X) {
   return(sum((X-mean(X))^4)/((length(X)-1) * var(X)^2));
}

kurtosis3 <-  function(X) {
    return(kurtosis(X) -3);
}



cm4pp <-  function(X) {
    n <-  length(X)
    nM1 <-  n -1
    maxX <- max(X)
    minX <- min(X)
    maxI <-  getIndex(maxX, X)
    minI <-  getIndex(minX, X)
    seq <- 1:n
    ## getIndex here returns an array potential error...
    tI <- seq != minI[1] & seq != maxI[1]
    tX <-  X[tI]
    tMean <-  mean(tX)
    tSD <- sd(tX)
    tVar <-  var(tX)
   
    mu <-  mean(X)
    sigma <- sd(X)
    sigma2 <- var(X)
    cm <- (X - mu)
    cm3 <-  cm^3
    cm4 <-  cm^4
    skew1 <-  sum(cm3)/(nM1 * sigma^3)
    kurt <-  sum(cm4)/(nM1 * sigma2^2)

    return(list( maxX=maxX, maxI=maxI, minX=minX,
                minI=minI, avg=mu,
                tAvg=tMean, sigma=sigma, tSigma=tSD,
                sigma2=sigma2, tSigma2=tVar, skew=skew1,
                kurt=kurt))
}

## do we need the edge values....
##stripWatFreqAtEdge(freq, 0, strictInEquality=FALSE)

stripWatFreqAtEdge <- function(freq, w, strictInEquality=TRUE) {
    ## fracW is not used... the idea is to multiply w by a fraction
    ## one can simply multiply nw by this fraction.
    ##stripWatFreqAtEdge(freq, 0, strictInEquality=FALSE) ##to prevent stripping
    nyquist <- freq[length(freq)]
    idx <- NULL
    if(strictInEquality) {
        idx <- (freq > w) & (freq < (nyquist - w ))
    } else {
        idx <- (freq >= w) & (freq <= (nyquist - w ))
    }
        
    list(idx=(1:length(freq))[idx], freq=freq[idx])
}

    

    

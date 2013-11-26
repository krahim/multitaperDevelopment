##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 
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
##     112 Jeffery Hall, Queen's University, Kingston Ontario
##     Canada, K7L 3N6


################################################################
##
##  .sphsed
##
##  Phase wrapping routine; takes phases and tracks violation
##  of +/-360 degree boundary, and wraps aliases. For use
##  by demod.dpss().
##
################################################################
.sphsed <-  function(ph,nfreq=length(ph)) {

    q <- 0.0
    pinc <- 0.0

    for(n in 1:nfreq) {
        t1 <- ph[n]
        d <- q - t1
        q <- t1
        if(abs(d) > 180.0) {
            pinc <-  pinc + sign(d)*360.0
        }
        ph[n] <- t1+pinc
    }
    return(ph)
}

################################################################
##  multitaper Development version -- adding stepsize
##  demod.dpss
##
##  Complex demodulation routine. Takes a series x, and 
##  demodulates the series around center frequency centreFreq,
##  using parameters NW, blockLen, and stepSize.
##  Step size has been added here, not added on original version
## Complex demodulation: 

## We begin with a frequency shift, 
## \begin{align}
##  z(t) = x(t) e^{-i 2 \pi f_0 t}. 
## \end{align}
## Then we take the inner product
## \begin{align}
##  y(t) &= \sum_{n=0}^{N-1} v_n^{(0)}(N,W) \, z(t+n) \\
##  &= \sum_{n=0}^{N-1} v_n^{(0)}(N,W) \, x(t+n) \, e^{-i 2 \pi f_0 (t+n)} \\
##  &=  e^{-i 2 \pi f_0 t} \underbrace{\sum_{n=0}^{N-1} v_n^{(0)}(N,W) \, x(t+n) \, e^{-i 2 \pi n f_0}}_{\text{Inner product}}. 
## \end{align}

## We note that the term prior to the under-brace is usually added back on; however, it need not always be.

## the term - 360*deltaT*centreFreq * (1:nResultVals) repesents
##  e^{-i 2 \pi f_0 t} 

## consider some latex in documentation.

##
################################################################


demod.dpss <- function(x, 
                       centreFreq, 
                       nw, 
                       blockLen, 
                       stepSize=1, 
                       wrapphase=TRUE,
                       ...) {

    stopifnot(blockLen %% stepSize == 0)  ## not implemented

    nwTmp <- list(...)$NW
    
    if(!is.null(nwTmp)) {
        warning("NW has been depreciated. Please use nw instead.")
        nw <- nwTmp
    }
    
    ndata <- length(x)

    deltaT <- deltat(x)
    v <- dpss(blockLen, 1, nw)$v
    U0 <- sum(v)

    ampScale <- 2.0/U0
    omegaDeltaT <- 2*pi*centreFreq*deltaT
    jSeq <- (1:blockLen) -1

    complexVal <- exp(-1i*omegaDeltaT*jSeq)
    complexVal <- complexVal*v*ampScale

    nResultVals <- ndata/stepSize - blockLen/stepSize +1

    complexDemod <- complex(nResultVals)

    for(i in 1:nResultVals) {
        ## modified for stepSize
        idx <- 1 + (i-1)*stepSize
        iSeq <- idx:(idx+blockLen-1)
        ## this is the same as t(x[iSeq]) %*%  complexVal
        complexDemod[i] <- crossprod(x[iSeq], complexVal)
    }

    phaseRaw <- Arg(complexDemod)*180/pi
    if(wrapphase) {
        phaseRaw <- multitaper:::.sphsed(phaseRaw)
    }

    ## this is equivlant to adding the following
    ## multitaper:::.sphsed(Arg(exp(-1i*2*pi*centreFreq*(1:nResultVals)))*180/pi)
    ## this term is not added in out dot product, and is generally required,
    ## however it may add a linear trend.

    ## Complex demodulation: 
    
    ## We begin with a frequency shift, 
    ## \begin{align}
    ##  z(t) = x(t) e^{-i 2 \pi f_0 t}. 
    ## \end{align}
    ## Then we take the inner product
    ## \begin{align}
    ##  y(t) &= \sum_{n=0}^{N-1} v_n^{(0)}(N,W) \, z(t+n) \\
    ##  &= \sum_{n=0}^{N-1} v_n^{(0)}(N,W) \, x(t+n) \, e^{-i 2 \pi f_0 (t+n)} \\
    ##  &=  e^{-i 2 \pi f_0 t} \underbrace{\sum_{n=0}^{N-1} v_n^{(0)}(N,W) \, x(t+n) \, e^{-i 2 \pi n f_0}}_{\text{Inner product}}. 
    ## \end{align}

    ## We note that the term prior to the under-brace is usually added back on; however, it need not always be.

    ## the term - 360*deltaT*centreFreq * (1:nResultVals) repesents
    ##  e^{-i 2 \pi f_0 t} 
    
    phase <- phaseRaw - 360*deltaT*centreFreq * (1:nResultVals) 

    ## rule of thumb, use corrected phase, however one does not want
    ## "important signals" obscured by a strong linear trend.
    
    list(amplitude=Mod(complexDemod), phase=phase, complexDemod=complexDemod,
         phaseWithoutCorrection=phaseRaw)
}

## resa <- demod.dpss(CETmonthly$temp, nw=4.5, centreFreq=1/12 ,blockLen=36)
## x <- CETmonthly$temp
## nw=4.5
## stepSize=12
## centreFreq=1/12
## blockLen=36 ## three years


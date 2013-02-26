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



## Utilities for and code for multitaper spectrograms 



## block mtm utils
## are there two uses for mtm, one for blocking the data, and two for looking at the different variates.
## utilities

.nextPowerOf2 <- function(x) {
   return(2^(ceiling(log2(x))));
}

tryAllOffSets <- function(N, blockLength) {
    stopifnot(length(N) == 1, length(blockLength) >= 1)
    resList <- NULL
    nBlocklength <- length(blockLength)
    for(i in 1:nBlocklength) {
        res <- NULL
        for(n in 1:blockLength[i]) {
            nGreaterThanN <-
                getNumberOfBlocks(N, n, blockLength[i])$greaterThanN
            if(nGreaterThanN == 0) {
                res <- c(n, res)
            }
        }
        resList <- c(resList, list(blocklen=blockLength[i], offsets=res))  
        
    }
    return(resList)
}

tryAllBiggerOffSets <- function(N, blockLength) {
    stopifnot(length(N) == 1, length(blockLength) >= 1)
    resList <- NULL
    nBlocklength <- length(blockLength)
    for(i in 1:nBlocklength) {
        res <- NULL
        for(n in blockLength[i]:N) {
            nGreaterThanN <-
                getNumberOfBlocks(N, n, blockLength[i])$greaterThanN
            if(nGreaterThanN == 0) {
                res <- c(n, res)
            }
        }
        resList <- c(resList, list(blocklen=blockLength[i], offsets=res))  
        
    }
    return(resList)
}

## getNumberOfBlocks <-  function(N, blockOffset, blockLength) {
##     num <- 1
##     curOffset <- blockOffset
##     overlap <- 0
##     overlapPercent <- 0
    
##     while(curOffset + blockLength < N) {
##         curOffset <- curOffset + blockOffset
##         num <- num +1
       
##     }
##     if(blockOffset != 0) {
##         num <- num + 1
##         overlapPercent <- (blockLength-blockOffset)/
##             blockLength
##     }
##     nFFT = 4*nextPowerOf2(blockLength)
##     return(list(numberOfBlocks=num,
##                 greaterThanN=curOffset+blockLength-N,
##                 overlapPercent=overlapPercent,
##                 nFFT_4XnextPowerOf2=nFFT))
## }

getNumberOfBlocks <-  function(N, blockOffset, blockLength) {
    idx <- 1:blockLength

    i=0
    currentEnd <- blockLength
    numBlocks <- 1
    while(N > currentEnd) {
        currentEnd <- currentEnd + blockOffset
        ##print(currentEnd)
        numBlocks <- numBlocks +1
        ##print(numBlocks)
    }
    overlapPercent <- (blockLength-blockOffset)/blockLength
    nFFT = 2*.nextPowerOf2(blockLength)
    
    return(list(numberOfBlocks=numBlocks,
                greaterThanN=currentEnd-N,
                overlapPercent=overlapPercent,
                nFFT_2XnextPowerOf2=nFFT))
}

##helper function used in setting up index for plotting spectrogram
meanTimeVectorForBlocks <- function(time1,
                                    nBlockLen,
                                    nBlocks,
                                    nOffset) {
    seqNBlocks <- 1:nBlocks
    seqNBlocksM1 <- seqNBlocks -1
   
    resTimes <- array(NA, nBlocks)

    for(i in seqNBlocksM1) {
        endVal <- (i*nOffset)+nBlockLen
        resTimes[i+1] <- mean(time1[(1+i*nOffset):
                                         endVal])

        ##print(endVal)
    }
    return(resTimes)

}

halfTimeVectorForBlocks <- function(time1,
                                    nBlockLen,
                                    nBlocks,
                                    nOffset) {
    
    halfTimes <- array(NA, nBlocks)
    idx <- 1:nBlockLen

    for( i in ((1:nBlocks) -1)) {
        ##print(i)
        halfTimes[i+1] <- median(time1[idx+nOffset*i])
    }

    return(halfTimes)

}


getIndexClosestVals <- function(v, a) {
    ##sequences cannot have duplicates.
    a <- as.array(a)
    nVals <- length(a)
    n <- length(v)
    seq1 <- 1:n
    res <- array(NA, nVals)

    for(i in 1:nVals) {
        res[i] <- seq1[abs(v-a[i]) == min(abs(v-a[i]))]
    }
    return(res)
}



sampleVarianceBlocks.bias <-  function (x,
                                        nBlockOffset,
                                        nBlockLen,
                                        nBlocks) {
    n <- length(x)
    seqNBlocksM1 <- 1:nBlocks -1
    sigma2 <- array(NA, nBlocks)
    
    for(iBlock in seqNBlocksM1) {
        iBlockP1 <- iBlock +1
        endVal <- iBlock*nBlockOffset + nBlockLen
        if( endVal > n) {
            print("WARNING: last block contains less values than other blocks...")
            ## pad the difference for the last block
            ## with extra zeros...
            len1 <- n - iBlock*nBlockOffset
            temp <- array(0.0, nBlockLen)
            temp[1:len1] <- x[(1+iBlock*nBlockOffset):
                                       n]
            sigma2[iBlockP1] <- .sampleVariance.bias(temp)
        } else {
            sigma2[iBlockP1] <-
                .sampleVariance.bias(x[(1+
                                       iBlock*nBlockOffset):
                                      endVal])         
        }
    }
    return(sigma2)
}

gmean <- function(X) {
    return(2**(mean(log2(X))))
}


spec.mtm.block <- function(timeSeries,
                           k,
                           nw,
                           nFFT,
                           nBlockLen = length(timeSeries),
                           nBlocks = 1,
                           nBlockOffset=nBlockLen, ##nBlockLen is zero overlap
                           deltaT=1.0,
                           dw=NULL, ## tapers$v * sqrt(deltaT)
                           ev=NULL,
                           centreData=TRUE,
                           centreWithSlepians=TRUE,
                           ##if 2, trims 2 extreme points in mean
                           ## used if not centering with slepians
                           numPointsTrimMean=2,
                           returnZeroFreq=T,
                           jackKnife=FALSE,
                           Ftest=FALSE, ## set this to false as it is not working yet
                           jkCIProb=.95,
                           maxAdaptiveIterations=100,
                           sdfTransformation=NULL,
                           returnAuxiliaryData=TRUE,
                           FtestComponents=FALSE){
    timeSeries <- as.double(timeSeries)
    n <- length(timeSeries)
    seqNBlocks <- 1:nBlocks
    seqNBlocksM1 <- seqNBlocks -1
    
    sigma2 <- sampleVarianceBlocks.bias(timeSeries,
                                        nBlockOffset,
                                        nBlockLen,
                                        nBlocks)
    
    receivedDW <- FALSE ## used in returned aux list
    if(is.null(dw)) {
        receivedDW <- TRUE
        ##changed to dpss from dpss dep Aug 2011
        ## this will use djt polarity and lapack
        ## there are problems with dpss_dep ****
        tapers <- dpss(nBlockLen, k, nw=nw, returnEigenvalues=TRUE)##,
                       ##useDJTPolarity=NULL, useLapack=NULL)
        dw <- tapers$v*sqrt(deltaT)
        ev <- tapers$eigen
    }
    
    nFreqs <- nFFT/2 + if(returnZeroFreq) 1 else 0 
    offSet <- if(returnZeroFreq) 0 else 1 
    fiddleFactorSDF <-  deltaT / as.double(nBlockLen) 
    fiddleFactorFreq <- 1 / as.double(nFFT * deltaT)
    
    if(centreData) {
        if(!centreWithSlepians) {
            if(numPointsTrimMean == 0 ) {
                timeSeries <- timeSeries - mean(timeSeries)
            } else if ( numPointsTrimMean > 1 )   {
                timeSeries <- timeSeries -
                    mean(timeSeries, trim=numPointsTrimMean/n)
            }
        } else {
            ## uses fortran 90 with djt polarity and lapack
            ## there are problems with dpss_dep do not use
            ldw <- dpss(n,k, nw=nw, returnEigenvalues=TRUE)##,
                        ##useDJTPolarity=NULL, useLapack=NULL)$v*sqrt(deltaT)
            ldw <- ldw$v
            lswz <- apply(ldw, 2, sum)
            ## zero swz where theoretically zero
            lswz[1:(k/2)*2] <- 0.0
            lssqswz <- sum(lswz**2)
            ## access internal package
            res <- multitaper:::.mweave(timeSeries, ldw, lswz,
                          n, k, lssqswz, deltaT)
            timeSeries <- timeSeries - res$cntr
        }
    }    

##     sdfs <- array(NA, dim=c(nBlocks,nFreqs))
##     dofs <- array(NA, dim=c(nBlocks,nFreqs))
##     wtcfts <- array(NA, dim=c(nBlocks,nFreqs,k))

    sdfs <- array(NA, dim=c(nFreqs,nBlocks))
    dofs <- array(NA, dim=c(nFreqs,nBlocks))
    wtcfts <- array(NA, dim=c(nFreqs,k,nBlocks))
    sas <- array(NA, dim=c(nFreqs,k,nBlocks))
    varjks <- NULL
    bcjks <- NULL
    sjks <- NULL
    UpperCIs <- NULL
    LowerCIs <- NULL
    
    if(jackKnife) {
        print("jackknife in block method needs double checks")
        ##         varjks <- array(NA, dim=c(nBlocks, nFreqs))
        ##         bcjks <- array(NA, dim=c(nBlocks, nFreqs))
        ##         sjks <- array(NA, dim=c(nBlocks, nFreqs))
        ##         UpperCIs <- array(NA, dim=c(nBlocks, nFreqs))
        ##         LowerCIs <- array(NA, dim=c(nBlocks, nFreqs))
        varjks <- array(NA, dim=c(nFreqs,nBlocks))
        bcjks <- array(NA, dim=c(nFreqs,nBlocks))
        sjks <- array(NA, dim=c(nFreqs,nBlocks))
        UpperCIs <- array(NA, dim=c(nFreqs, nBlocks))
        LowerCIs <- array(NA, dim=c(nFreqs, nBlocks))        
    }
    
    FtestResults <- NULL
    FtestResComponents <- NULL
    if(FtestComponents) {
        Ftest <- FALSE
        FtestResComponents <- array(NA, dim=c(nFreqs, 3, nBlocks))
    }
    
    if(Ftest) {
        ##ftestResults <-  array(NA, dim=c(nBlocks, nFreqs))
        FtestResults <-  array(NA, dim=c(nFreqs,nBlocks))
    }
    
    resultFreqs <- ((0+offSet):(nFreqs+offSet-1))*fiddleFactorFreq
    
    for(iBlock in seqNBlocksM1) {
        
        iBlockP1 <-  iBlock +1
        
        endVal <- (iBlock*nBlockOffset)+nBlockLen
        
        taperedData <- NULL
        if( endVal > n) {
            print("WARNING: last block contains less values than other blocks...")
            ## pad the difference for the last block
            ## with extra zeros...
            len1 <- n - iBlock*nBlockOffset
            taperedData <- matrix(0.0, nrow=nBlockLen, ncol=k)
            taperedData[1:len1,] <- dw[1:len1,] *
                timeSeries[(1+iBlock*nBlockOffset):n]
        } else {
            taperedData <- dw*timeSeries[(1+iBlock*nBlockOffset):
                                         endVal]
        }
       
        nPadLen <- nFFT - nBlockLen
        paddedTaperedData <- rbind(taperedData, matrix(0, nPadLen, k))
        cft <- mvfft(paddedTaperedData)
        cft <- cft[(1+offSet):(nFreqs+offSet),]
        sa <- abs(cft)**2
        
        adaptive <-  NULL
        jk <- NULL

        ##need block sigmas  
        if(!jackKnife) {
            adaptive <- multitaper:::.mw2wta(sa, nFreqs, k,
                                             sigma2[iBlockP1], deltaT, ev)
        } else {        
            adaptive <- multitaper:::.mw2jkw(sa, nFreqs, k,
                                             sigma2[iBlockP1], deltaT, ev)
            scl <- exp(qt(jkCIProb,adaptive$dofs)*
                       sqrt(adaptive$varjk))
            ##         varjks[iBlockP1,] <- adaptive$varjk
            ##             bcjks[iBlockP1,] <- adaptive$bcjk
            ##             sjks[iBlockP1,] <- adaptive$sjk
            ##             upperCIs[iBlockP1,] <- adaptive$s*scl
            ##             lowerCIs[iBlockP1,] <- adaptive$s/scl
            varjks[,iBlockP1] <- adaptive$varjk
            bcjks[,iBlockP1] <- adaptive$bcjk
            sjks[,iBlockP1] <- adaptive$sjk
            upperCIs[,iBlockP1] <- adaptive$s*scl
            lowerCIs[,iBlockP1] <- adaptive$s/scl
            
        }

        
        
        swz <- apply(dw, 2, sum)
        ## zero swz where theoretically zero
        swz[1:(k/2)*2] <- 0.0
        ssqswz <- sum(swz**2)
        ## so we average ftests


        Ftest <- FALSE ## this is an issues right now Jan 2013
        ## if(Ftest) {
        ##     ##             if(is.null(swz)) {
        ##     ##                 swz <- apply(dw, 2, sum)
        ##     ##             }
        ##     ## test fortran ftest.
        ##     ##print("using fortran"
        ##     ## optional fortran version of ftest _for not used
        ##     FtestResults[,iBlockP1] <- .HF4mp1_mod(cft,
        ##                                            swz,1,nFreqs,k,##ks=0,
        ##                                            ssqswz=ssqswz)$F_
        ## }

        if(FtestComponents) {
            FtestResComponents[,,iBlockP1] <- 
                .HF4mp1_mod(cft,
                            swz,1,nFreqs,k,
                            ssqswz=ssqswz)
        }
        
        if(is.function(sdfTransformation)) {
            adaptive$s <- sdfTransformation(adaptive$s) 
        }
        
        sdfs[,iBlockP1] <- adaptive$s
        dofs[,iBlockP1] <- adaptive$dofs
        wtcfts[,,iBlockP1] <- cft*sqrt(adaptive$wt)
        sas[,,iBlockP1] <- sa
    }

    if(jackKnife) {
        jk <- list(varjks=varjks,
                   bcjks=bcjks,
                   sjks=sjks,
                   upperCIs=upperCIs,
                   lowerCIs=lowerCIs)
    }
    
    auxiliary <- NULL 
    if(returnAuxiliaryData) {
        if(receivedDW) {
            dw <- NULL
            ev <- NULL
        }
        auxiliary <- list(dw=dw, 
                          ev=ev,
                          wtcfts=wtcfts,
                          sas=sas,
                          nfreqs=nFreqs,
                          nFFT=nFFT,
                          FtestRes=if(Ftest) FtestResults else
                          FtestResComponents) 
    }
    
    minSdf <- apply(sdfs,1,min)
    maxSdf <- apply(sdfs,1,max)
    gmeanSdf <- apply(sdfs,1,gmean)
    
    sdf <- apply(sdfs,1, mean)
    avgLogF <- NULL
    if(Ftest) {
        avgLogF <- apply(log(FtestResults),1,mean)
    }
    
    ##need to average ftests and spectrums....
    ##place individual ftests in aux

    ## Feb 25 2013 change freqs to freq
    ## and aux to mtm to match package code
    return(list(freq=resultFreqs,
                sdf=sdf,
                minSdf=minSdf,
                maxSdf=maxSdf,
                gmeanSdf=gmeanSdf,
                sdfs=sdfs,
                avgLogF=avgLogF,
                dofs=dofs,
                jk=jk,
                mtm=auxiliary))
}


##modified for return ssqres and to loose ks
## this sets ftbase of required for plotting
.HF4mp1_mod<- function(cft,swz,nfrlo,nfrhi,
                      nord ,ssqswz,
                      ftbase=0) { ##ftbase=1.01) {
    ## expect nfrlo=1 and nfrhi=nfft or issues
    ## need to ensure data is in matrix format
    ## for crossprod, tcrossprod.
    swz <- as.matrix(swz) ##hopefully this will do it
    nn <- nord
    fnrm <- as.double(nn-1)
    cmv <- (cft %*% swz) /ssqswz
    ssqave <-  abs(cmv)**2*ssqswz
    ##cftEst <- t(swz %*% t(cmv))
    ##swz <- as.matrix(swz)
    cftEst <- t(tcrossprod(swz,cmv))
        
    ##ssqres <- apply( abs(cft - (cmv %*% t(swz)))**2,
    ##                    1, sum)
    ssqres <- apply( abs(cft - (tcrossprod(cmv,swz)))**2,
                    1, sum)
                
    F_<- fnrm*ssqave/ssqres
    if(ftbase != 0) {
        ftBaseI <- (ftbase > F_)
        F_[ftBaseI] <- ftbase
    }

    ##cols 2 and 3 are complex
    res <- as.matrix(data.frame(F_, cmv, ssqres ))
    colnames(res) <- c("F_", "cmv", "ssqres")
    return(res)
}


#S_{N} is what is used not S_{N-1} in the SAPA Lisp code
.sampleVariance.bias <-function(timeSeries) {
   N <- length(timeSeries);
   return((as.double(t(timeSeries)%*%timeSeries) - sum(timeSeries)^2/N)/N);
}


##no overlaps for ci
bartlettM <- function(sdfs, k, J=dim(sdfs)[2],  nu=2*k) {
    ##nu = k*2
    am <- apply(sdfs, 1, mean)
    gm <- apply(log(sdfs), 1, mean)
    M <- J * nu * (log(am) - gm)
    C <- 1+ (J+1)/(3*J*nu)
    return(list(M=M, C=C, MdivC=M/C))
}

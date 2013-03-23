library(multitaper)

## nord <- k ## number of tapers
## crl <- nw ## time bandwidth
## ndat <- n

## dpss1 <- dpss(ndat, nord, crl )
## vMat <- dpss$v

## Pmat3d <-  array(0, dim=c(nord,nord,nord))

## for( j in 1:nord) {
##     for( k in 1:nord) {
##         for( l in 1:nord) {
##             for( n in 1:ndat) {
##                 Pmat3d[j,k,l]  <- Pmat3d[j,k,l] + vMat[n, j] *  vMat[n, k] *
##                     vMat[n, l]
##             }
##         }
##     }
## }


## default yf

##basicfft <- function(dat, length)

fourierCoef <- function(data, nfft=length(data)) {
    len <- length(data)
    data <- c(data, rep(0, nfft-len))
    fft(data)
}

pgram2 <- function(data, nfft=length(data)) {
    ## we consider the fft to have
    ## [0, f1, f1, .... ,-f2, -f1]
    ## we will need to consider the third frequency for this
    ## not a problem
    originalLength <- length(data)
    fourierCoef1 <- fourierCoef(data, nfft=nfft)
    len <- length(fourierCoef1)
    resLen <- floor(len / 2) + 1
    res <- array(NA, resLen)
    res[1] <- abs(fourierCoef1[1])**2
    for( i in 2:resLen) {
        res[i] <-  fourierCoef1[i] * fourierCoef1[len - i +2]
    }
    res <- res/originalLength
    res
}

## biperiodogram
biPeriodogramNegFreq <- function(data, nfft=length(data)) {
    
    ##for testing
    ##data=1:9
    ##nfft=length(data)
    
    
    originalLength <- length(data)
    fourierCoef1 <- fourierCoef(data, nfft=nfft)
    len <- length(fourierCoef1)
    resLen <- floor(len / 2) + 1
    res <- array(NA, dim=c(len, len))
    
    lenD2 <- as.integer(len/2)
    
    ##seq1 <- NULL
    seq1 <-  c((resLen+1):len, 1:resLen)

    ## freq1 indicates the numerator of frequencies,
    ## actual frequencies are freq1/l
    freq1 <- NULL
    if(len %% 2 == 0 ) { ## even case
    ##seq1 <- c((resLen+1):len, 1:resLen)
        ##fourierCoef1[c((resLen+1):len, 1:resLen)]
        ##lenD2 <- as.integer(len/2)
        freq1 <- -(lenD2-1):lenD2
        ##xSeq <- 1:resLen
        
    } else { ## odd case
        ##seq1 <-  c((resLen+1):len, 1:resLen)
        freq1 <- -(lenD2):lenD2
    }
    
    ## reorder freqs
    fourierCoef1 <-  fourierCoef1[seq1]
    
    for( i in 1:len) {
        for(j in 1:len) {
            freqSum <- freq1[i] + freq1[j]
            thirdIdx <- freqSum %% len
            if(thirdIdx  > lenD2) {
                thirdIdx = len - thirdIdx
            }
            
            res[i,j] <-  fourierCoef1[i] *
                fourierCoef1[j] * fourierCoef1[thirdIdx+lenD2]
        }
        ##res <- res/len
    }
    res <- res/len
    list(bispec=res, freq=freq1)
}


biPeriodogram <- function(data, nfft=length(data)) {
    ## real case only dropping unneeded Fourier coefficients and using conj
    ##for testing
    ##data=1:9
    ##nfft=length(data)
    
    
    originalLength <- length(data)
    fourierCoef1 <- fourierCoef(data, nfft)
    len <- length(fourierCoef1)
    resLen <- as.integer(len / 2) + 1
    fourierCoef1 <- fourierCoef1[1:resLen]
    ## resLen is sufficient
    res <- array(NA, dim=c(resLen, resLen))
    
    lenD2 <- as.integer(len/2)
    
    ##seq1 <- NULL

    ## we only need the top right quadrant
    ## so we need to sort the complex conj
    seq1 <-  c((resLen+1):len, 1:resLen)
    seq1 <- 1:resLen

    ## freq1 indicates the numerator of frequencies,
    ## actual frequencies are freq1/len
    freq1 <- 0:lenD2 
    ## if(len %% 2 == 0 ) { ## even case
    ## ##seq1 <- c((resLen+1):len, 1:resLen)
    ##     ##fourierCoef1[c((resLen+1):len, 1:resLen)]
    ##     ##lenD2 <- as.integer(len/2)
    ##     freq1 <- -(lenD2-1):lenD2
    ##     ##xSeq <- 1:resLen
        
    ## } else { ## odd case
    ##     ##seq1 <-  c((resLen+1):len, 1:resLen)
    ##     freq1 <- -(lenD2):lenD2
    ## }
    
    ## reorder freqs
    ##fourierCoef1 <-  fourierCoef1[seq1]
    
    for( i in 1:resLen) {
        for(j in 1:resLen) {
            ## this time we are only working with postitive
            ## frequencies
            freqSum <- freq1[i] + freq1[j]
            thirdIdx <- freqSum ##%% len
            setConj <- TRUE
            if(thirdIdx  > lenD2) {
                thirdIdx = len - thirdIdx
                setConj <- FALSE
            }

            fCoef12 <- Conj(fourierCoef1[thirdIdx+1])
            if(!setConj) {
                fCoef12 <- fourierCoef1[thirdIdx+1]
            }
            
            res[i,j] <-  fourierCoef1[i] *
                fourierCoef1[j] * fCoef12
        }
        ##res <- res/len
    }
    res <- res/len
    list(bispec=res, freq=freq1)
}

## how to proceede
## see djt's paper
## it appears we have to complex demodulate accross the frquency grid.
## and sum a tripple product like above. 

## we need to do the complex demodulates across the frequency grid and then combine them

## nFFT <- 512
## nFreqs <- as.integer(nFFT/2) +1
## freq1 <- (1:nFreqs) -1
## ## frequency grid for complex demodulates
## freq2 <- freq1/nFFT

## so we need a block length if the block length is the entire series, then we get one value
##demod.dpss(1:10, freq2[1], 4, 2)

data <- 1:10
nord <- 2
crl <- 4


mtm.bispectrum <- function(data, NW, k, dT=1, nFFT=length(data)) {
    ##len <- length(data)

    ## data <- 1:10
    ## NW <- 2
    ## k <-  4
    ## dT <- 1

    ## nFFT <- 32

    
    crl <- NW
    nord <- k
    res <- spec.mtm(data, crl, nord, dT=dT, plot=FALSE,
                    nFFT=nFFT, returnInternals=TRUE)

    vk <- res$mtm$dpss$v
    actualFreqs <- res$freq
    ##res <- spec.mtm(1:10, crl, nord, dT=1, plot=FALSE, returnInternals=TRUE)
    cft <- sqrt(res$mtm$eigenCoefWt) * res$mtm$eigenCoefs
    nFreq <- res$mtm$nfreqs
    nfft <- res$mtm$nFFT
    
    ## the demodMat is djt's \hat(x)(f;t) from eqn 8 in the bispectrum paper
    demodMat <- array(complex(1), dim=c(nFreq, len))
    ##demod each freq
    for(n1 in 1:len) {
        ##sum1 <- t(cft[n1,]) %*% vk[n1,]
        for(f1 in 1:nFreq) {
            demodMat[f1, n1] <- t(cft[f1,]) %*% vk[n1,]
        }
    }
    
    
    

    ##Pmat3d is equation (9) and gamma0 is eqn (12)
    Pmat3d <-  array(0, dim=c(nord,nord,nord))
    gamma0 <- 0
    for( j in 1:nord) {
        for( k in 1:nord) {
            for( l in 1:nord) {
                for( n in 1:len) {
                    Pmat3d[j,k,l]  <- Pmat3d[j,k,l] + vk[n, j] *  vk[n, k] *
                        vk[n, l]
                }
                gamma0 <- gamma0 +  (Pmat3d[j,k,l])^2
            }
        }
    }
    
    ## gamma0 <- 0
    ## for( j in 1:nord) {
    ##     for( k in 1:nord) {
    ##         for( l in 1:nord) {
    ##             gamma0 <- gamma0 + (Pmat3d[j,k,l])^2
                
    ##         }
    ##     }
    ## }
    
    
    ##lenD2 <- as.integer(len/2)
    
    ## freq index
    freq1 <- (1:nFreq-1)
    demodMat1 <- t(demodMat)
    res2 <- array(0i, dim=c(nFreq, nFreq))
    for(i in 1:nFreq) {
        for(j in 1:nFreq) {
            ## this time we are only working with postitive
            ## frequencies
            freqSum <- freq1[i] + freq1[j]
            thirdIdx <- freqSum ##%% len
            setConj <- TRUE
            ## was lenD2
            if(thirdIdx  > (nFreq -1)) {
                thirdIdx = nfft - thirdIdx
                setConj <- FALSE
            }
            
            for(n in 1:len) {
                
                fCoef12 <- Conj(demodMat1[n, thirdIdx+1])
                if(!setConj) {
                    fCoef12 <- demodMat1[n, thirdIdx+1]
                }
                
                res2[i,j] <- res2[i,j] + demodMat1[n, i] *
                    demodMat1[n, j] * fCoef12
            }
        }
        ##res <- res/len
    }
    
    res2 <- res2/gamma0
    list(bispec=res2, freq=actualFreqs)
}

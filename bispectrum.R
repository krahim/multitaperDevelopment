library(multitaper)

nord <- k ## number of tapers
crl <- nw ## time bandwidth
ndat <- n

dpss1 <- dpss(ndat, nord, crl )
vMat <- dpss$v

Pmat3d <-  array(0, dim=c(k,k,k))

for( j in 1:nord) {
    for( k in 1:nord) {
        for( l in 1:nord) {
            for( n in 1:ndat) {
                Pmat3d[j,k,l]  <- Pmat3d[j,k,l] + vMat[n, j] *  vMat[n, k] *
                    vMat[n, l]
            }
        }
    }
}


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
bispectrum <- function(data, nfft=length(data)) {
    
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


## how to proceede
## see djt's paper
## it appears we have to complex demodulate accross the frquency grid.
## and sum a tripple product like above. 

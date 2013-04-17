reshapeMultitaperSpectrumAtFreqIndex <- function(mtmSpec, freqIdx, width=2.5) {
    ## reshape is flatten multitaper spectrum. 
    ##     P and W use 2, djt suggest 2.5
    ## P and W eqn 500 and djt eqn 13.11
    ## P and W H(f) is the same as djt U(f),
    ## in practice these are fourier transforms of the dpss'
    ## in the case of f=0, these are a sum of dpss, this f=0
    ## case is commonly used 

    ## this is only designed for even order fft!!

    stopifnot(mtmSpec$mtm$nFFT %% 2 ==0)
    k <- mtmSpec$mtm$k

    resSpec <- mtmSpec$spec
    halfWidth <- width/2

    ## we use weighted coefficients
    eigenCoefsWt <- sqrt(mtmSpec$mtm$eigenCoefWt) * mtmSpec$mtm$eigenCoefs
    ##eigenCoefsWt <- mtmSpec$mtm$eigenCoefs
    ## pad the dpsses so that he frequencies match
    ## we could use slow ffts or fftw .

    dpss1 <- mtmSpec$mtm$dpss$v
    nDpss <- dim(dpss1)[1]
    ## centre resutant fft by multiplying data by seq 1, -1, 1, -1, ...
    seqM1Centre <- (-1)**(0:(nDpss-1))
    dpss1 <- dpss1 * seqM1Centre

    paddedDPSSCent <- rbind(dpss1,
                        matrix(0, nrow=mtmSpec$pad, ncol=k))

    fftDpssCent <- mvfft(paddedDPSSCent)
    ## using P and W notation 
    ## correction H(f) = deltaT * fft(h_t)
    sqrtDeltaT <- sqrt(deltaT)
    fftDpssCent <- deltaT * fftDpssCent

    nw <- mtmSpec$mtm$nw
    n <- mtmSpec$origin.n
    deltaT <- mtmSpec$mtm$deltaT
    

    ## nFFT is even and centred
    zeroBin <- mtmSpec$mtm$nFFT/2 +1 
    
    cmv1 <- mtmSpec$mtm$cmv[freqIdx]
    
    w <- nw/(n*deltaT)
    ## we expect the zeroth freq but...
    f0 <- mtmSpec$freq[2] - mtmSpec$freq[1]
    
    ## see page 512 P & W
    npts <- as.integer(halfWidth*w/f0)  ## round?? P & W have <= w
    seq1 <- (freqIdx-npts):(freqIdx+npts) ## double check endpoints
    ## we only flatten to zero
    seq1 <- seq1[seq1>0]

    for(j in seq1) {
                
        tempSpec <- 0
        freqDiffIdx <- j - freqIdx
        
        freqDiffIdx <- freqDiffIdx + zeroBin
        
        for(k1 in 1:k) {
            tempSpec <- tempSpec + (abs(eigenCoefsWt[j, k1] -
                                        cmv1 * fftDpssCent[freqDiffIdx, k1] /sqrtDeltaT))**2
        }
        tempSpec <- tempSpec/k
        ##if(j == freqIdx) {
        ##tempSpec <- abs(cmv1)^2 + tempSpec
        ##}
        resSpec[j] <- tempSpec
    }
    resSpec
}

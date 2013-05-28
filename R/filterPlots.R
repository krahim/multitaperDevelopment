plotFilterTime <- function(filt) {
    plot(1:length(filt), filt, ylab="Amplitude", xlab="Samples",
         type="l", xaxs="i", main="Time Domain")
    grid()
}

plotFilterFreq <- function(filt, deltaT=1,
                           nfft=max(1024, 2*2^(ceiling(log2(length(filt)))))) {
    res <- pgram(filt, nfft=nfft, scale=1)
    ## drop last freq and convert to db
    res <- 10*log10(res[-length(res)])
    freq <- (0:(length(res)-1))/(nfft * deltaT)
    plot(freq, res, type="l", xlab="Frequency ", ylab="Magnitude (dB)",
         xaxs="i", main="Frequency Domain")
    grid()
    invisible(list(freq=freq, spec=res))
}
    
plotFilter <- function(filt, deltaT=1,
                       nfft=max(1024, 2*2^(ceiling(log2(length(filt)))))) {
    op <- par(mfrow=c(1,2))
    plotFilterTime(filt)
    res <- plotFilterFreq(filt, deltaT=deltaT, nfft=nfft)
    par(op)
    invisible(res)
}

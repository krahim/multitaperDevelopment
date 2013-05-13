## this gives us essentially the matlab hilbert transform
## and we need the Im to get what we want.

hilbert <- function(dat, n=length(dat)) {
    
    dat <- c(dat, rep(0, n-length(dat)))
    res <- fft(dat)
    h <- rep(0, n)
    nDiv2 <- (n%/%2)
    seq1 <- 2:nDiv2
    h[seq1] <- 2
    h[c(1, nDiv2+1)] <- 1

    res <- res*h
    ## is this n or length(dat) ...
    fft(res, inverse=TRUE)*1/n
}

## plot(Im(hilbert(sin(seq1))))
## lines(-cos(seq1))

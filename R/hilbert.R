## this gives us essentially the matlab hilbert transform
## and we need the Im to get what we want.
## this only works for even values of n only


hilbert <- function(dat, n=length(dat)) {
    stopifnot(n %% 2 == 0)
    ## even values of n only at this time.
    
    dat <- c(dat, rep(0, n-length(dat)))
    res <- fft(dat)
    h <- rep(0, n)
    nDiv2 <- (n%/%2)
    seq1 <- 2:nDiv2
    h[seq1] <- 2
    h[c(1, nDiv2+1)] <- 1

    res <- res*h
    ## is this n or length(dat) ?
    ## no it is not.
    fft(res, inverse=TRUE)*1/n
}

## why /n and not length(dat)
## fft(fft(c(1:10), inverse=TRUE))
## fft(fft(c(1:10,rep(0,6)), inverse=TRUE))/10
## fft(fft(c(1:10,rep(0,6)), inverse=TRUE))/16

## plot(Im(hilbert(sin(seq1))))
## lines(-cos(seq1))

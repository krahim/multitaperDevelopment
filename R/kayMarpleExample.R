## see P and W page 480
## Kay and Marple (1981) page 1386

kayMarpleExample <- function(t=1:16) {
    ## we expect t=1,2,..., N
    res <- 0.9 * cos(2 * pi * t / 7.5) +
        0.9 * cos((2 * pi *t / 5.3) + pi/2) +
        cos(2 * pi * t / 3)
    res
}


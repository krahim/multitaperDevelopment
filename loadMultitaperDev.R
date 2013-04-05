## this provides a quick and dirty way to load the development libraries.
## the path is required
## example of how it is used.
## load.mtm.dev("~/PWLisp/multitaperDevelopment")

## If you are reading this and curious, yes PWlisp refers to the folder I initially downloaded
## and ran the Perciaval and Walden Lisp code.
## The mulititaper package contains modified Fortran code code written by Thomson (1982)
## and at times R code based on the LISP code provided with Percival and Walden (1993).

load.mtm.dev <- function(absolutePathToDir) {

    library("multitaper")
    library("fields")

    source(paste(absolutePathToDir, "bispectrum.R"))
    source(paste(absolutePathToDir, "imagePlot2.R"))
    source(paste(absolutePathToDir, "loeveSpectrum.R"))
    source(paste(absolutePathToDir, "multitaperBlock.R"))
    source(paste(absolutePathToDir, "multitaperHelper2.R"))
    source(paste(absolutePathToDir, "spectraldist.R"))
    source(paste(absolutePathToDir, "quadraticInverse.R"))
    source(paste(absolutePathToDir, "utils2.R"))
}


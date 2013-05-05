## this provides a quick and dirty way to load the development libraries.
## the path is required
## example of how it is used.
## in this case the multitaper development code is in the folder
## ~/PWLisp/multitaperDevelopment/ (on a Linux machine)
## You can use the approprite Windows ar Mac path

## source("~/PWLisp/multitaperDevelopment/loadMultitaperDev.R")
## load.mtm.dev("~/PWLisp/multitaperDevelopment/")

## or you can load the files individually

## If you are reading this and curious, yes PWlisp refers to the folder I initially downloaded
## and ran the Perciaval and Walden Lisp code.
## The mulititaper package contains modified Fortran code code written by Thomson (1982),
## and some of the R code based on the LISP code provided with Percival and Walden (1993).

## load.mtm.dev <- function(absolutePathToDir) {

##     library("multitaper")
##     library("fields")

##     source(paste(absolutePathToDir, "bispectrum.R", sep=""))
##     source(paste(absolutePathToDir, "imagePlot2.R", sep=""))
##     source(paste(absolutePathToDir, "loeveSpectrum.R", sep=""))
##     source(paste(absolutePathToDir, "multitaperBlock.R", sep=""))
##     source(paste(absolutePathToDir, "multitaperHelper2.R", sep=""))
##     source(paste(absolutePathToDir, "spectraldist.R", sep=""))
##     source(paste(absolutePathToDir, "quadraticInverse.R", sep=""))
##     source(paste(absolutePathToDir, "utils2.R", sep=""))
## }

## required to not break my local code don't build ;)
##library("multitaperDevelopment")

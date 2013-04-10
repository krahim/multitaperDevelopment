library("multitaper")


## pushed to multitaper package on github

## the comment had been updated here.
## Please, remember to update this comment regarding the weights in the package
## once verified.
## .mw2wta <- function(sa, nfreq, nord,
##                     var, dt_, ev, evp=(1-ev),
##                     tol=.03, maxadaptiveiteration=100) {

##     ## this is equation (5.3) and (5.4) form 
##     ##   Thomson, D.J. Spectrum Estimation and Harmonic Analysis,
##     ##   Proceedings of the IEEE, 1982.

##     ## note that the weights are squared, they are |d_k(f)^2 from equation
##     ## (5.4)

##     out <- .Fortran("mw2wta", as.double(sa),
##                     wt=matrix(as.double(0), nfreq, nord),
##                     as.integer(nfreq), as.integer(nord),
##                     s=double(nfreq), as.double(ev), as.double(evp),
##                     dofs=double(nfreq), dofav=double(1),
##                     as.double(var), as.double(dt_),
##                     as.double(tol),
##                     as.integer(maxadaptiveiteration),
##                     mxiter=integer(1), aviter=double(1),
##                     PACKAGE='multitaper')
    
##     return(list(s=out$s, wt=out$wt, dofs=out$dofs, dofav=out$dofav,
##                 mxiter=out$mxiter, aviter=out$aviter))
## }

## change in comment has been made in packaged and is on Cran now March 21, 2013

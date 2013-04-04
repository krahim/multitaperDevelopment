##     The multitaper R package (development code)
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2013 Karim Rahim 
##
##     Written by Karim Rahim.
##
##     This file is part of the multitaper package for R.
##     http://cran.r-project.org/web/packages/multitaper/index.html
## 
##     The multitaper package is free software: you can redistribute it and 
##     or modify it under the terms of the GNU General Public License as 
##     published by the Free Software Foundation, either version 2 of the 
##     License, or any later version.
##
##     The multitaper package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
##
##     If you wish to report bugs please contact the author:
## 
##     Karim Rahim
##     karim.rahim@gmail.com
##     Jeffery Hall, Queen's University, Kingston Ontario
##     Canada, K7L 3N6

if(!require("multitaper")){
    stop("multitaper package must be installed")
}

## Code for multitaper non-stationary quadratic inverse. Spectral
## derivitives.





qinssu <- function(L, ndata, nord, crl, deltaT, dw) {
    ## L number of dimensions to keep
    ## ndata data langth ndata
    ## nord number of tapers
    ## crl NW parameter
    ## deltaT
    ## dw vector of dpss multiplied by sqrt(deltaT)
    ## ex dw <- sqrt(deltaT) * dpss(ndata, k=nord, nw=crl)$v
    
    quadraticInverse(L, ndata, nord, crl, deltaT, dw)
}

quadraticInverse <- function (L, ndata, nord, crl, deltaT, dw) {
    ## L number of dimensions to keep
    ## ndata data langth ndata
    ## nord number of tapers
    ## crl NW parameter
    ## deltaT
    ## dw vector of dpss multiplied by sqrt(deltaT)
    ## ex
    ## L <- 2
    ## ndata <- 45
    ## nord <- 6  ## djt is using nord=6 and crel=8 with pinot noir note that (crel > nord)
    ## crl <- 8 ## nw
    ## deltaT <- 1
    ## dw <- sqrt(deltaT) * dpss(ndata, k=nord, nw=crl)$v
    #quadraticInverse(L, ndata, nord, crl, deltaT, dw)
    
    
    ## function
    nTOlim <- round(4*nord)
    NTOM <- L ## input number of dimensions to keep
    cent <- (ndata+1)/2.0
    A1 <- array(0.0, dim=c(ndata,nTOlim))
    F3 <- array(0,0, dim=c(nord,nord,nTOlim))
    FrobP <- array(0.0, dim=c(nTOlim, nTOlim))
    ##Afn <- array(0.0, dim=c(ndata, nTOlim))
    
    ##emprically derived constants.
    ext1 <- 1.045
    ext2 <- 0.045
    
    for( k in 1:nTOlim) {
        span <- k*pi
        ##!c              Extension of about 9 percent to zeroes
        ##! generates and scales  the initial A matrices
        da <- span/(ext1*ndata + ext2*ndata/k)
        
        ##da
        pwr <- 0.0
        for( n in  1:ndata) {
            ang <- da*(n-cent)
            if((k %% 2) == 1) {
                A1[n,k] <- cos(ang)
            }
            else{
                A1[n,k] = sin(ang)
            }
            pwr <- pwr + A1[n,k]**2
        }
        scl = sqrt(ndata/pwr)
        ##for( n in 1:ndata} {
        ##    A1[n,k] = A1[n,k]*scl
        ##}
        A1[,k] <- A1[,k]*scl
    }

    ## need dw = sqrt(eigencoefficients)
    ## we will use internals from a simple set.
    ## keep going, don't give up
    
    ##data(willamette)
    
    ## need dw
    
    ##dpss(ndata, k=nord, nw=crl)
    dw <- sqrt(deltaT) * dpss(ndata, k=nord, nw=crl)$v
    
    ## slepian eigen values near 1 so omitted, just delta t is taken care of.
    
    ## this tripple can be considered V^T A1 V where A1 is an nxn diagonal matrix
    ##and V is nxk
    ## of course we have nTOlim such products.
    
    ##c             Covariance Matrices 
    ## this is eqn 29 from the referenced paper
    ## slepian eigen values near 1 so omitted, just delta t is taken care of.
    ## dw = sqrt(deltaT*V)
    
    ## this tripple can gives us multiple covaraince matrices.
    ## V  nxk but we are using different k in each
    ## A1 is n x NTOlim
    ## if we drop for(  m in 1:nTOlim) then we have V^T V and to give us a k x k matrix
    
    ## so if we matrix out of each column in A, Al where l represents a column.
    ## Al would be matrix of repeated columns, then one F3 matrix would be
    ##  F3 = V^T Al V, now we have nTO lim of these matrices
    ## 
    for( j in 1:nord) {
        for( k in 1:nord) {
            for(  m in 1:nTOlim) {
                F3[k,j,m] = 0.0
                for( n in 1:ndata) {
                    F3[k,j,m] <- F3[k,j,m] + dw[n,j]*dw[n,k]*A1[n,m]
                }
                F3[k,j,m] = F3[k,j,m]/deltaT
            }
        }
    }
    
    
    ## construction with three dimensions.
    ##!  c   3d Frobenius Products of Covariance Matrices
    ##! http://en.wikipedia.org/wiki/Matrix_multiplication#Frobenius_product
    
    ## we vectorize the first two dimensions of F3 and then do an inner product
    ## F3 is kxkxnTOlim  and F3A is k^2 x nTOlim and we get F3A^T F3A,
    ## the loops are okay for this
    
    ## this is also the sum of the Hadamard, Schur, or entrywise product.
    for( m in 1:nTOlim) {
        for( n in 1:nTOlim) {
            FrobP[n,m] <- 0.0
            for(  j in 1:nord) {
                for( k in 1:nord) {
                ##!c   400 FrobP(n,m) = sdot(KbK,F3(1,1,n),1,F3(1,1,m),1)
                    FrobP[n,m] <- FrobP[n,m] + F3[j,k,n] * F3[j,k,m]
                }
            }
        }
    }
    
    
    eigenRes <- eigen(FrobP)
    FPeval <- eigenRes$values
    FPevec <- eigenRes$vectors
    
    ## eigenvalues of FrobP give you equation 18 in the spie paper.
    
    ## there is an error check for the number greater than zero
    
    nGnz <- 0 ## number greater than zero
    for( k in 1:nTOlim) {
        if(FPeval[k] > 1*10^(-3)) {
            nGnz <- k
        }
    }

    ##sprintf("%s %d.", "num evals greater than zero", nGnz)
    
    
    ## Afn is ndata x NTOM
    
    ##c     Orthonogal basis Functions,
    ## this is a stanard matrix multiplication of
    ## A1 by FPevec dropping columns in the latter
    ## after the matrix multiplication a standardization is implemented
    
    ##Afn <- array(0.0, dim=c(ndata, nTOlim))
    nsc <- ndata/2
    
    ## used in slope normalization.
    npc <- as.integer((ndata+3)/2)
    
    
    Afn <- A1 %*% FPevec[,1:2]
    
    for( m in 1:NTOM) { ## this drops columns in FPevec
        ## do Afn <- A1 %*% FPevec[,1:2]
        ## instead ...
        
        ##for( n in 1:ndata) {
        ##    Afn[n,m] <- 0.0
        ##    for(j in 1:nTOlim) {
        ##        ## messy
        ##        Afn[n,m] = Afn[n,m] + A1[n,j]*FPevec[j,m]
        ##    }
        ##}
        
        ## Afn[,m] is a one column vector
        ##fnrm <-  sqrt(ndata/(Afn[,m] %*% Afn[,m]))
        
        Afn <- A1 %*% FPevec[,1:2]
        
        fnrm <- sqrt(ndata/crossprod(Afn[,m]))
        
        if(Afn[npc,m] < 0.0) {
            fnrm <- -fnrm
        }
        Afn[,m] <- Afn[,m] * fnrm
    }
    
    
    
    ## some declarations required
    
    ##c              Generate Trace-Orthogonal Matrices
    ##    write(6,*) " nTOlim,nord,NTOM ",nTOlim,nord,nGnz

    Atr <- array(0, dim=c(NTOM,1))
    Gev <- FPeval[1:NTOM]
    ## complex array, imaginary values will still be set to zero.
    A3m <- array(0i, dim=c(nord, nord, NTOM))
    Fnorm <- array(0, dim=c(NTOM, 1))
    xBcf <- array(0, dim=c(NTOM, 1))
    
    ## this is the same tripple product reversed, except we are now using dw
    
    ## this is why we use the eigen function   Nowwer are doing the same tripple product
    ## as above but using AFN
    for( l in 1:NTOM) {
        ##Atr[l] = 0..
        ##Gev[l] = FPeval[l]
        for( j in 1:nord) {
            for( k in 1:nord) {
                F3[j,k,l] = 0.0
                for( n in 1:ndata)  {
                    F3[j,k,l] <- F3[j,k,l] + dw[n,j]*dw[n,k]*Afn[n,l]
                }
                F3[j,k,l] <- F3[j,k,l]/deltaT
                A3m[j,k,l] <- complex(real= F3[j,k,l], imaginary= 0.0)
            }
            Atr[l] = Atr[l] + F3[j,j,l]
        }
        ## how is the matrix stored?? likely the correct order for this
        Fnorm[l] = as.vector(F3[,,l]) %*% as.vector(F3[,,l])
        xBcf[l] <- sum(Afn[,l])/ndata
    }
    
    ## Afn is the function...
    ## 
    ## this should be the end of the function.
    ##
    ## we may print out data and lea
    
    
    ##return
    
    list(Gev=Gev, Fnorm=Fnorm, Trace=Atr, A.bar=xBcf, A3m=A3m, Afn=Afn)
    
    
}

## inputs
## L <- 2
## ndata <- 45
## nord <- 6  ## djt is using nord=6 and crel=8 with pinot noir note that (crel > nord)
## crl <- 8 ## nw
## deltaT <- 1
## dw <- sqrt(deltaT) * dpss(ndata, k=nord, nw=crl)$v

## quadraticInverse(L, ndata, nord, crl, deltaT, dw)

## ##This generates the AI



  ## call eenter('qicoef',ncseqn)
  ##     call djta(numarg(),18)
  ##     if(ndecs.lt.1) call djtf('ndecs < 1 ')
  ##     if(nfslo*ndecs.lt.nfrlo) call djtf('nfslo*ndecs < nfrlo ')
  ##     if(nfshi*ndecs.gt.nfrhi) call djtf('nfshi*ndecs > nfrhi ')
  ##     tsum = ssum(L,xBcf,1)
  ##     trnrm = sdot(L,xBcf,1,xBcf,1)

## we do one section at time
## djt goes over the different sections here
## cft dimensins nfreqs, nord, nblocks for me
## we will set have cft <- cft[,,sect] 
## B3m is A3m 


## A3m and B3m are the same, xBcf is Z.bar

mwmbsd <- function(L, cft, A3m, xBcf, Gev, curblk, nord, curs=TRUE) {
    getQICcoefficients1Block(L, cft, A3m, xBcf, Gev, curblk, nord, curs)
}


## this is very slow right now, and will need a little tweaking and thought
##
getQICcoefficients <- function(L, cfts, A3m, xBcf, Gev, nord, curs=TRUE) {
    nFreqs <- dim(cfts)[1]
    nBlocks <- dim(cfts)[3]
    QIC <- array(0, dim=c(nFreqs, nBlocks, L))
    urs <- array(0, dim=c(nFreqs, nBlocks))
    for(j in 1:nBlocks) {
        res <- getQICcoefficients1Block(L, cfts[,,j], A3m, xBcf, Gev,
                                        curblk, nord, curs=TRUE)
        QIC[,j,] <- res$QIC
        urs[,j] <- res$urs        
    }
    list(QIC=QIC, urs=urs)
        
}

    
getQICcoefficients1Block <- function(L, cft, A3m, xBcf, Gev, curblk, nord, curs=TRUE) {
    tsum <- sum(xBcf)
    trnrm <- crossprod(xBcf) ##t(xBcf) %*% xBcf
    fnord <- nord
    nFreqs <- dim(cft)[1]

    sb <- 0.0
    sbar <- array(0, dim=c(nFreqs))
    csm <- array(0i, dim=c(L))
    QIC <- array(0, dim=c(nFreqs, L))
    urs <- array(0, dim=c(nFreqs))
    
    for( nn in 1:nFreqs) { ##1111 david tosses away lower freqs
        ##sb <- 0.0
        ##sbar[nn] <- 0.0
        n <- nn ## I don't decimate at this point
        for( k in 1:nord) { ## 100
            sbar[nn] = sbar[nn] + Re(cft[n,k])**2 + Im(cft[n,k])**2
        } ##100
        sbar[nn] <- sbar[nn]/fnord
        for( m  in  1:L) { ## 777
            ##csm[m] = 0i
            for(j in 1:nord) { ## 480
                for( k in 1:nord) { ## 480
                    csm[m] <- csm[m] + Conj(cft[n,j])*A3m[k,j,m]*cft[n,k]
                }
            }
            QIC[nn,m] <- Re(csm[m])/Gev[m]
            sb <- sb + QIC[nn,m]*xBcf[m]
        } ##777 continue
        if(curs)  {
            sb <- sb/trnrm
            ##urs[nn] = 0.0
            for( m in 1:L) { ## 850
                urs[nn] <- urs[nn] + Gev[m]*((QIC[nn,m]/sb - xBcf[m])**2)
            }
        }
    } ##1111 continue

    list(QIC=QIC, urs=urs)
    ##  call eleave(ncseqn)
}



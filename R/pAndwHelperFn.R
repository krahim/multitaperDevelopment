##Helper function's used in P and W (1993) lisp source code 
################################################################################
##library(bitops)

getNDFT <- function(nNonZeroFreq, sampleSize) {
   if(is.character(nNonZeroFreq)) {
      return (switch(tolower(nNonZeroFreq),
                  "halfnextpowerof2" =  nextPowerOf2(sampleSize),
                  "nextpowerof2" =  2*nextPowerOf2(sampleSize),
                  "twicenextpowerof2" =  4*nextPowerOf2(sampleSize),
                  "fourier" = sampleSize));
   }
   if(is.numeric(nNonZeroFreq) &&
                  powerOf2(nNonZeroFreq) &&
                  2*nNonZeroFreq >= sampleSize) {
      return(2*nNonZeroFreq);
   }
   return();
}

getNFreqs <- function(nNonZeroFreq, sampleSize, returnEstFor0FreqP) {
   if(!is.logical(returnEstFor0FreqP)) {
      return();
   }
   if(returnEstFor0FreqP) {
      return(1+ getNDFT(nNonZeroFreq, sampleSize) %/% 2);
   }
   return(getNDFT(nNonZeroFreq, sampleSize) %/% 2);
}

convertTodB <- function(X) {
   if(!identical(X[X<=0], numeric(0))) {
      return();
   }
   return(10*log10(X));
}

convertFromdB <- function(X) {
   return(10^(X/10));
}

#library(bitops)

powerOf2 <- function(x) {
   if(as.integer(x)!=x) {
      return();
   }
   if(bitAnd(x,x-1) == 0) {
      return(log2(x));
   }
   else {
      return();
   }
}

nextPowerOf2 <- function(x) {
   return(2^(ceiling(log2(x))));
}

NyquistFrequency  <- function(sampleTime) {
   return(1/(2*sampleTime));
}

sampleTime <- function(NyquistFrequency) {
   return(1/(2*NyquistFrequency));
}

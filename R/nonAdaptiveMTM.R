blockMTMevWeights <- function(blockMTMRes) {
    ev <- blockMTMRes$mtm$ev
    sumEv <- sum(ev)
    nBlocks <- dim(blockMTMRes$mtm$sas)[3]
    nFreq <-  dim(blockMTMRes$mtm$sas)[1]
    sas <- blockMTMRes$mtm$sas
    specRes <- array(NA, dim=c(nFreq, nBlocks))
    for(i in 1:nBlocks) {
        specRes[,i] <- (blockMTMRes$mtm$sas[,,i] %*% ev)/sumEv
    }
    specRes
}

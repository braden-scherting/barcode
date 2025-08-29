library(Rcpp)
suppressPackageStartupMessages(library(tidyverse))
library(sp)
library(stringr)
library(doFuture, quietly = TRUE)
library(future.apply)
library(poismf)

loadData <- function(noFactors=6, initialization=F){
  # Loads and packages data for input to sampler

  dat <- readRDS("../data/birds2.rds"); names(dat)[2] <- "X"

  ij <- which(dat$Y!=0, arr.ind = T)
  dat$RspY <- Matrix::sparseMatrix(i=ij[,1], j=ij[,2], x=dat$Y[dat$Y!=0])


  dims <- list(N=nrow(dat$X), P=ncol(dat$Y), L=noFactors, Q=ncol(dat$X), K=ncol(dat$des01))

  Sigs <- array(dim=c(dims$K, dims$K, dims$L))
  Rinv <- array(dim=c(dims$K, dims$K, dims$L))
  condVar <- array(dim=c(dims$K, dims$K, dims$L))
  condVarChol <- array(dim=c(dims$K, dims$K, dims$L))

  for (spat in 1:dims$L){
    Sigs[,,spat] <- 1 * (exp(-dat$D / 0.5) + diag(sqrt(.Machine$double.eps), dims$K))
    Rinv[,,spat] <- solve(Sigs[,,spat])
    condVar[,,spat] <- solve(Rinv[,,spat] + t(dat$des01) %*% dat$des01)
    condVarChol[,,spat] <- t(chol(condVar[,,spat]))
  }

  Ypresent <- !is.na(dat$Y)

  storeInit <- list(condVar=condVar, condVarChol=condVarChol, Ypresent=Ypresent)

  blockSize <- min(dims$L, 3)
  storeInit$Cstar <- matrix(0, ncol=blockSize, nrow=2^blockSize)
  for (i in 1:(2^blockSize - 1)){
    storeInit$Cstar[i+1,] <- binaryLogic::as.binary(i, n=blockSize)
  }


  if (initialization) {
    mf <- poismf(dat$Y, k=dims$L - 1, l2_reg=0, niter=100, nthreads=1)
    storeInit$Ginit <- cbind(rep(1, dims$P), t(mf$B) * (rep(1, dims$P) %*% t(rowSums(mf$A))))
    storeInit$Pinit <- cbind(rep(1/dims$N, dims$N), t(mf$A) / (rep(1, dims$N) %*% t(rowSums(mf$A))))
  } else {
    storeInit$Pinit <- matrix(1/dims$N, nrow=dims$N, ncol=dims$L)
    storeInit$Ginit <- matrix(rgamma(dims$P*dims$L, 2,1), nrow=dims$P, ncol=dims$L)
  }

  return(list(dims=dims, dat=dat, storeInit=storeInit, coords=dat$long.lat))
}

runOneChain <- function(dims, dat, storeInit, workerSeed=sample(10:1e7,1),
                        workerNo=0, niter=200, nthin=1, nburn=200, pred=0, blockSize=5){
  sourceCpp("../src/barcodeSampler.cpp")

  # blockSize <- min(dims$L, blockSize)
  # storeInit$Cstar <- matrix(0, ncol=blockSize, nrow=2^blockSize)
  # for (i in 1:(2^blockSize - 1)){
  #   storeInit$Cstar[i+1,] <- binaryLogic::as.binary(i, n=blockSize)
  # }

  return(sampler_camellia(dat, dims, storeInit, workerSeed, niter, nthin, nburn,
                          x=workerNo))
}

runParallel <- function(inputs, nParallelChains=4, seedSeed=172024, prefix="out", ...){

  registerDoFuture()
  plan(multisession, workers = nParallelChains)

  predTested <- ifelse(any(is.na(inputs$dat$Y)) && pred==TRUE, 1, 0)

  tmpFit <- foreach(i = 1:4, .options.future = list(seed = TRUE)) %dofuture% {
    runOneChain(dims=inputs$dims, dat=inputs$dat, storeInit=inputs$storeInit,
                workerSeed=round(seedSeed * i / (i+1)), workerNo=i, pred=predTested,
                ...)
  }

  saveRDS(tmpFit, paste0("../output/results/", prefix, format(Sys.time(), "%H%M%a%d%b%Y"), ".rds"))
  return(NULL)
}

stackChains <- function(list4, name){
  # Collates samples when multiple chains are run (1/2)

  dims <- dim(list4[[1]][[name]])
  ndim <- length(dims[dims!=1])


  niter <- dims[ndim]
  dims[ndim] <- dims[ndim]*4

  tmpArr <- array(dim=dims)

  if (ndim==3){
    for (i in 1:4){
      tmpArr[,,((i-1)*niter + 1):(i*niter)] <- list4[[i]][[name]]
    }
  } else if (ndim==2){
    for (i in 1:4){
      tmpArr[,((i-1)*niter + 1):(i*niter)] <- list4[[i]][[name]]
    }
  } else {
    for (i in 1:4){
      tmpArr[((i-1)*niter + 1):(i*niter)] <- list4[[i]][[name]]
    }
  }
  return(tmpArr)
}

stackedList <- function(list4, pars=c("Gamma","Phi", "C", "S", "Beta", "Xi", "bGamma", "psi")){
  # Collates samples when multiple chains are run (2/2)
  newList <- list()

  for (name in pars){
    newList[[name]] <- stackChains(list4, name)
  }
  return(newList)
}


# new comment

library(Rcpp)
suppressPackageStartupMessages(library(tidyverse))
library(sp)
library(stringr)
library(NMF)
library(doFuture, quietly = TRUE)
library(future.apply)

loadData <- function(nYears=11, nSpecies=137, noFactors=10, initialization=T, yearSetFlag=F, yearSet=F){
  load("data/unfitted_models.RData")
  
  # number of years to consider
  minusYear <- -nYears
  
  if (yearSetFlag){
    obsTF <-  -models$abu$XData$year %in% yearSet
  } else {
    obsTF <- models$abu$XData$year > minusYear
  }
  
  XData <- models$abu$X[obsTF,] |> as.data.frame() 
  Y <- models$abu$Y[rownames(models$abu$XData[obsTF,]), ]
  
  # number of species to consider
  P <- min(c(nSpecies, ncol(Y)))
  Y <-  Y[,1:P]
  
  nonzeroRows <- which(rowSums(Y)>0)
  
  XData <- XData[nonzeroRows, ]
  Y <- Y[nonzeroRows, ]
  
  colnames(XData) <- plyr::mapvalues(colnames(XData), 
                                     colnames(XData)[startsWith(colnames(XData), "poly")],
                                     c("polyAprMay1", "polyAprMay2",
                                       "polyDecFeb1", "polyDecFeb2",
                                       "polyJunJul1", "polyJunJul2"))
  
  XFormula <- models$abu$XFormula
  XFormula <- update.formula(XFormula, 
                             new = "~.-poly(AprMay, degree = 2, raw = TRUE) -
                                        poly(DecFeb, degree = 2, raw = TRUE) -
                                        poly(JunJul, degree = 2, raw = TRUE) -
                                        VMI_PC2 - VMI_PC3 - VMI_PC4 - VMI_PC5 - 
                                        year -1 +
                                        # as.factor(year) +
                                        polyAprMay1 +
                                        polyDecFeb1 +
                                        polyJunJul1")
  
  XData <- model.matrix(XFormula, XData)
  XData <- scale(XData)
  XData <- cbind(rep(1, nrow(XData)), XData)
  colnames(XData)[1] <- "(Intercept)"
  
  recentDesign <- models$abu$studyDesign[obsTF,]
  recentDesign <- recentDesign[nonzeroRows,]
  studyDesign <- as.character(recentDesign$route)
  routes <- unique(studyDesign)
  
  N <- nrow(XData)
  Q <- ncol(XData)
  K <- length(routes)
  
  des <- as.numeric(plyr::mapvalues(studyDesign, routes, 1:length(routes))) 
  des01 <- matrix(0, nrow=N, ncol=K)
  des01[cbind(1:N, des)] <- 1
  
  long.lat <- cbind(as.data.frame(models$abu$rL$route$s)[routes,]$long, 
                    as.data.frame(models$abu$rL$route$s)[routes,]$lat)
  
  D <- plgp::distance(long.lat)
  
  dat <- list(X=as.matrix(XData), Y =Y, des01=des01, D=D, XFormula=XFormula)
  dims <- list(N=N, P=P, L=noFactors, Q=Q, K=K)
  
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
  
  if (initialization){
    nmfInit <- NMF::nmf(dat$Y, noFactors, nrun=10, seed=441503)
    PhiInit <- NMF::basis(nmfInit)
    GammaInit <- t(NMF::coef(nmfInit))
    # ab <- sqrt(mean(PhiInit) / mean(GammaInit))
    # PhiInit <- as.matrix(PhiInit / ab)
    # GammaInit <- t(as.matrix(GammaInit * ab))
    GammaInit <- GammaInit * colSums(PhiInit)
    PhiInit <- PhiInit / colSums(PhiInit)
    
    storeInit$PhiInit <- PhiInit
    storeInit$GammaInit <- GammaInit
  } else {
    storeInit$PhiInit <- matrix(1/dims$N, nrow=dims$N, ncol=dims$L)
    storeInit$GammaInit <- matrix(rgamma(dims$P*dims$L, 2,1), nrow=dims$P, ncol=dims$L)
  }
  
  return(list(dims=dims, dat=dat, storeInit=storeInit, coords=long.lat))
}

runOneChain <- function(dims, dat, storeInit, workerSeed=sample(10:1e7,1), 
                        workerNo=0, niter=200, nthin=1, nburn=200, pred=0, blockSize=5){
  sourceCpp("models/sampler_snapdragon.cpp")
  
  blockSize <- min(dims$L, blockSize)
  storeInit$Cstar <- matrix(0, ncol=blockSize, nrow=2^blockSize)
  for (i in 1:(2^blockSize - 1)){
    storeInit$Cstar[i+1,] <- binaryLogic::as.binary(i, n=blockSize)
  }
  
  return(sampler_snapdragon(niter, nthin, nburn, dims, dat, storeInit, 
                          workerSeed, x=workerNo, pred=pred))
}

runParallel <- function(inputs, nParallelChains=4, seedSeed=172024, ...){
  
  registerDoFuture()
  plan(multisession, workers = nParallelChains)
  
  predTested <- ifelse(any(is.na(inputs$dat$Y)) && pred==TRUE, 1, 0)
  
  tmpFit <- foreach(i = 1:4, .options.future = list(seed = TRUE)) %dofuture% {
    runOneChain(dims=inputs$dims, dat=inputs$dat, storeInit=inputs$storeInit, 
                workerSeed=round(seedSeed * i / (i+1)), workerNo=i, pred=predTested, 
                ...)
  } 
  
  saveRDS(tmpFit, paste0("results/", format(Sys.time(), "%H%M%a%d%b%Y"), ".rds"))
  return(NULL)
}

stackChains <- function(list4, name){
  dims <- dim(list4[[1]][[name]])
  ndim <- length(dims)
  
  
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

stackedList <- function(list4){
  newList <- list()
  for (name in names(list4[[1]])){
    newList[[name]] <- stackChains(list4, name)
  }
  return(newList)
}


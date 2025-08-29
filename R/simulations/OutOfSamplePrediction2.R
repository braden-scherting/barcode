library(gllvm)
source("setup_functions.R")
sourceCpp("../src/barcodeSampler.cpp")

getCVdata <- function(dat, dims){
  ij <- which(dat$Y>0 & !is.na(tmp$dat$Y), arr.ind = T)
  dat$RspY <- Matrix::sparseMatrix(i=ij[,1], j=ij[,2], x=dat$Y[dat$Y>0 & !is.na(tmp$dat$Y)])

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
  Ymissing <- which(is.na(dat$Y), arr.ind = T)
  
  storeInit <- list(condVar=condVar, condVarChol=condVarChol, Ypresent=Ypresent,
                    Ymissing=Ymissing)
  
  blockSize <- min(dims$L, 3)
  storeInit$Cstar <- matrix(0, ncol=blockSize, nrow=2^blockSize)
  for (i in 1:(2^blockSize - 1)){
    storeInit$Cstar[i+1,] <- binaryLogic::as.binary(i, n=blockSize)
  }
  
  mf <- poismf(dat$Y, k=dims$L - 1, l2_reg=0, niter=100, nthreads = 1)
  storeInit$Ginit <- cbind(rep(1, dims$P), t(mf$B) * (rep(1, dims$P) %*% t(rowSums(mf$A))))
  storeInit$Pinit <- cbind(rep(1/dims$N, dims$N), t(mf$A) / (rep(1, dims$N) %*% t(rowSums(mf$A))))
  
  return(list(dims=dims, dat=dat, storeInit=storeInit, coords=dat$long.lat))
}


set.seed(331825)
Full <- loadData()

#----Create Folds-----#
N <- Full$dims$N; P <- Full$dims$P; L <- 5+1
allInds <- sample(1:(P*N), P*N, F)
#---------------------#

des <- apply(Full$dat$des01, 1, function(x) which(x==1))


nFolds <- 3

preds <- list()

foldRMSE <- matrix(NA, nFolds, 5)

for (fold in 2:nFolds){
  fInds <- sort(allInds[(((fold-1)*(N*P/nFolds) + 1):((fold)*(N*P/nFolds)))])
  Ycv <- as.matrix(Full$dat$Y); Ycv[fInds] <- NA
  
  outPois <- gllvm(Ycv, Full$dat$X, formula=Full$dat$XFormula, family = "poisson",
                     studyDesign=des,
                     dist = Full$coords, starting.val="res", seed = 6196041 %*% fold)
    
  print(paste("Pois: ", fold, " done"))
  
  outNB <- gllvm(Ycv, Full$dat$X, formula=Full$dat$XFormula, family = "negative.binomial", studyDesign=des,
                   dist = Full$coords, starting.val="res", seed = 6196041 %*% fold)
  print(paste("NB: ", fold, " done"))
  #--------------------------------------#
  
  tmp <- loadData(noFactors = 4); tmp$dat$Y <- Ycv
  ls <- getCVdata(tmp$dat, tmp$dims)
  niter=2500; nburn=2500
  outBarcode4 <- sampler_camellia_missing(ls$dat, ls$dims, ls$storeInit, 6196043 %*% fold, niter, 1, nburn, 1)
  
  Ebarcode4 <- numeric(length(fInds))
  
  for (it in 1:niter){
    Ebarcode4 <- Ebarcode4 + ((outBarcode4$Phi[,,it] * outBarcode4$C[,,it]) %*% t(outBarcode4$Gamma[,,it] * outBarcode4$S[,,it]))[(fInds)]
  }
  Ebarcode4 <- Ebarcode4 / niter
  
  p04 <- numeric(length(fInds))
  
  for (it in 1:niter){
    p04 <- p04 + dpois(0,((outBarcode4$Phi[,,it] * outBarcode4$C[,,it]) %*% t(outBarcode4$Gamma[,,it] * outBarcode4$S[,,it]))[fInds])
  }
  p04 <- p04 / niter
  
  
  print(paste("bar4: ", fold, " done"))
  
  #--------------------------------------#
  
  tmp <- loadData(noFactors = 7); tmp$dat$Y <- Ycv
  ls <- getCVdata(tmp$dat, tmp$dims)
  niter=2500; nburn=2500
  outBarcode7 <- sampler_camellia_missing(ls$dat, ls$dims, ls$storeInit, 6196040 %*% fold, niter, 1, nburn, 1)
  
  Ebarcode7 <- numeric(length(fInds))
  
  for (it in 1:niter){
    Ebarcode7 <- Ebarcode7 + ((outBarcode7$Phi[,,it] * outBarcode7$C[,,it]) %*% t(outBarcode7$Gamma[,,it] * outBarcode7$S[,,it]))[(fInds)]
  }
  Ebarcode7 <- Ebarcode7 / niter
  
  p07 <- numeric(length(fInds))
  
  for (it in 1:niter){
    p07 <- p07 + dpois(0,((outBarcode7$Phi[,,it] * outBarcode7$C[,,it]) %*% t(outBarcode7$Gamma[,,it] * outBarcode7$S[,,it]))[fInds])
  }
  p07 <- p07 / niter
  
  
  print(paste("bar7: ", fold, " done"))
  
  #--------------------------------------#
  
  tmp <- loadData(noFactors = 10); tmp$dat$Y <- Ycv
  ls2 <- getCVdata(tmp$dat, tmp$dims)
  niter=2500; nburn=2500
  outBarcode10 <- sampler_camellia_missing(ls2$dat, ls2$dims, ls2$storeInit, 6196044 %*% fold, niter, 1, nburn, 1)
  Ebarcode10 <- numeric(length(fInds))
  for (it in 1:niter){
    Ebarcode10 <- Ebarcode10 + ((outBarcode10$Phi[,,it] * outBarcode10$C[,,it]) %*% 
                                  t(outBarcode10$Gamma[,,it] * outBarcode10$S[,,it]))[fInds]
  }
  Ebarcode10 <- Ebarcode10 / niter
  
  p010 <- numeric(length(fInds))
  
  for (it in 1:niter){
    p010 <- p010 + dpois(0,((outBarcode10$Phi[,,it] * outBarcode10$C[,,it]) %*% t(outBarcode10$Gamma[,,it] * outBarcode10$S[,,it]))[fInds])
  }
  p010 <- p010 / niter
  
  print(paste("bar10: ", fold, " done"))
  #--------------------------------------#
  
  foldRMSE[fold,1] <- sqrt(mean((Full$dat$Y[is.na(Ycv)] - Ebarcode4)^2))
  foldRMSE[fold,2] <- sqrt(mean((Full$dat$Y[is.na(Ycv)] - Ebarcode7)^2))
  foldRMSE[fold,3] <- sqrt(mean((Full$dat$Y[is.na(Ycv)] - Ebarcode10)^2))
  foldRMSE[fold,4] <- sqrt(mean((Full$dat$Y[is.na(Ycv)] - exp(predict(outPois))[is.na(Ycv)])^2))
  foldRMSE[fold,5] <- sqrt(mean((Full$dat$Y[is.na(Ycv)] - exp(predict(outNB))[is.na(Ycv)])^2))
  
  preds[[fold]] <- array(cbind(
    Full$dat$Y[is.na(Ycv)],
    Ebarcode4,
    Ebarcode7,
    Ebarcode10,
    exp(predict(outPois))[is.na(Ycv)],
    exp(predict(outNB))[is.na(Ycv)],
    p04,
    p07,
    p010,
    dpois(0, exp(predict(outPois))[is.na(Ycv)]),
    dpois(0, exp(predict(outNB))[is.na(Ycv)])
  ), dim=c(length(fInds), 11))
}


saveRDS(foldRMSE, "../output/results/foldRMSE.rds")
saveRDS(foldRMSE, "../output/results/OOSPEpreds.rds")


c(sqrt(mean((preds[[1]][preds[[1]][,5]<1890,5] - preds[[1]][preds[[1]][,5]<1890,1])^2)),
  sqrt(mean((preds[[2]][preds[[2]][,5]<1890,5] - preds[[2]][preds[[2]][,5]<1890,1])^2)),
  sqrt(mean((preds[[3]][preds[[3]][,5]<1890,5] - preds[[3]][preds[[3]][,5]<1890,1])^2)) ) |> mean()

c(sqrt(mean((preds[[1]][preds[[1]][,6]<1890,6] - preds[[1]][preds[[1]][,6]<1890,1])^2)),
  sqrt(mean((preds[[2]][preds[[2]][,6]<1890,6] - preds[[2]][preds[[2]][,6]<1890,1])^2)),
  sqrt(mean((preds[[3]][preds[[3]][,6]<1890,6] - preds[[3]][preds[[3]][,6]<1890,1])^2)) ) |> mean()
library(NMF)
library(Rcpp)
library(tidyverse)
library(latex2exp)
library(pROC)

source("SetupFunctions.R")
sourceCpp("./models/sampler_snapdragonXonly.cpp")

getStoreInit <- function(dat, dims){
  storeInit <- list(Ypresent=1 * !is.na(dat$Y))
  
  blockSize <- min(dims$L, 5)
  storeInit$Cstar <- matrix(0, ncol=blockSize, nrow=2^blockSize)
  for (i in 1:(2^blockSize - 1)){
    storeInit$Cstar[i+1,] <- binaryLogic::as.binary(i, n=blockSize)
  }
  
  storeInit$PhiInit <- matrix(1/dims$N, nrow=dims$N, ncol=dims$L)
  storeInit$GammaInit <- matrix(rgamma(dims$P*dims$L, 100, 0.1), nrow=dims$P, ncol=dims$L)
  
  return(storeInit)
}

set.seed(18704)

nreps <- 100
N <- 100
P <- 30
Q <- 10
nAx <- 5


RoGauRMSE <- matrix(NA, nreps, 3)
RoGauMAE <- matrix(NA, nreps, 3)
RoGauProc <- matrix(NA, nreps, 3)

for (rep in 1:nreps){
  X <- matrix(rnorm(Q*N), N, Q) |> scale()
  X <- cbind(rep(1, N), X)
  colnames(X) <- 1:(Q+1)
  B <- matrix(0, Q, nAx)
  while (any(colSums(B)==0)){
    B <- matrix(rnorm(Q*nAx, sd=2), Q, nAx)
    B[sample(1:(Q*nAx), round(0.3*(Q*nAx)), F)] <- 0
  }
  B <- rbind(rnorm(nAx, 10, 5), B)
  Eta <- X%*%B
  Lambda <- matrix(rnorm(P*nAx, 0.25, 2), P, nAx)
  Y <- Eta%*%t(Lambda) + matrix(rnorm(N*P), N, P) %*% diag(runif(P, 0.5, 10))
  Y <- ceiling(Y); Y[Y<0] <- 0
  
  
  Ycv <- Y; Ycv[as.matrix(reshape::melt(matrix(0,N,P))[sample(1:(P*N), round(0.2*P*N), F),1:2])] <- NA
  
  outGLLVMpois <- try(gllvm(y=Ycv, X=X, num.lv.c=nAx, 
                        num.lv = 0, family="poisson", starting.val="zero", maxit=1e7,
                        sd.errors=F, optimizer="alabama"))
  if(inherits(outGLLVMpois,"try-error")){outGLLVMpois <- try(gllvm(y=Ycv, X=X, num.lv.c=nAx, 
                              num.lv = 0, family="poisson", starting.val="zero", maxit=1e7,
                              sd.errors=F, optimizer="alabama"))}
  while(inherits(outGLLVMpois,"try-error")){outGLLVMpois <- try(gllvm(y=Ycv, X=X, num.lv.c=nAx, 
                              num.lv = 0, family="poisson", starting.val="zero", maxit=1e7,
                              sd.errors=F, optimizer="alabama"))}
  
  # outGLLVMnb <- try(gllvm(y=Ycv, X=X, num.lv.c=nAx, 
  #                   num.lv = 0, family="negative.binomial", starting.val="random", maxit=1e7,
  #                   sd.errors=F, optimizer="alabama"))
  # if(inherits(outGLLVMnb,"try-error")){outGLLVMnb <- try(gllvm(y=Ycv, X=X, num.lv.c=nAx, 
  #                                                              num.lv = 0, family="negative.binomial", starting.val="random", maxit=1e7,
  #                                                              sd.errors=F, optimizer="alabama"))}
  # while(inherits(outGLLVMnb,"try-error")){outGLLVMnb <- try(gllvm(y=Ycv, X=X, num.lv.c=nAx, 
  #                                                                 num.lv = 0, family="negative.binomial", starting.val="random", maxit=1e7,
  #                                                                 sd.errors=F, optimizer="alabama"))}
  
  dat <- list(Y=Ycv, X=X)
  dims <- list(N=nrow(dat$Y), P=ncol(dat$Y), L=nAx, Q=ncol(dat$X))
  storeInit <- getStoreInit(dat, dims)
  
  outBarcode <- sampler_snapdragonXonly(1000, 3, 3000, dims, dat, storeInit, 1870 %*% rep)
  Ebarcode <- apply(outBarcode$C*outBarcode$Phi, c(1,2), mean) %*% t(apply(outBarcode$S*outBarcode$Gamma, c(1,2), mean))
  
  RoGauRMSE[rep, 1] <- sqrt(mean((Ebarcode[is.na(Ycv)] - Y[is.na(Ycv)])^2))
  RoGauRMSE[rep, 2] <- sqrt(mean((exp(predict(outGLLVMpois))[is.na(Ycv)] - Y[is.na(Ycv)])^2))
  RoGauRMSE[rep, 3] <- sqrt(mean((exp(predict(outGLLVMnb))[is.na(Ycv)] - Y[is.na(Ycv)])^2))
  
  RoGauMAE[rep, 1] <- mean(abs(Ebarcode[is.na(Ycv)] - Y[is.na(Ycv)]))
  RoGauMAE[rep, 2] <- mean(abs(exp(predict(outGLLVMpois))[is.na(Ycv)] - Y[is.na(Ycv)]))
  RoGauMAE[rep, 3] <- mean(abs(exp(predict(outGLLVMnb))[is.na(Ycv)] - Y[is.na(Ycv)]))
  
  RoGauProc[rep, 1] <- procrustes(apply(outBarcode$C*outBarcode$Phi, c(1,2), mean), Eta)$ss
  RoGauProc[rep, 2] <- procrustes(exp(getLV(outGLLVMpois)), Eta)$ss
  RoGauProc[rep, 3] <- procrustes(exp(getLV(outGLLVMnb)), Eta)$ss
  
  print(paste("Rep:", rep))
}

saveRDS(list(RoGauRMSE=RoGauRMSE, 
             RoGauMAE=RoGauMAE,
             RoGauProc=RoGauProc),
        "./results/RoundGaussConcurrent")
library(NMF)
library(Rcpp)
library(tidyverse)
library(latex2exp)

source("SetupFunctions.R")
sourceCpp("./models/sampler_snapdragonStatic.cpp")

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

set.seed(18654)

nreps <- 100
N <- 100
P <- 30
nAx <- 5

pCopRMSE <- matrix(NA, nreps, 2)
pCopMAE <- matrix(NA, nreps, 2)

for (rep in 1:nreps){
  Lambda <- matrix(rnorm(nAx*P), nrow=P)
  Ycop <- matrix(rnorm(N*P), nrow=N) %*% (chol(cov2cor(Lambda%*%t(Lambda) + diag(P))))
  Y <- qpois(pnorm(Ycop), rgamma(N*P, 1, 1/10))
  Ycv <- Y; Ycv[as.matrix(reshape::melt(matrix(0,N,P))[sample(1:(P*N), round(0.2*P*N), F),1:2])] <- NA
  outGLLVM <- try(gllvm(y=Ycv, num.lv.c=0, 
               num.lv = 5, family="poisson", starting.val="res", maxit=1e7,
               sd.errors=F, optimizer="optim"))
  if (inherits(outGLLVM,"try-error")){try(gllvm(y=Ycv, num.lv.c=0, 
                                                num.lv = 5, family="poisson", starting.val="res", maxit=1e7,
                                                sd.errors=F, optimizer="optim"))}
  while (inherits(outGLLVM,"try-error")){try(gllvm(y=Ycv, num.lv.c=0, 
                                                   num.lv = 5, family="poisson", starting.val="res", maxit=1e7,
                                                   sd.errors=F, optimizer="optim"))}
  
  dat <- list(Y=Ycv)
  dims <- list(N=nrow(dat$Y), P=ncol(dat$Y), L=nAx)
  storeInit <- getStoreInit(dat, dims)
  
  outBarcode <- sampler_snapdragonStatic(1000, 3, 3000, dims, dat, storeInit, 4223 %*% rep)
  Ebarcode <- apply(outBarcode$C*outBarcode$Phi, c(1,2), mean) %*% t(apply(outBarcode$S*outBarcode$Gamma, c(1,2), mean))
  
  pCopRMSE[rep, 1] <- sqrt(mean((Ebarcode[is.na(Ycv)] - Y[is.na(Ycv)])^2))
  pCopRMSE[rep, 2] <- sqrt(mean((exp(predict(outGLLVM))[is.na(Ycv)] - Y[is.na(Ycv)])^2))
  
  pCopMAE[rep, 1] <- mean(abs(Ebarcode[is.na(Ycv)] - Y[is.na(Ycv)]))
  pCopMAE[rep, 2] <- mean(abs(exp(predict(outGLLVM))[is.na(Ycv)] - Y[is.na(Ycv)]))
  
  print(paste("Rep:", rep))
}

saveRDS(list(pCopRMSE=pCopRMSE, pCopMAE=pCopMAE), 
        "./results/PoissonCopulaUnconstrained.rds")
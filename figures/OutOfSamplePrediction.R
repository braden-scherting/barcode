library(gllvm)
source("SetupFunctions.R")
sourceCpp("sampler_snapdragonXonly.cpp")

getStoreInit <- function(dat, dims){
  storeInit <- list(Ypresent=1 * !is.na(dat$Y))
  
  blockSize <- min(dims$L, 5)
  storeInit$Cstar <- matrix(0, ncol=blockSize, nrow=2^blockSize)
  for (i in 1:(2^blockSize - 1)){
    storeInit$Cstar[i+1,] <- binaryLogic::as.binary(i, n=blockSize)
  }
  Ynmf <- dat$Y
  Ynmf[is.na(dat$Y)] <- colMeans(dat$Y, na.rm=T)[which(is.na(dat$Y), arr.ind = T)[,2]]
  
  PhiInit <- matrix(1/dims$N, nrow=dims$N, ncol=dims$L)
  GammaInit <- matrix(rgamma(dims$P*dims$L, 2,1), nrow=dims$P, ncol=dims$L)
  
  nonZrows <- which(rowSums(Ynmf)>0)
  nonZcols <- which(colSums(Ynmf)>0)
  
  nmfInit <- NMF::nmf(Ynmf[nonZrows, nonZcols], dims$L, nrun=10, seed=20114)
  PhiInit[nonZrows,] <- NMF::basis(nmfInit)
  GammaInit[nonZcols,] <- t(NMF::coef(nmfInit))
  # ab <- sqrt(mean(PhiInit) / mean(GammaInit))
  # PhiInit <- as.matrix(PhiInit / ab)
  # GammaInit <- t(as.matrix(GammaInit * ab))
  GammaInit <- GammaInit * matrix(colSums(PhiInit), dims$P, dims$L, byrow = T)
  PhiInit <- PhiInit / matrix(colSums(PhiInit), dims$N, dims$L, byrow = T)
  
  
  storeInit$PhiInit <- PhiInit
  storeInit$GammaInit <- GammaInit
  return(storeInit)
}


set.seed(331825)
nreps <- 25

tmp <- loadData(initialization = F)
N <- 300; P <- 30; L <- 5; holdoutFrac <- 0.25


finnRMSE <- matrix(NA, nreps, 3)
finnMAE <- matrix(NA, nreps, 3)

for (rep in 1:nreps){
  samps <- sample(1:tmp$dims$N, N, F)
  spec <- sample(1:tmp$dims$P, P, F)
  Y <- tmp$dat$Y[samps,spec]
  Xx <- tmp$dat$X[samps,]
  
  Ycv <- Y; Ycv[as.matrix(reshape::melt(matrix(0,N,P))[sample(1:(P*N), round(holdoutFrac*P*N), F),1:2])] <- NA  
  
  outPois <- gllvm(y=Ycv, X=Xx, num.lv.c=L, 
                   num.lv = 0, family="poisson", starting.val="res", maxit=1e7,
                   sd.errors=F, optimizer="alabama", control = list(reltol.c=1e-6))
  
  outNB <- gllvm(y=Ycv, X=Xx, num.lv.c=L, 
                 num.lv = 0, family="negative.binomial", starting.val="random", maxit=1e7,
                 sd.errors=F, optimizer="alabama", control = list(reltol.c=1e-6))
  
  dat <- list(Y=Ycv, X=Xx)
  dims <- list(N=nrow(dat$Y), P=ncol(dat$Y), L=L, Q=ncol(dat$X))
  storeInit <- getStoreInit(dat, dims)
  
  outBarcode <- sampler_snapdragonXonly(1000, 3, 3000, dims, dat, storeInit, 1870)
  Ebarcode <- apply(outBarcode$C*outBarcode$Phi, c(1,2), mean) %*% t(apply(outBarcode$S*outBarcode$Gamma, c(1,2), mean))
  
  finnRMSE[rep,1] <- sqrt(mean((Y[is.na(Ycv)] - Ebarcode[is.na(Ycv)])^2))
  finnRMSE[rep,2] <- sqrt(mean((Y[is.na(Ycv)] - exp(predict(outPois))[is.na(Ycv)])^2))
  finnRMSE[rep,3] <- sqrt(mean((Y[is.na(Ycv)] - exp(predict(outNB))[is.na(Ycv)])^2))
  
  finnMAE[rep,1] <- mean(abs(Y[is.na(Ycv)] - Ebarcode[is.na(Ycv)]))
  finnMAE[rep,2] <- mean(abs(Y[is.na(Ycv)] - exp(predict(outPois))[is.na(Ycv)]))
  finnMAE[rep,3] <- mean(abs(Y[is.na(Ycv)] - exp(predict(outNB))[is.na(Ycv)]))
  
  print(paste("Rep:", rep))
}

# saveRDS(list(finnRMSE=finnRMSE,
#              finnMAE=finnMAE), "results/OOPSEappN300P30R25.rds")


library(ggplot2)
library(ggpattern)
library(tidyverse)

OOPSE <- readRDS("results/OOPSEappN300P30R25.rds")

data.frame(reshape::melt(OOPSE)) %>%
  filter(X2<3) |> 
  ggplot() +
  geom_violin_pattern(aes(x=L1, y=log(value), pattern=as.factor(X2)), width=1, 
                      position=position_dodge(1), pattern_density=0.01, pattern_spacing=0.015, fill="white") + 
  geom_point(aes(x=L1, y=log(value), fill=as.factor(X2)), 
             position=position_jitterdodge(dodge.width = 1, jitter.width = 0.45), 
             shape="*", alpha=0.75, size=4, color="red", show.legend = F) + 
  scale_pattern_manual(values=c('wave', 'stripe'), labels=c("Barcode", "Pois. GLLVM")) + 
  ylab("log error") + xlab("") + ggtitle("Predicting community data") + 
  scale_x_discrete(labels= c("MAE", "RMSE")) + 
  labs(pattern=NULL) + theme_bw() 

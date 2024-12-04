library(NMF)
library(Rcpp)
library(tidyverse)
library(latex2exp)
library(pROC)

source("awwwSnap.R")
sourceCpp("sampler_snapdragonStatic.cpp")
# sourceCpp("sampler_snapdragon.cpp")

# ---FIXED VALUES--- # 
set.seed(1161140)

N <- 250
P <- 50
L <- 5

StrueArr <- array(dim=c(P, L, 3))
CtrueArr <- array(0, dim=c(N, L, 3))

pSC <- c(0.35, 0.5, 0.65)

for (xx in 1:3){
  StrueArr[,,xx] <- rbind(diag(L), diag(L), 
                          matrix(0, ncol=L, nrow=P-(2*L)))
  StrueArr[,,xx][cbind((2*L + 1):P, sample(1:L, P-(2*L), T))] <- 1
  StrueArr[(2*L + 1):P,,xx][sample(which(StrueArr[(2*L + 1):P,,xx]==0), P*L*(pSC[xx])-P)] <- 1
  
  CtrueArr[,,xx][cbind(1:N, sample(1:L, N, T))] <- 1
  CtrueArr[,,xx][sample(which(CtrueArr[,,xx]==0), N*L*(pSC[xx])-N)] <- 1
}

# structural sparsity: 0.12, 0.23, 0.48
nreps <- 25
PE00 <- matrix(nrow=nreps, ncol=3)
PE01 <- matrix(nrow=nreps, ncol=3)

for (xx in 1:3){
  for (rep in 1:nreps){
    Phixx <- matrix(rgamma(N*L, 1, 1/3), nrow=N)
    Gammaxx <- matrix(rgamma(P*L, 1, 1/5), nrow=P)
    Yxx <- matrix(nrow=N, ncol=P)
    
    ETxx <- (Phixx*CtrueArr[,,xx]) %*% t(Gammaxx*StrueArr[,,xx])
    
    for (i in 1:N){for (j in 1:P){ 
      Yxx[i,j] <- rpois(1, ETxx[i,j]) 
    }}
    
    dims <- list(N=N, P=P, L=L)
    dat <- list(Y=Yxx)
    
    storeInit <- list(Ypresent=matrix(1, nrow=dims$N, ncol=dims$P))
    
    blockSize <- min(dims$L, 5)
    storeInit$Cstar <- matrix(0, ncol=blockSize, nrow=2^blockSize)
    for (i in 1:(2^blockSize - 1)){
      storeInit$Cstar[i+1,] <- binaryLogic::as.binary(i, n=blockSize)
    }
    
    storeInit$PhiInit <- matrix(1/dims$N, nrow=dims$N, ncol=dims$L)
    storeInit$GammaInit <- matrix(rgamma(dims$P*dims$L, 100, 0.1), nrow=dims$P, ncol=dims$L)

    nmfInit <- NMF::nmf(dat$Y[rowSums(dat$Y)>0,colSums(dat$Y)>0], dims$L, nrun=10, seed=441503)
    PhiInit <- NMF::basis(nmfInit)
    GammaInit <- NMF::coef(nmfInit)
    ab <- colSums(PhiInit)
    PhiInit <- as.matrix(PhiInit / ab)
    GammaInit <- t(as.matrix(GammaInit * ab))
    
    storeInit$PhiInit[rowSums(dat$Y)>0,] <- PhiInit
    storeInit$GammaInit[colSums(dat$Y)>0] <- GammaInit
    
    
    out <- sampler_snapdragonStatic(2000, 1, 5000, dims, dat, storeInit, 612243 + (xx * rep^2))
    perm <- apply(apply(out$S, c(1,2), mean)[1:L,], 1, which.max)
    Ehat <- matrix(0, N, P)
    for (it in 1:2000){
      Ehat <- Ehat + ( (out$Phi[,,it]*out$C[,,it]) %*% t(out$Gamma[,,it]*out$S[,,it]) ==0 )
    }
    PE00[rep, xx] <- mean((Ehat/2000)[ETxx==0])
    PE01[rep, xx] <- mean((Ehat/2000)[ETxx!=0])
    # plot(apply(out$S[, perm, ], 3, function(x) mean(x == 1)), type="l", ylim=c(0.0,1))
    # lines(apply(out$C[, perm, ], 3, function(x) mean(x == 1)), type="l", ylim=c(0.0,01))
    # abline(h=pSC[xx], col="red")
  }
}

saveRDS(list(PE00, PE01), "results/PE.rds")

PE <- readRDS("results/PE.rds")
PE00 <- PE[[1]]; PE01 <- PE[[2]]

data.frame(y=c(PE00, PE01),
           spars=rep(sort(rep(c("a", "b", "c"), 25)), 2),
           E01=sort(rep(c("a0:0", "b0:+"), 75))) %>%
  ggplot() + 
  geom_boxplot(aes(x=spars, y=y, color=E01), fill="black", linewidth=0.45) + 
  theme_bw() + ylab("") + 
  # xlab(TeX(r"(Sparsity in C, S ($\mu$, induced) )")) +
  scale_x_discrete(name =TeX(r"(Sparsity in C, S ($\mu$, induced) )"), 
                   labels=c("a"="65% (48%)","b"="50% (23%)","c"="35% (9%)")) +
  ylim(-0.01, 1.01) +
  scale_color_manual(labels=c(TeX(r"( $\hat{Pr}(\mu_{ij}=0 | \mu^{*}_{ij}=0)$ )"), 
                              TeX(r"( $\hat{Pr}(\mu_{ij}=0 | \mu^{*}_{ij}>0)$ )")),
                     values=c("#009E73", "orange")) +
  labs(color='') + ggtitle("Recovering structural sparsity")



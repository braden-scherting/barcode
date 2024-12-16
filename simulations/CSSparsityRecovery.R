library(NMF)
library(Rcpp)
library(tidyverse)
library(latex2exp)
# library(pROC)

source("SetupFunctions.R")
sourceCpp("./models/sampler_snapdragonStatic.cpp")


# ---FIXED VALUES--- # 
set.seed(1171051)

Ns <- c(50, 100, 500, 1000)
Ps <- c(15, 30, 50, 75)
L <- 5

nreps <- 25


# rocC <- list(k1=list(),
#              k2=list(),
#              k3=list(),
#              k4=list())

avgMissC <- matrix(NA, nreps, 4)

for (k in 1:4){
  P <- Ps[k]
  N <- 500
  for (rep in 1:nreps){
    Srep <- rbind(diag(L), diag(L), 
                  matrix(0, ncol=L, nrow=P-(2*L)))
    while(min(rowSums(Srep[(2*L+1):P,] )) == 0){
      psiS <- rbeta(1, 5,5)
      Srep[(2*L+1):P,] <- sample(0:1, (P-(2*L))*5, T, c(1-psiS,psiS)) 
      # print(paste("S:", k , ",", rep) )
    }
    Crep <- matrix(0, N, L)
    while(min(rowSums(Crep )) == 0){
    etaC <- rbeta(L, 3, 3)
    etaC[1] <- rbeta(1, 10, 1)
    Crep <- cbind(sample(0:1, N, T, c(1-etaC[1], etaC[1])),
                  sample(0:1, N, T, c(1-etaC[2], etaC[2])),
                  sample(0:1, N, T, c(1-etaC[3], etaC[3])),
                  sample(0:1, N, T, c(1-etaC[4], etaC[4])),
                  sample(0:1, N, T, c(1-etaC[5], etaC[5])))
    }
    
    Phixx <- matrix(rgamma(N*L, 1, 1/3), nrow=N)
    Gammaxx <- matrix(rgamma(P*L, 1, 1/5), nrow=P)
    Yxx <- matrix(nrow=N, ncol=P)
    
    ETxx <- (Phixx*Crep) %*% t(Gammaxx*Srep)
    
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

    
    out <- sampler_snapdragonStatic(2000, 1, 3000, dims, dat, storeInit, 34289 + (k * rep^2))
    perm <- apply(apply(out$S, c(1,2), mean)[1:L,], 1, which.max)
    
    print((k-1)*nreps + rep )
    
    # rocC[[k]][[rep]] <- roc(c(Crep), c(apply(out$C[,perm,], c(1,2), mean)))$auc
    avgMissC[rep, k] <- mean(apply(out$C[,perm,], 3, function(x) mean((x!=Crep))))
  }
}

saveRDS(avgMissC, "./results/recoverC.rds")


# rocS <- list(k1=list(),
#              k2=list(),
#              k3=list(),
#              k4=list())

avgMissS <- matrix(NA, nreps, 4)

for (k in 1:4){
  P <- 50
  N <- Ns[k]
  for (rep in 1:nreps){
    Srep <- rbind(diag(L), diag(L), 
                  matrix(0, ncol=L, nrow=P-(2*L)))
    while(min(rowSums(Srep[(2*L+1):P,] )) == 0){
      psiS <- rbeta(1, 5,5)
      Srep[(2*L+1):P,] <- sample(0:1, (P-(2*L))*5, T, c(1-psiS,psiS)) 
      # print(paste("S:", k , ",", rep) )
    }
    Crep <- matrix(0, N, L)
    while(min(rowSums(Crep )) == 0){
      etaC <- rbeta(L, 3, 3)
      etaC[1] <- rbeta(1, 10, 1)
      Crep <- cbind(sample(0:1, N, T, c(1-etaC[1], etaC[1])),
                    sample(0:1, N, T, c(1-etaC[2], etaC[2])),
                    sample(0:1, N, T, c(1-etaC[3], etaC[3])),
                    sample(0:1, N, T, c(1-etaC[4], etaC[4])),
                    sample(0:1, N, T, c(1-etaC[5], etaC[5])))
    }
    
    Phixx <- matrix(rgamma(N*L, 1, 1/3), nrow=N)
    Gammaxx <- matrix(rgamma(P*L, 1, 1/5), nrow=P)
    Yxx <- matrix(nrow=N, ncol=P)
    
    ETxx <- (Phixx*Crep) %*% t(Gammaxx*Srep)
    
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
    
    
    out <- sampler_snapdragonStatic(1000, 2, 3000, dims, dat, storeInit, 4822243 + (k * rep^2))
    perm <- apply(apply(out$S, c(1,2), mean)[1:L,], 1, which.max)
    
    print((k-1)*nreps + rep )
    
    # rocS[[k]][[rep]] <- roc(c(Srep), c(apply(out$S[,perm,], c(1,2), mean)))$auc
    avgMissS[rep, k] <- mean(apply(out$S[,perm,], 3, function(x) mean((x!=Srep))))
  }
}

saveRDS(avgMissS, "./results/recoverS.rds")

avgMissC <- readRDS("./results/recoverC.rds")
avgMissS <- readRDS("./results/recoverS.rds")

data.frame(cm=colMeans(avgMissS), 
           cmax=apply(avgMissS,2,quantile, probs=0.95), 
           cmin=apply(avgMissS,2,quantile, probs=0.05), ss=c(50, 100, 250, 500)) %>%
  ggplot() + 
  geom_jitter(data=data.frame(cm=c(avgMissS), ss=sort(rep(c(50, 100, 250, 500),nreps))), aes(y=cm, x=as.factor(ss)), width=0.075, alpha=0.5, col="red", shape=4) +
  geom_segment(aes(x=as.factor(ss), y=cmin, yend=cmax), col="black", linewidth=0.75) +
  geom_line(aes(x=as.factor(ss), y=cm, group=1)) +
  geom_point(aes(x=as.factor(ss), y=cm), shape=23, color="black", fill="red", size=2) + 
  geom_point(aes(x=as.factor(ss), y=cmax), shape="—") + 
  geom_point(aes(x=as.factor(ss), y=cmin), shape="—") + 
  theme_bw() + ylab(TeX(r"( RMSE $\hat{S}$ )")) + xlab(TeX(r"( Sample size: $n$ )")) +
  ylim(-0.01, 0.25) + 
  ggtitle(TeX(r"( $$ )")) 


data.frame(cm=colMeans(avgMissC), 
           cmax=apply(avgMissC,2,quantile, probs=0.95), 
           cmin=apply(avgMissC,2,quantile, probs=0.05), ss=c(15, 30, 50, 75)) %>%
  ggplot() + 
  geom_jitter(data=data.frame(cm=c(avgMissC), ss=sort(rep(c(15, 30, 50, 75),nreps))), aes(y=cm, x=as.factor(ss)), width=0.075, alpha=0.5, col="red", shape=4) +
  geom_segment(aes(x=as.factor(ss), y=cmin, yend=cmax), col="black", linewidth=0.75) +
  geom_line(aes(x=as.factor(ss), y=cm, group=1)) +
  geom_point(aes(x=as.factor(ss), y=cm), shape=23, color="black", fill="red", size=2) + 
  geom_point(aes(x=as.factor(ss), y=cmax), shape="—") + 
  geom_point(aes(x=as.factor(ss), y=cmin), shape="—") + 
  theme_bw() + ylab(TeX(r"( RMSE $\hat{C}$ )")) + xlab(TeX(r"( Sample size: $p$ )")) +
  ylim(-0.01, 0.25) + 
  ggtitle(TeX(r"( $$ )")) 

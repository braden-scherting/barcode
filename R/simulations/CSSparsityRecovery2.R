library(tidyverse)
library(latex2exp)
library(scales)

source("setup_functions.R")
sourceCpp("../src/barcodeSampler.cpp")


# ---FIXED VALUES--- # 
set.seed(1171051)

Ns <- c(50, 100, 500, 1000)
Ps <- c(15, 30, 50, 75)
L <- 5+1

nreps <- 25

avgMissC <- matrix(NA, nreps, 4)
medMissC <- matrix(NA, nreps, 4)

PC10 <- matrix(NA, nreps, 4)
PC11 <- matrix(NA, nreps, 4)

for (k in 1:4){
  P <- Ps[k]
  N <- 500
  for (rep in 1:nreps){
    Srep <- rbind(diag(L), diag(L), 
                  matrix(0, ncol=L, nrow=P-(2*L)))
    while(min(rowSums(Srep[(2*L+1):P,] )) == 0){
      psiS <- rbeta(1, 5,5)
      Srep[(2*L+1):P,] <- sample(0:1, (P-(2*L))*L, T, c(1-psiS,psiS)) 
    }
    Crep <- matrix(0, N, L)
    while(min(rowSums(Crep )) == 0){
    etaC <- rbeta(L, 3, 3)
    etaC[1] <- rbeta(1, 10, 1)
    Crep <- cbind(rep(1, N),
                  sample(0:1, N, T, c(1-etaC[1], etaC[1])),
                  sample(0:1, N, T, c(1-etaC[2], etaC[2])),
                  sample(0:1, N, T, c(1-etaC[3], etaC[3])),
                  sample(0:1, N, T, c(1-etaC[4], etaC[4])),
                  sample(0:1, N, T, c(1-etaC[5], etaC[5])))
    }
    
    Phixx <- cbind(rep(1/N, N),matrix(rgamma(N*(L-1), 1, 1/3), nrow=N))
    Gammaxx <- matrix(rgamma(P*L, 1, 1/5), nrow=P)
    Yxx <- matrix(nrow=N, ncol=P)
    
    ETxx <- (Phixx*Crep) %*% t(Gammaxx*Srep)
    
    for (i in 1:N){for (j in 1:P){ 
      Yxx[i,j] <- rpois(1, ETxx[i,j]) 
    }}
    dims <- list(N=N, P=P, L=L, Q=1)
    
    ij <- which(Yxx!=0, arr.ind = T)
    RspYxx <- Matrix::sparseMatrix(i=ij[,1], j=ij[,2], x=Yxx[Yxx!=0], dims=c(N, P))
    dat <- list(RspY=RspYxx, X=matrix(1, nrow=N, ncol=1))
    
    storeInit <- list()
    
    blockSize <- min(dims$L, 5)
    storeInit$Cstar <- matrix(0, ncol=blockSize, nrow=2^blockSize)
    for (i in 1:(2^blockSize - 1)){
      storeInit$Cstar[i+1,] <- binaryLogic::as.binary(i, n=blockSize)
    }
    
    mf <- poismf(Yxx, k=dims$L - 1, l2_reg=0, niter=100, nthreads = 1)
    storeInit$Ginit <- cbind(rep(1, dims$P), t(mf$B) * (rep(1, dims$P) %*% t(rowSums(mf$A))))
    storeInit$Pinit <- cbind(rep(1/dims$N, dims$N), t(mf$A) / (rep(1, dims$N) %*% t(rowSums(mf$A))))

    out <- sampler_camellia_ordination(dat, dims, storeInit, 612243 + (k * rep^2),
                            25, 1, 500)
    
    perm <- apply(apply(out$S, c(1,2), mean)[2:L,], 1, which.max)
    
    print((k-1)*nreps + rep )
    
    avgMissC[rep, k] <- mean(apply(out$C[,perm,], 3, function(x) mean((x!=Crep[,2:L]))))
    medMissC[rep, k] <- mean(apply(out$C[,perm,], c(1,2), median)!=Crep[,2:L])
    
    PC10[rep,k] <- mean(apply(out$C[,perm,], c(1,2), mean)[Crep[,2:L]==0])
    PC11[rep,k] <- mean(apply(out$C[,perm,], c(1,2), mean)[Crep[,2:L]==1])
  }
}

saveRDS(avgMissC, "../output/results/recoverC.rds")

avgMissS <- matrix(NA, nreps, 4)
medMissS <- matrix(NA, nreps, 4)

for (k in 1:4){
  P <- 50
  N <- Ns[k]
  for (rep in 1:nreps){
    Srep <- rbind(diag(L), diag(L), 
                  matrix(0, ncol=L, nrow=P-(2*L)))
    while(min(rowSums(Srep[(2*L+1):P,] )) == 0){
      psiS <- rbeta(1, 5,5)
      Srep[(2*L+1):P,] <- sample(0:1, (P-(2*L))*L, T, c(1-psiS,psiS)) 
      # print(paste("S:", k , ",", rep) )
    }
    Crep <- matrix(0, N, L)
    while(min(rowSums(Crep )) == 0){
      etaC <- rbeta(L, 3, 3)
      etaC[1] <- rbeta(1, 10, 1)
      Crep <- cbind(rep(1, N),
                    sample(0:1, N, T, c(1-etaC[1], etaC[1])),
                    sample(0:1, N, T, c(1-etaC[2], etaC[2])),
                    sample(0:1, N, T, c(1-etaC[3], etaC[3])),
                    sample(0:1, N, T, c(1-etaC[4], etaC[4])),
                    sample(0:1, N, T, c(1-etaC[5], etaC[5])))
    }
    Crep[,1] <- 1
    
    Phixx <- cbind(rep(1/N, N), matrix(rgamma(N*(L-1), 1, 1/3), nrow=N))
    Gammaxx <- matrix(rgamma(P*L, 1, 1/5), nrow=P)
    Yxx <- matrix(nrow=N, ncol=P)
    
    ETxx <- (Phixx*Crep) %*% t(Gammaxx*Srep)
    
    for (i in 1:N){for (j in 1:P){ 
      Yxx[i,j] <- rpois(1, ETxx[i,j]) 
    }}
    dims <- list(N=N, P=P, L=L, Q=1)
    
    ij <- which(Yxx!=0, arr.ind = T)
    RspYxx <- Matrix::sparseMatrix(i=ij[,1], j=ij[,2], x=Yxx[Yxx!=0], dims=c(N, P))
    dat <- list(RspY=RspYxx, X=matrix(1, nrow=N, ncol=1))
    
    storeInit <- list()
    
    blockSize <- min(dims$L, 5)
    storeInit$Cstar <- matrix(0, ncol=blockSize, nrow=2^blockSize)
    for (i in 1:(2^blockSize - 1)){
      storeInit$Cstar[i+1,] <- binaryLogic::as.binary(i, n=blockSize)
    }
    
    mf <- poismf(Yxx, k=dims$L - 1, l2_reg=0, niter=100, nthreads = 1)
    storeInit$Ginit <- cbind(rep(1, dims$P), t(mf$B) * (rep(1, dims$P) %*% t(rowSums(mf$A))))
    storeInit$Pinit <- cbind(rep(1/dims$N, dims$N), t(mf$A) / (rep(1, dims$N) %*% t(rowSums(mf$A))))
    
    out <- sampler_camellia_ordination(dat, dims, storeInit, 612243 + (k * rep^2),
                                       2000, 1, 2000)
    perm <- apply(apply(out$S, c(1,2), mean)[1:L,], 1, which.max)
    
    print((k-1)*nreps + rep )
    
    avgMissS[rep, k] <- mean(apply(out$S[,perm,], 3, function(x) mean((x!=Srep))))
    medMissS[rep, k] <- mean(apply(out$S[,perm,], c(1,2), median)!=Srep)
  }
}

saveRDS(avgMissS, "../output/results/recoverS.rds")

avgMissC <- readRDS("../output/results/recoverC.rds")
avgMissS <- readRDS("../output/results/recoverS.rds")

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
  scale_y_continuous(labels=number_format(0.01), limits=c(-0.01, 0.3))

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
  scale_y_continuous(labels=number_format(0.01), limits=c(-0.01, 0.1)) +
  theme(axis.ticks.y=)


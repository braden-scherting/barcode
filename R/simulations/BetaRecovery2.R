library(tidyverse)
library(latex2exp)
library(scales)

source("setup_functions.R")
sourceCpp("../src/barcodeSampler.cpp")

# ---FIXED VALUES--- # 
set.seed(1097329)

P <- 50
L <- 6+1
Q <- 5

Strue <- rbind(diag(L), 
               diag(L), 
               matrix(0, 
                      ncol=L, 
                      nrow=P-(2*L)))
pS <- rbeta(L, 10, 9)
for (l in 1:L){
  Strue[(2*L + 1):P, l] <- sample(0:1, P-(2*L), T, c(1-pS[l], pS[l]))
}

tmpX <- cbind(1, matrix(rnorm((Q-1)*1000, sd=0.5), nrow=1000))

TBeta <- matrix(runif((L-1)*(Q), -1., 1.), ncol=(L-1))
PrCfull <-  pnorm(tmpX %*% TBeta)

# ------------------ # 

Ns <- c(100, 250, 500, 1000)
nreps <- 25

BigBeta <- array(dim=c(3, Q, L-1, length(Ns), nreps))

for (k in 1:length(Ns)){
  N <- Ns[k]
  nset <- sample(1:1000, N, F)
  for (rep in 1:nreps){
    TC <- 1 * (matrix(runif(N*(L-1)), nrow=N) < PrCfull[nset,])
    TC <- cbind(rep(1/N, N), TC)
    TGamma <- matrix(rgamma(P * L, 1, 1/5), nrow=P)
    TPhi <- cbind(rep(1/N, N), matrix(rgamma(N * (L-1), 1, 1/3), nrow=N))
    
    ET <- (TPhi*TC) %*% t(TGamma * Strue)
    
    Y <- matrix(nrow=N, ncol=P)
    for (i in 1:N){for (j in 1:P){Y[i,j] <- rpois(1,ET[i,j])}}
    
    ij <- which(Y!=0, arr.ind = T)
    RspYxx <- Matrix::sparseMatrix(i=ij[,1], j=ij[,2], x=Y[Y!=0], dims=c(N, P))
    
    dims <- list(N=N, P=P, Q=Q, L=L)
    dat <- list(X=tmpX[nset,], RspY=RspYxx)
    
    storeInit <- list()
    
    blockSize <- min(dims$L, 5)
    storeInit$Cstar <- matrix(0, ncol=blockSize, nrow=2^blockSize)
    for (i in 1:(2^blockSize - 1)){
      storeInit$Cstar[i+1,] <- binaryLogic::as.binary(i, n=blockSize)
    }
    
    storeInit$Pinit <- matrix(1/dims$N, nrow=dims$N, ncol=dims$L)
    storeInit$Ginit <- matrix(rgamma(dims$P*dims$L, 1000, 0.1), nrow=dims$P, ncol=dims$L)
    
    # out <- sampler_snapdragonXonly(5000, 1, 2500, dims, dat, storeInit, 8115201 + (k * rep^2))
    out <- sampler_camellia_ordination(dat, dims, storeInit, 8115201 + (k * rep^2),
                                       2000, 1, 2000)
    perm <- apply(apply(out$S, c(1,2), mean)[2:L,], 1, which.max)
    
    BigBeta[,,,k,rep] <- apply(out$Beta[,perm,], c(1,2), 
                               quantile, probs=c(0.025, 0.5, 0.975))
    print(paste(k,":", rep))
  }
}

for (l in 1:(L-1)){
  plot(out$Beta[1,perm[l],], type="l", col=2, ylim=c(-3,3), main=paste(l))
  abline(h=TBeta[1,l], col=2)
  for (q in 2:Q){
    lines(out$Beta[q,perm[l],], col=q+1)
    abline(h=TBeta[q,l], col=q+1)
  }
}

saveRDS(BigBeta, "../output/results/BigBeta2.rds")

BigBeta <- readRDS("../output/results/BigBeta2.rds")
# BigBeta <- readRDS("~/Documents/Duke/Research/barcodeRevision/output/results/BigBeta2.rds")

myl=2; myq=2; data.frame(med=c(BigBeta[2,myq, myl,,]),
                         low=c(BigBeta[3,myq, myl,,]),
                         high=c(BigBeta[1,myq, myl,,]),
                         ss=c(rep(Ns, nreps)),
                         reps=rep( seq(1, 50, 2), 4) +  (rep(c(-0.25,-0.1, 0.1, 0.25), nreps))) %>% 
  ggplot(aes(x=reps, group=reps)) + 
  geom_hline(yintercept = TBeta[myq, myl]) + 
  geom_point(aes(color=as.factor(ss), y=med), shape=18) +
  geom_point(aes(color=as.factor(ss), y=high), shape="—") +
  geom_point(aes(color=as.factor(ss), y=low), shape="—") +
  geom_segment(aes(y = low, yend=high, color=as.factor(ss))) + 
  ylim(-1.5,1.5) + scale_color_manual(values=c("blue", "#D55E00", "black", "orange")) + 
  theme_bw() + ggtitle("Across replicates")


# myl <- l
for (myl in 1:(L-1)){
  pllt <- data.frame(med=c(BigBeta[2,, myl,,6]),
                     low=c(BigBeta[3,, myl,,6]),
                     high=c(BigBeta[1,, myl,,6]),
                     qid=rep(1:Q, 4) + sort(rep(c(-0.3, -0.1,0.1,0.3),Q)),
                     ss=sort(rep(Ns, Q))) %>% 
    ggplot(aes(x=qid, color=as.factor(ss))) + 
    annotate("segment", x=1-0.35, xend=1+0.35, y=TBeta[1,myl], color="black") +
    annotate("segment", x=2-0.35, xend=2+0.35, y=TBeta[2,myl], color="black") +
    annotate("segment", x=3-0.35, xend=3+0.35, y=TBeta[3,myl], color="black") +
    annotate("segment", x=4-0.35, xend=4+0.35, y=TBeta[4,myl], color="black") +
    annotate("segment", x=5-0.35, xend=5+0.35, y=TBeta[5,myl], color="black") +
    geom_point(aes(y=med), show.legend=T) + 
    geom_segment(aes(y = low, yend=high), show.legend=F) +
    geom_point(aes(color=as.factor(ss), y=high), shape="—", show.legend=F) +
    geom_point(aes(color=as.factor(ss), y=low), shape="—", show.legend=F) +
    scale_color_manual(values=c("blue", "#D55E00", "#009E73", "orange")) +
    ylim(-1.7,1.6) + theme_bw() + ylab("") + xlab("") + coord_fixed() + 
    theme(legend.box.background = element_rect(colour = "black") ) +
    labs(color='Sample Size') +
    ggtitle(paste("(",myl,")", sep=""))
  print(pllt)
}



BBcover <- matrix(0, nrow=Q-1, ncol=length(Ns))

BBdev <- array(0, dim=c(nreps, length(Ns)))
BBdevNoInt <- array(0, dim=c(nreps, length(Ns)))

for (k in 1:length(Ns)){
  for (rep in 1:nreps){
    BBdev[rep, k] <-  sqrt(mean((BigBeta[2, 1:Q, ,k, rep] - TBeta[1:Q,])^2))
    BBdevNoInt[rep, k] <-  sqrt(mean((BigBeta[2, 2:Q, ,k, rep] - TBeta[2:Q,])^2))
  }
}

data.frame(cm=colMeans(BBdev), cmax=apply(BBdev,2,max), cmin=apply(BBdev,2,min), ss=c(100, 250, 500, 1000)) %>%
  ggplot() + 
  geom_segment(aes(x=as.factor(ss), y=cmin, yend=cmax), col="black", linewidth=0.75) +
  geom_jitter(data=data.frame(cm=c(BBdev), ss=sort(rep(c(100, 250, 500, 1000),nreps))), aes(y=cm, x=as.factor(ss)), width=0.075, alpha=0.5, col="red", shape=4) +
  geom_line(aes(x=as.factor(ss), y=cm, group=1)) +
  geom_point(aes(x=as.factor(ss), y=cm), shape=23, color="black", fill="red", size=2) + 
  geom_point(aes(x=as.factor(ss), y=cmax), shape="—") + 
  geom_point(aes(x=as.factor(ss), y=cmin), shape="—") + 
  theme_bw() + ylab(TeX(r"( RMSE $\hat{B}$ )")) + xlab(TeX(r"( Sample size: $n$ )")) +
  scale_y_continuous(labels=number_format(0.01), limits=c(-0.01, 0.4)) +
  ggtitle(TeX(r"( $$ )")) 

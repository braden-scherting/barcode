library(tidyverse)
library(latex2exp)

source("SetupFunctions.R")

chains <- readRDS("results/birds6factors05Dec2024.rds")
tmp <- loadData(initialization = F, noFactors = 6)

# chains <- readRDS("results/spatial10factors06Dec2024INCOMPLETE.rds")
# tmp <- loadData(initialization = F, noFactors = 10)

# chains <- readRDS("results/spatial0424Wed11Dec2024.rds")
# tmp <- loadData(initialization = F, noFactors = 8)

n <- nrow(chains[[1]]$C)
p <- nrow(chains[[1]]$S)
niter <- ncol(chains[[1]]$bPhi)

for (chain in 1:4){
  chains[[chain]]$Ypart <- NULL
  # chains[[chain]]$E <- array(-1, dim=c(n, p, niter))
  
  chains[[chain]]$LL <- array(0, dim=c(niter))
  
  for (it in 1:niter){
    chains[[chain]]$LL[it] <- sum(dpois(
      tmp$dat$Y,
      (chains[[chain]]$C[,,it]*chains[[chain]]$Phi[,,it]) %*% 
      t(chains[[chain]]$S[,,it]*chains[[chain]]$Gamma[,,it]), log=T ))
  }
}

coda::effectiveSize(cbind(chains[[1]]$LL,
                          chains[[2]]$LL,
                          chains[[3]]$LL,
                          chains[[4]]$LL))

posterior::rhat(cbind(chains[[1]]$LL,
                          chains[[2]]$LL,
                          chains[[3]]$LL,
                          chains[[4]]$LL))

Bhat <- matrix(nrow=tmp$dims$Q, ncol=tmp$dims$L)
for (q in 1:tmp$dims$Q){
  for (l in 1:tmp$dims$L){
    Bhat[q,l] <- posterior::rhat(cbind(chains[[1]]$Beta[q,l,],
                          chains[[2]]$Beta[q,l,],
                          chains[[3]]$Beta[q,l,],
                          chains[[4]]$Beta[q,l,]))
  }
}

GShat <- matrix(nrow=tmp$dims$P, ncol=tmp$dims$L)
for (p in 1:tmp$dims$P){
  for (l in 1:tmp$dims$L){
    GShat[p,l] <- posterior::rhat(cbind(chains[[1]]$Gamma[p,l,]*chains[[1]]$S[p,l,],
                                        chains[[2]]$Gamma[p,l,]*chains[[2]]$S[p,l,],
                                        chains[[3]]$Gamma[p,l,]*chains[[3]]$S[p,l,],
                                        chains[[4]]$Gamma[p,l,]*chains[[4]]$S[p,l,]))
  }
}

PChat <- matrix(nrow=tmp$dims$N, ncol=tmp$dims$L)
for (n in 1:tmp$dims$N){
  for (l in 1:tmp$dims$L){
    PChat[n,l] <- posterior::rhat(cbind(chains[[1]]$Phi[n,l,]*chains[[1]]$C[n,l,],
                                        chains[[2]]$Phi[n,l,]*chains[[2]]$C[n,l,],
                                        chains[[3]]$Phi[n,l,]*chains[[3]]$C[n,l,],
                                        chains[[4]]$Phi[n,l,]*chains[[4]]$C[n,l,]))
  }
}

Xihat <- matrix(nrow=tmp$dims$K, ncol=tmp$dims$L)
for (k in 1:tmp$dims$K){
  for (l in 1:tmp$dims$L){
    Xihat[k,l] <- posterior::rhat(cbind(chains[[1]]$Xi[k,l,],
                                       chains[[2]]$Xi[k,l,],
                                       chains[[3]]$Xi[k,l,],
                                       chains[[4]]$Xi[k,l,]))
  }
}

combined <- stackedList(chains)

coda::effectiveSize(combined$LL)

# E.ess <- apply(combined$E, c(1,2), coda::effectiveSize)
B.ess <- apply(combined$Beta, c(1,2), coda::effectiveSize)

GS.ess <- apply(combined$Gamma*combined$S, c(1,2), coda::effectiveSize)
PC.ess <- apply(combined$Phi*combined$C, c(1,2), coda::effectiveSize)



data.frame(ll = c(chains[[1]]$LL,
                  chains[[2]]$LL,
                  chains[[3]]$LL,
                  chains[[4]]$LL), 
      chain = as.factor(sort(rep(1:4, 1500))),
      it = rep(1:1500, 4)) %>% 
  ggplot() +
  geom_line(aes(x=it, y=ll, color=chain), linetype="dotted") + 
  scale_color_manual(values=c("blue", "#D55E00", "#009E73", "orange")) +
  xlab("Thinned iter.") + ylab("Log-likelihood") + 
  ggtitle("Two modes") + theme_bw() 

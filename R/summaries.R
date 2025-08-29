library(tidyverse)
library(latex2exp)

source("setup_functions.R")
# samples <- stackedList(readRDS("../output/results/rank7_1948Wed25Jun2025.rds"))
samples <- stackedList(readRDS("../output/results/MYRESULTS"))
data <- loadData(noFactors = ncol(samples$Gamma))

# Reorder factors
oldOrder <- c(1, 2, 3, 4, 5, 6, 7)
newOrder <- c(1, 2, 5, 6, 7, 3, 4)

Cmed <- apply(samples$C, c(1,2), median)[,newOrder]
Smed <- apply(samples$S, c(1,2), median)[,newOrder]

round(colMeans(Cmed)*100)
round(colMeans(Smed)*100)

# Site fidelity
onOff <- (t(data$dat$des01) %*% Cmed) / rowSums(t(data$dat$des01))
round(colMeans(onOff==0 | onOff==1)*100)

# Site clusters
siteCmed <- round((t(data$dat$des01) %*% Cmed[,2:7]) / rowSums(t(data$dat$des01)))

siteCmed <- array(dim=c(555, 6, dim(samples$C)[3]))
for (it in 1:dim(samples$C)[3]){
  siteCmed[,,it] <- (t(data$dat$des01) %*% samples$C[,2:7,it]) / rowSums(t(data$dat$des01))
}

siteCmed2 <- round(apply(siteCmed, c(1,2), median))
siteCmed3 <- 1.0*(apply(siteCmed, c(1,2), median)==1)
siteCmed4 <- apply(apply(siteCmed, c(1,2), mean), 1, which.max)

siteCPhimed <- array(dim=c(555, 6, dim(samples$C)[3]))
for (it in 1:dim(samples$C)[3]){
  siteCPhimed[,,it] <- (t(data$dat$des01) %*% (samples$C[,2:7,it]*samples$Phi[,2:7,it])) 
}
siteMax <- apply(apply(siteCPhimed, c(1,2), mean), 1, which.max)
siteMaxClust <- apply(apply(siteCPhimed, c(1,2), median), c(1), function(x)paste(1.*(x>0), collapse = ""))



table(apply(siteCmed2, 1, function(x) paste(x, collapse = "")))

# Factors 3 and 6
round(mean(Cmed[,3]==Cmed[,6] & Cmed[,6]==1)*100)
sum(Smed[,3]==Smed[,6] & Smed[,6]==1)
round(mean(Smed[,3]==Smed[,6] & Smed[,6]==1)*100)

sum(Smed[,3]==Smed[,6] & Smed[,6]==1 & Smed[,6]==Smed[,7])


colnames(data$dat$Y)[which(rowSums(Smed)==1 & Smed[,1]==1)]
colnames(data$dat$Y)[which(rowSums(Smed)==1 & Smed[,2]==1)]
colnames(data$dat$Y)[which(rowSums(Smed)==1 & Smed[,3]==1)]
colnames(data$dat$Y)[which(rowSums(Smed)==1 & Smed[,4]==1)]
colnames(data$dat$Y)[which(rowSums(Smed)==1 & Smed[,5]==1)]
colnames(data$dat$Y)[which(rowSums(Smed)==1 & Smed[,6]==1)]
colnames(data$dat$Y)[which(rowSums(Smed)==1 & Smed[,7]==1)]

colnames(data$dat$Y)[which(rowSums(Smed)==1 & Smed[,1]==1)]
colnames(data$dat$Y)[which(rowSums(Smed[,2:7])==1 & Smed[,2]==1)]
colnames(data$dat$Y)[which(rowSums(Smed[,2:7])==1 & Smed[,3]==1)]
colnames(data$dat$Y)[which(rowSums(Smed[,2:7])==1 & Smed[,4]==1)]
colnames(data$dat$Y)[which(rowSums(Smed[,2:7])==1 & Smed[,5]==1)]
colnames(data$dat$Y)[which(rowSums(Smed[,2:7])==1 & Smed[,6]==1)]
colnames(data$dat$Y)[which(rowSums(Smed[,2:7])==1 & Smed[,7]==1)]


colnames(data$dat$Y)[which(rowSums(Smed)==2 & Smed[,3]==1 & Smed[,4]==1)]
colnames(data$dat$Y)[which(rowSums(Smed)==2 & Smed[,2]==1 & Smed[,5]==1)]
colnames(data$dat$Y)[which(rowSums(Smed)==2 & Smed[,2]==1 & Smed[,6]==1)]
colnames(data$dat$Y)[which(rowSums(Smed)==2 & Smed[,7]==1)]

length(unique(apply(Smed, 1, function(x) paste(as.character(x), collapse = ""))))
table(apply(Smed, 1, function(x) paste(as.character(x), collapse = "")))

colnames(data$dat$Y)[which(Smed[,1]==Smed[,2] & Smed[,3]==Smed[,4] & Smed[,2] + Smed[,3]==2)]
Smed[which(Smed[,1]==Smed[,2] & Smed[,3]==Smed[,4] & Smed[,2] + Smed[,3]==2),]

colnames(data$dat$Y)[which(rowSums(Smed)==7)]

# MCMC diagnostics
# chains <- readRDS("../output/results/rank7_1948Wed25Jun2025.rds")
chains <- readRDS("../output/results/MYRESULTS")

n <- nrow(chains[[1]]$C)
p <- nrow(chains[[1]]$S)
niter <- nrow(chains[[1]]$psi)

for (chain in 1:4){
  chains[[chain]]$Ypart <- NULL
  chains[[chain]]$E <- array(-1, dim=c(n, p, niter))

  chains[[chain]]$LL <- array(0, dim=c(niter))

  for (it in 1:niter){
    chains[[chain]]$LL[it] <- sum(dpois(
      data$dat$Y,
      (chains[[chain]]$C[,,it]*chains[[chain]]$Phi[,,it]) %*%
        t(chains[[chain]]$S[,,it]*chains[[chain]]$Gamma[,,it]), log=T ))
  }
}

coda::effectiveSize(cbind(chains[[1]]$LL,
                          chains[[2]]$LL,
                          chains[[3]]$LL,
                          chains[[4]]$LL))

coda::effectiveSize(c(chains[[1]]$LL,
                          chains[[2]]$LL,
                          chains[[3]]$LL,
                          chains[[4]]$LL))

posterior::rhat(cbind(chains[[1]]$LL,
                      chains[[2]]$LL,
                      chains[[3]]$LL,
                      chains[[4]]$LL))

Bhat <- matrix(nrow=data$dims$Q, ncol=data$dims$L-1)
for (q in 1:data$dims$Q){
  for (l in 2:data$dims$L){
    Bhat[q,l-1] <- posterior::rhat(cbind(chains[[1]]$Beta[q,l,],
                                       chains[[2]]$Beta[q,l,],
                                       chains[[3]]$Beta[q,l,],
                                       chains[[4]]$Beta[q,l,]))
  }
}
quantile(Bhat, c(0.025,0.975))

GShat <- matrix(nrow=data$dims$P, ncol=data$dims$L)
for (p in 1:data$dims$P){
  for (l in 1:data$dims$L){
    GShat[p,l] <- posterior::rhat(cbind(chains[[1]]$Gamma[p,l,]*chains[[1]]$S[p,l,],
                                        chains[[2]]$Gamma[p,l,]*chains[[2]]$S[p,l,],
                                        chains[[3]]$Gamma[p,l,]*chains[[3]]$S[p,l,],
                                        chains[[4]]$Gamma[p,l,]*chains[[4]]$S[p,l,]))
  }
}
quantile(GShat, c(0.025,0.975))

PChat <- matrix(nrow=data$dims$N, ncol=data$dims$L-1)
for (n in 1:data$dims$N){
  for (l in 2:data$dims$L){
    PChat[n,l-1] <- posterior::rhat(cbind(chains[[1]]$Phi[n,l,]*chains[[1]]$C[n,l,],
                                        chains[[2]]$Phi[n,l,]*chains[[2]]$C[n,l,],
                                        chains[[3]]$Phi[n,l,]*chains[[3]]$C[n,l,],
                                        chains[[4]]$Phi[n,l,]*chains[[4]]$C[n,l,]))
  }
}
quantile(PChat, c(0.025,0.975), na.rm = T)

Xihat <- matrix(nrow=data$dims$K, ncol=data$dims$L-1)
for (k in 1:data$dims$K){
  for (l in 2:data$dims$L){
    Xihat[k,l-1] <- posterior::rhat(cbind(chains[[1]]$Xi[k,l,],
                                        chains[[2]]$Xi[k,l,],
                                        chains[[3]]$Xi[k,l,],
                                        chains[[4]]$Xi[k,l,]))
  }
}
quantile(Xihat, c(0.025,0.975))

data.frame(ll = c(chains[[1]]$LL,
                  chains[[2]]$LL,
                  chains[[3]]$LL,
                  chains[[4]]$LL), 
           chain = as.factor(sort(rep(1:4, 2500))),
           it = rep(1:2500, 4)) %>% 
  ggplot() +
  geom_line(aes(x=it, y=ll, color=chain), linetype="dotted") + 
  scale_color_manual(values=c("blue", "#D55E00", "#009E73", "orange")) +
  xlab("Thinned iter.") + ylab("Log-likelihood") + 
  theme_bw() 

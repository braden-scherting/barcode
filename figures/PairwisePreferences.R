library(corrplot)

source("SetupFunctions.R")
outStack <- stackedList(readRDS("./results/birds6factors05Dec2024.rds"))

Cmed <- apply(outStack$C[,2:6,], c(1,2), median)
Cstr <- apply(Cmed, 1, function(x) paste(x, collapse=""))

Smed <- apply(outStack$S, c(1,2), median)
Sstr <- apply(Smed, 1, function(x) paste(x, collapse=""))

P <- tmp$dims$P
SstrMat <- matrix(0, P,P)
for (j in 1:P){
  SstrMat[j,] <- Sstr[j]==Sstr
}

colnames(SstrMat) <- colnames(tmp$dat$Y)
rownames(SstrMat) <- colnames(tmp$dat$Y)

c <- corrplot::corrplot(SstrMat, method="color", order="hclust", 
                        is.corr = F, col = COL1('YlOrRd', 6))

for (j in 1:(P-1)){
  for (k in (j):(P)){
    tmpjk <- 0
    for (it in 1:(dim(outStack$S)[3])){
      tmpjk <- tmpjk + sum(outStack$S[j,,it]==outStack$S[k,,it]) 
    }
    SstrMat[j,k] <- tmpjk / dim(outStack$S)[3]
    SstrMat[k,j] <- SstrMat[j,k]
  }
}

corrplot(SstrMat[rownames(c$corr)[100:125],rownames(c$corr)[100:125]], 
         method="circle",order="hclust", 
                   is.corr = F, col = COL1('Blues', 6), 
                   tl.srt = 45,  tl.offset = 0.25, tl.cex=0.6,
                    tl.col = "black",
                   type ="lower", bg = 'white', col.lim=c(3,6),
                   addgrid.col = 'lightgrey') 

#  corrRect(index=c(1,5,10,14,22,27), col="red")



# corrplot::corrplot(SstrMat, order="hclust")

N <- tmp$dims$N
CstrMat <- matrix(0, N,N)
for (i in 1:(N-1)){
  for (k in (i):(N)){
  CstrMat[i,k] <- sum(Cmed[i,]==Cmed[k,])
  CstrMat[k,i] <- CstrMat[i,k]
  }
}

library(corrplot)
# COL2(diverging = c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "RdYlBu"), n = 200)

corrplot::corrplot(method="color", CstrMat[1100:1200,1100:1200], order="FPC", 
                   is.corr = F, col = COL1('YlOrRd', 5))


data.frame(reshape::melt(CstrMat[1100:1200,1100:1200])) %>% 
  ggplot() + 
  geom_tile(aes(x=X2, y=X1, fill=as.factor(value))) + 
  scale_fill_manual(labels=1:6, values=c("white", "yellow", "orange",
                                           "red", "darkred", "black"))
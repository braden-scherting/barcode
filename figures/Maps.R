library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(latex2exp)

source("SetupFunctions.R")
outStack <- stackedList(readRDS("./results/birds6factors05Dec2024.rds"))

tmp <- loadData(noFactors = dim(outStack$C)[2], initialization = F)
dims <- tmp$dims

# D <- plgp::distance(cbind(as.data.frame(models$abu$rL$route$s)$long, 
#                           as.data.frame(models$abu$rL$route$s)$lat))
D <- plgp::distance(tmp$coords)

world <- st_as_sf(maps::map("world", fill=TRUE, plot =FALSE))
FinSf <- world[6,]

# candidateCoords <- expand.grid(
#   seq(min(as.data.frame(models$abu$rL$route$s)$long), 
#       max(as.data.frame(models$abu$rL$route$s)$long), length.out=80/2),
#   seq(min(as.data.frame(models$abu$rL$route$s)$lat), 
#       max(as.data.frame(models$abu$rL$route$s)$lat), length.out=150/2) )

candidateCoords <- expand.grid(
  seq(min(tmp$coords[,1]), 
      max(tmp$coords[,1]), length.out=100/2),
  seq(min(tmp$coords[,2]), 
      max(tmp$coords[,2]), length.out=300/2) )

colnames(candidateCoords) <- c("long", "lat")

candidateCoords_sf <- st_as_sf(candidateCoords, coords = c(1,2))
st_crs(candidateCoords_sf) <- st_crs(FinSf)

in_Finland <- st_intersects(FinSf, candidateCoords_sf)[[1]]
candidateCoords <- candidateCoords[in_Finland,]

obsCoords <- data.frame(long=tmp$coords[,1], 
                        lat=tmp$coords[,2])

Sall <- rbind(obsCoords, candidateCoords)
bigSigma <- (exp(-plgp::distance(Sall) / 0.5) + diag(sqrt(.Machine$double.eps), nrow(Sall)))
S11 <- bigSigma[1:nrow(obsCoords), 1:nrow(obsCoords)]
S22 <- bigSigma[(nrow(obsCoords)+1):nrow(Sall), (nrow(obsCoords)+1):nrow(Sall)]
S12 <- bigSigma[1:nrow(obsCoords), (nrow(obsCoords)+1):nrow(Sall)]
S21 <- t(S12)

niter <- dim(outStack$Xi)[3]
newK <- nrow(candidateCoords)

# Cmat <- matrix(nrow=newK, ncol=dims$L)
# 
# for (l in 1:dims$L){
#   Xil <- (outStack$Xi[,l,]  )
#   # Betal <- (rep(1, newK) %*% t(outStack$Beta[1,l,]))
#     # 
#             # (rep(1, dims$K) %*% t(outStack$Beta[20,l,]))
#   # Cmat[,l] <- rowMeans(pnorm(S21 %*% solve(S11) %*% (Xil + Betal)))
#   Cmat[,l] <- (pnorm(S21 %*% solve(S11) %*% (Xil)))
# }
# 
# LambdaMean <- apply(outStack$Gamma*outStack$S, c(1,2), mean)
# # Especies <- (Cmat*matrix(apply(allPhi, 2, mean, na.rm=T), nrow=newK, ncol=dims$L, byrow = T)) %*% t(LambdaMean)
# Especies <- (Cmat/colSums(Cmat)) %*% t(LambdaMean)

PrCest <- array(1, dim=c(dims$K, dims$L))
for (l in 2:dims$L){
  Xil <- outStack$Xi[,l,]
  XBetal <- t(tmp$dat$des01 / rowSums(tmp$dat$des01)) %*%  tmp$dat$X %*% outStack$Beta[,l,]
  PrCest[,l] <- pnorm(rowMeans(Xil + XBetal))
}

dat <- data.frame(X=obsCoords$long, Y=obsCoords$lat, Z=(PrCest[,2]))
ggplot(dat) + 
  geom_point(aes(x=X, y=Y, fill=Z), size=2.5, shape=21, show.legend = F) + 
  geom_sf(data=world, fill="transparent") +
  coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5) +
  theme_bw() + 
  theme(axis.text.x=element_text(colour = "transparent")) + 
  ylab("") + xlab("") + ggtitle(label=paste("Factor 1"))

for (l in 3:dims$L){
  dat <- data.frame(X=obsCoords$long, Y=obsCoords$lat, Z=(PrCest[,l]))
  
  plt <- ggplot(dat) + 
    geom_point(aes(x=X, y=Y, fill=Z), size=2.5, shape=21, show.legend = F) + 
    geom_sf(data=world, fill="transparent") +
    coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5) +
    theme_bw() + 
    theme(axis.text.x=element_text(colour = "transparent")) + 
    theme(axis.text.y=element_text(colour = "transparent")) + 
    ylab("") + xlab("") + ggtitle(label=paste("Factor ", l-1))
  
  print(plt)
}


XiPred <- array(dim=c(newK, dims$L, niter))
PrC <- array(1, dim=c(newK, dims$L, niter))

for (l in 2:dims$L){
  Xil <- outStack$Xi[,l,]
  Betal <- (rep(1, newK) %*% t(outStack$Beta[1,l,]))
  # 
  # (rep(1, dims$K) %*% t(outStack$Beta[20,l,]))
  # Cmat[,l] <- rowMeans(pnorm(S21 %*% solve(S11) %*% (Xil + Betal)))
  # PrC[,l,] <- ((Betal + S21 %*% solve(S11) %*% (Xil)) > 0) * 1
  PrC[,l,] <- pnorm((Betal + S21 %*% solve(S11) %*% (Xil))) 
}

dat <- data.frame(X=candidateCoords$long, Y=candidateCoords$lat, Z=rowMeans(PrC[,2,]))

ggplot(dat) + 
  geom_point(aes(x=X, y=Y, color=Z), size=3, shape=18, show.legend = F) + 
  geom_sf(data=world, fill="transparent") +
  coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0.5) +
  theme_bw() + 
  ylab("") + xlab("") + ggtitle("")

for (l in 3:dims$L){
  dat <- data.frame(X=candidateCoords$long, Y=candidateCoords$lat, Z=rowMeans(PrC[,l,]))

  plt <- ggplot(dat) + 
    geom_point(aes(x=X, y=Y, color=Z), size=3, shape=18, show.legend = F) + 
    geom_sf(data=world, fill="transparent") +
    coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0.5) +
    theme_bw() + 
    theme(axis.text.y=element_text(colour = "transparent")) + 
    ylab("") + xlab("") + ggtitle("")
  
  print(plt)
}

# PhiOut <- outStack$Phi; PhiOut[outStack$C==0] <- NaN;
# PhiCols <- matrix(rep(colMeans(apply(PhiOut, c(1,2), mean, na.rm=T), na.rm = T) , niter), nrow=dims$L, byrow = F)
# PhiCols <- apply(PhiOut, c(2,3), mean, na.rm=T)

Csum <- apply(PrC, c(2,3), function(x) ifelse(sum(x)==0, Inf, sum(x)) )


specInds <- c(1:3, which(startsWith(colnames(tmp$dat$Y), "Anas")),
              which(startsWith(colnames(tmp$dat$Y), "Cygnus")), 
              which(startsWith(colnames(tmp$dat$Y), "Gallinago")))

specInds <- c(which(startsWith(colnames(tmp$dat$Y), "Anthus")),
              which(startsWith(colnames(tmp$dat$Y), "Sylvia")))

specInds <- c(which(startsWith(colnames(tmp$dat$Y), "Anthus t")),
              which(startsWith(colnames(tmp$dat$Y), "Sylvia bo")),
              which(startsWith(colnames(tmp$dat$Y), "Sylvia cur")))

a <- ggplot() + 
  geom_rect(aes(xmin=0,xmax=6,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect(aes(xmin=0,xmax=1,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect(aes(xmin=1,xmax=2,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect(aes(xmin=2,xmax=3,ymin=0,ymax=1), color="black", fill="darkgrey") + 
  geom_rect(aes(xmin=3,xmax=4,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect(aes(xmin=4,xmax=5,ymin=0,ymax=1), color="black", fill="darkgrey") + 
  geom_rect(aes(xmin=5,xmax=6,ymin=0,ymax=1), color="black", fill="darkgrey") + 
  coord_fixed() + theme_void()

b <- ggplot() + 
  geom_rect(aes(xmin=0,xmax=6,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect(aes(xmin=0,xmax=1,ymin=0,ymax=1), color="black", fill="darkgrey") + 
  geom_rect(aes(xmin=1,xmax=2,ymin=0,ymax=1), color="black", fill="darkgrey") + 
  geom_rect(aes(xmin=2,xmax=3,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect(aes(xmin=3,xmax=4,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect(aes(xmin=4,xmax=5,ymin=0,ymax=1), color="black", fill="darkgrey") + 
  geom_rect(aes(xmin=5,xmax=6,ymin=0,ymax=1), color="black", fill="darkgrey") + 
  coord_fixed() + theme_void()

# c <- ggplot() + 
#   geom_rect(aes(xmin=0,xmax=6,ymin=0,ymax=1), color="black", fill="white") + 
#   geom_rect_pattern(aes(xmin=0,xmax=1,ymin=0,ymax=1), color="black", fill="white") + 
#   geom_rect_pattern(aes(xmin=1,xmax=2,ymin=0,ymax=1), color="black", fill="white") + 
#   geom_rect(aes(xmin=2,xmax=3,ymin=0,ymax=1), color="black", fill="white") + 
#   geom_rect(aes(xmin=3,xmax=4,ymin=0,ymax=1), color="black", fill="white") + 
#   geom_rect_pattern(aes(xmin=4,xmax=5,ymin=0,ymax=1), color="black", fill="white", pattern="stripe") + 
#   geom_rect_pattern(aes(xmin=5,xmax=6,ymin=0,ymax=1), color="black", fill="white", pattern="stripe") + 
#   coord_fixed() + theme_void()

# abc <- c(a,b,c)

for (specInd in specInds){
  Emap <- numeric(newK)
  for (k in 1:newK){
    Emap[k] <- mean(colSums((PrC[k,,]/Csum) * outStack$Gamma[specInd,,] * outStack$S[specInd,,]))
  }
  
  
  # dat <- data.frame(X=candidateCoords$long, Y=candidateCoords$lat, Z=rowMeans(PrC[,2,]))
  dat <- data.frame(X=candidateCoords$long, Y=candidateCoords$lat, Z=Emap)
  
  plt <- ggplot(dat) + 
    geom_point(aes(x=X, y=Y, color=Z), size=3, shape=18) + 
    geom_sf(data=world, fill="transparent") +
    coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=mean(dat$Z)) +
    theme_bw() + 
    labs(color=TeX(r"( $\mu$)")) + 
    ylab("") + xlab("") + 
    ggtitle(paste(colnames(tmp$dat$Y)[specInd]))
  
  print(plt / a)
  print(plt / b)
} 



 ggplot() + 
  geom_rect(aes(xmin=0,xmax=6,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect(aes(xmin=0,xmax=1,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect(aes(xmin=1,xmax=2,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect_pattern(aes(xmin=2,xmax=3,ymin=0,ymax=1), 
                    color="black", fill="white", pattern="stripe") + 
  geom_rect(aes(xmin=3,xmax=4,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect_pattern(aes(xmin=4,xmax=5,ymin=0,ymax=1), color="black", fill="white", pattern="stripe") + 
  geom_rect_pattern(aes(xmin=5,xmax=6,ymin=0,ymax=1), color="black", fill="white", pattern="stripe") + 
  coord_fixed() + theme_void()
  
  # geom_rect(aes(xmin=0,xmax=4,ymin=0,ymax=0.5), color="black", fill="white") + 
  geom_rect_pattern(aes(xmin=0,xmax=4*Cfrac[myl],ymin=0.5,ymax=1), color="black", fill="white", pattern="stripe") +
  geom_rect_pattern(aes(xmin=0,xmax=4*Sfrac[myl],ymin=0,ymax=0.5), color="black", fill="white",pattern="stripe") +
  theme_void() + coord_fixed() +
  xlim(c(-1,4.05)) + 
  annotate("text", x=-0.65, y=0.25, label="Species") +
  annotate("text", x=-0.65, y=0.75, label="Samples")


















Emap <- numeric(newK)
for (k in 1:newK){
  Emap[k] <- mean(colSums((PrC[k,,]/Csum) * outStack$Gamma[specInd,,] * outStack$S[specInd,,]))
}
  

# dat <- data.frame(X=candidateCoords$long, Y=candidateCoords$lat, Z=rowMeans(PrC[,2,]))
dat <- data.frame(X=candidateCoords$long, Y=candidateCoords$lat, Z=Emap)

ggplot(dat) + 
  geom_point(aes(x=X, y=Y, color=Z), size=3, shape=18) + 
  geom_sf(data=world, fill="transparent") +
  coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=mean(dat$Z)) +
  theme_bw() + 
  ylab("") + xlab("") + 
  ggtitle("")

mean(tmp$dat$Y[, specInd])


# specInd <- 2
# dat <- data.frame(X=candidateCoords$long, Y=candidateCoords$lat, Z=Especies[,specInd])

# Cmat <- apply(PrC, c(1,2), mean)
# lind <- 6
# dat <- data.frame(X=candidateCoords$long, Y=candidateCoords$lat, Z=Cmat[,lind])

dat <- data.frame(X=obsCoords$long, Y=obsCoords$lat, Z=rowMeans(pnorm(outStack$Xi[,4,])))

ggplot(dat) + 
  geom_point(aes(x=X, y=Y, color=Z), size=3, shape=18) + 
  geom_sf(data=world, fill="transparent") +
  coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0.75) +
  theme_bw() + 
  ylab("") + xlab("")


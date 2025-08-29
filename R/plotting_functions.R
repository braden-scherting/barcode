library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(patchwork)
library(latex2exp)
library(ggpattern)
library(corrplot)
library(gplots)
library(paletteer)
library(see)

source("setup_functions.R")

plot_factor_maps <- function(data, samples){
  dims <- data$dims
  
  D <- plgp::distance(data$coords)
  
  world <- st_as_sf(maps::map("world", fill=TRUE, plot =FALSE))
  FinSf <- world[6,]
  
  obsCoords <- data.frame(long=data$coords[,1], 
                          lat=data$coords[,2])
  
  c01Phi <- array(1, dim=c(dims$K, dims$L))
  c01Phi[,1] <- (t(data$dat$des01) %*% apply(samples$Phi[,1,], 1, mean)) / colSums(data$dat$des01)
  
  for (l in 2:dims$L){
    Xil <- samples$Xi[,l,]
    XBetal <- t(data$dat$des01) %*%  data$dat$X %*% samples$Beta[,l,]
    cltmp <- (rowMeans(XBetal + Xil) > 0.5) * 1
    c01Phi[,l] <- (t(data$dat$des01) %*% apply(samples$Phi[,l,], 1, mean)) / colSums(data$dat$des01) * cltmp
  }
  
  mapList <- list()
  dat <- data.frame(X=obsCoords$long, Y=obsCoords$lat, Z=(c01Phi[,1]))
  mapList[[1]] <- ggplot(dat) +
    geom_point(aes(x=X, y=Y, fill=log(Z)), size=2., shape=21, show.legend = T ) +
    geom_sf(data=world, fill="transparent") +
    coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
    scale_fill_gradient2(low="white", mid="pink", high="darkred", midpoint=-8, limits=c(-9.25,-4.9)) +
    # scale_color_manual(values=c("black", "black")) +
    theme_bw() +
    theme(axis.text.x=element_text(colour = "transparent"), legend.title = element_blank()) +
    ylab("") + xlab("") + ggtitle(label=paste0("Factor 1/",dims$L," (Fixed)"))
  
  for (l in 2:dims$L){
    dat <- data.frame(X=obsCoords$long, Y=obsCoords$lat, Z=(c01Phi[,l]))
    
    mapList[[l]] <- ggplot() + 
      geom_point(data=dat[dat$Z==0,], mapping=aes(x=X, y=Y), color="lightgray", fill="white", size=1.25, shape=22, show.legend = F) + 
      geom_point(data=dat[dat$Z>0,], mapping=aes(x=X, y=Y, fill=log(Z)), color="black", size=1.25, shape=22, show.legend = T) + 
      geom_sf(data=world, fill="transparent") +
      coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
      scale_fill_gradient2(low="white", high="darkred", midpoint=-8, limits=c(-9.25,-4.9)) +
      theme_bw() + 
      labs(title = paste0("Factor ", l,"/",dims$L)) + 
      theme(legend.title = element_blank(), plot.title = element_text(family = "serif"),
            legend.text = element_text(family = "serif"), 
            axis.text.y=element_text(colour = "transparent"),
            axis.text.x=element_text(colour = "transparent")) +
      ylab("") + xlab("")
    
    # if (l != 2){
    #   mapList[[l]] <- mapList[[l]] + theme(axis.text.y=element_text(colour = "transparent"))
    # }
    
    # if (l < 4){ 
    # mapList[[l]] <- mapList[[l]] + theme(axis.text.x=element_text(colour = "transparent"))  
    # }
  }
  return(mapList)
}

plot_factor_proportions <- function(data, samples){
  year <- data$dat$year
  
  year01 <- model.matrix(~factor(year)-1)
  CPhiMean <- apply(samples$C * samples$Phi, c(1,2), mean)
  
  yCPhi <- t(year01) %*% CPhiMean[,2:data$dims$L]
  traj <- yCPhi / matrix(rep(rowSums(yCPhi),data$dims$L-1), nrow=11)
  rownames(traj) <- paste(2006:2016)
  colnames(traj) <- paste(2:data$dims$L)
  
  plt <- reshape::melt(traj) %>% 
        ggplot() + 
        geom_bar(aes(x=as.factor(X1), y=value, fill=as.factor(X2)), stat = "identity", color="black") + 
        scale_fill_okabeito(name="Factor", order=2:data$dims$L) + 
        theme_minimal() + ylab("Factor relevance") + xlab("") +
        theme(axis.text = element_text(family="serif"),
              axis.title = element_text(family="serif"))
  return(plt)
}

plot_covariate_effects <- function(data, samples){
  EB <- apply(samples$Beta, c(1,2), mean)
  EL <- apply(samples$Beta, c(1,2), quantile, probs=0.025)
  EU <- apply(samples$Beta, c(1,2), quantile, probs=0.975)
  
  ES <- sign(EL) == sign(EU)
  
  dirMat <- sign(EB[,2:data$dims$L])*ES[,2:data$dims$L]
  magMat <- sign(EB[,2:data$dims$L])*((abs(EB[,2:data$dims$L])^0.5))*ES[,2:data$dims$L]
  
  xlabs <- paste(colnames(data$dat$X))
  xlabs <- stringr::str_replace_all(xlabs, "_", " ")
  xlabs <- stringr::str_replace_all(xlabs, "poly", "")
  xlabs <- stringr::str_replace_all(xlabs, "1", "")
  rownames(magMat) <- xlabs
  
  plt <- reshape::melt(magMat) %>% 
          ggplot() + 
          geom_tile(aes(y=as.factor(X2+1), x=X1, fill=value), color="black", size=0.1, show.legend = T) + 
          # scale_fill_viridis(option="Magma") + 
          scale_fill_gradient2(low=viridis::plasma(10)[2], mid="white", high=viridis::plasma(10)[8], 
                               midpoint=0, name="", limits=c(-1.5,2)) + 
          geom_point(data=data.frame(x=1:(data$dims$L-1), y=apply(magMat[2:21,], 2, which.max)+1),
                     mapping=aes(x=y, y=x), color="black", shape="+", size=5) +
          geom_point(data=data.frame(x=1:(data$dims$L-1), y=apply(magMat[2:21,], 2, which.min)+1),
                     mapping=aes(x=y, y=x), color="black", shape="-", size=5) +
          # scale_color_manual(values=c("transparent", "black"), guide="none") +
          theme_classic() + ylab("") + xlab("") +
          theme(axis.text.x=element_text(angle = 70, hjust=0.95, vjust=0.95, family="serif"),
                legend.text=element_text(family="serif"),
                axis.text.y=element_text(family="serif"))+
          coord_fixed()
  return(plt)
}

plot_clustered_loadings <- function(data, samples, retClust=F){
  Smed <- apply(samples$S, c(1,2), median)
  SGmed <- apply(samples$S * samples$Gamma, c(1,2), median)
  
  # row.names(SGmed) <- colnames(data$dat$Y)
  row.names(Smed) <- colnames(data$dat$Y)
  
  pmv <- matrix(NA, data$dims$P, 1000)
  pme <- matrix(NA, data$dims$P, 1000)
  for (it in 1:1000){
    E <- (samples$Phi[,,it*10] * samples$C[,,it*10]) %*% t(samples$Gamma[,,it*10] * samples$S[,,it*10])
    Ypred <- rpois(data$dims$N * data$dims$P, c(E)) %>% matrix(.,nrow=data$dims$N, ncol=data$dims$P)  
    pmv[,it] <- apply(Ypred, 2, var)
    pme[,it] <- apply(Ypred, 2, mean)
  } 
  
  empMarVar <- apply(data$dat$Y, 2, var)
  empMarE <- apply(data$dat$Y, 2, mean)
  
  height <- apply(pmv, 1, mean) / empMarVar
  heightE <- apply(pme, 1, mean) / empMarE
  
  hm2 <- heatmap.2(Smed, Colv = NA,  dendrogram = "row", key=F, margins=c(7.5,10),
                   keysize = 0.1, trace = "none", colsep=0:11, rowsep=1:10, sepcolor = "black", 
                   sepwidth = 0.01, cexRow=0.5, lmat = rbind(4:3,2:1), lwid = c(.25,1), cexCol = 1.25)
  
  gooeySG <- reshape::melt(
    (SGmed / matrix(rep(rowSums(SGmed), data$dims$L), nrow=data$dims$P, ncol=data$dims$L, byrow = F))[hm2$rowInd,]
  )
  gooeySG$segment_height <- gooeySG$value * height
  
  splitNames <- str_split(colnames(data$dat$Y)[hm2$rowInd], " ", simplify = T)
  G <- str_sub(splitNames[,1], 1, 1)
  namesShort <- paste0("italic('", str_glue("{G}. {splitNames[,2]}"), "')")
  
  plt <- ggplot() +
    geom_segment(aes(x=1:data$dims$P, xend=1:data$dims$P, y=rep(0, data$dims$P), 
                     yend=c(rep(c(1.02, 1.145), 66))), color="gray") +
    geom_col(data=gooeySG, mapping=aes(x = X1, y = segment_height, fill = as.factor(X2)), 
             color="black", size=0.2) + 
    geom_segment(aes(x=1, xend=data$dims$P, y=1, yend=1), color="black") +
    theme_void() + ylab("") + xlab("") +
    theme(axis.text.x = element_blank(),
          legend.position = c(0.5,0.07), legend.direction = "horizontal") +
    scale_fill_okabeito(name="Factor", order=1:data$dims$L) + 
    labs(fill = "Factor") +
    ylim(0,1.3) + xlim(-10,138) +
    annotate("text", x=1:data$dims$P, y=c(rep(c(1.025, 1.15), 66)), label=namesShort, 
             parse = TRUE, hjust = 0., size=1.7, angle=c(rep(c(0, 0), 66))) +
    coord_flip()
  
  gooeySG$segment_height <- gooeySG$value * heightE
  pltE <- ggplot() +
    geom_segment(aes(x=1:data$dims$P, xend=1:data$dims$P, y=rep(0, data$dims$P), 
                     yend=c(rep(c(1.02, 1.095), 66))), color="gray") +
    geom_col(data=gooeySG, mapping=aes(x = X1, y = segment_height, fill = as.factor(X2)), 
             color="black", size=0.2) + 
    geom_segment(aes(x=1, xend=data$dims$P, y=1, yend=1), color="black") +
    theme_void() + ylab("") + xlab("") +
    theme(axis.text.x = element_blank(),
          legend.position = c(0.5,0.07), legend.direction = "horizontal") +
    scale_fill_okabeito(name="Factor", order=1:data$dims$L) + 
    labs(fill = "Factor") +
    ylim(0,1.3) + xlim(-10,138) +
    annotate("text", x=1:data$dims$P-0.2, y=c(rep(c(1.025, 1.2), 66)), label=namesShort, 
             parse = TRUE, hjust = 0., size=1.7, angle=c(rep(c(0, 0), 66))) +
    coord_flip()
  if (retClust) {
    return(list(plt, pltE, hm2))
  } else{
    return(list(plt, pltE))
  }
}

plot_species_correlations <- function(data, samples, clust=NULL){
  GamS <- samples$Gamma*samples$S
  PhiC <- samples$Phi*samples$C
  niter <- dim(samples$Gamma)[3]
  
  postCov <- matrix(0, dim(samples$Gamma)[1], dim(samples$Gamma)[1])
  for (it in 1:(niter)){
    EV <- diag(c((samples$Gamma[,,it] * samples$S[,,it]) %*% colMeans(samples$Phi[,,it] * samples$C[,,it])))
    VE <- (samples$Gamma[,,it] * samples$S[,,it]) %*% var(samples$Phi[,,it] * samples$C[,,it])%*% t(samples$Gamma[,,it] * samples$S[,,it])
    
    # E <- (samples$Phi[,,it*10] * samples$C[,,it*10]) %*% t(samples$Gamma[,,it*10] * samples$S[,,it*10])
    # Ypred <- rpois(data$dims$N * data$dims$P, c(E)) %>% matrix(.,nrow=data$dims$N, ncol=data$dims$P)  
    # postCov <- postCov + var(Ypred)
    postCov <- postCov + (EV + VE)
  }
  postCorr <- cov2cor(postCov / (niter))
  colnames(postCorr) <- colnames(data$dat$Y)
  rownames(postCorr) <- colnames(data$dat$Y)
  
  if (is.null(clust)){
    df <- reshape2::melt(postCorr) 
  } else {
    df <- reshape2::melt(postCorr[rev(clust), c(clust)]) 
  }
  
  plt <- ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    # geom_text(aes(label = round(value, 2)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-1, 1)) +
    theme_minimal() +
    coord_fixed() +
    labs(x = "", y = "", fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank())
  # theme(axis.text.y = element_text(size = 5)) +
  # theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1))
  return(plt)
}

plot_obs_correlations <- function(data,  clust=NULL){
  postCorr <- cor(data$dat$Y)
  colnames(postCorr) <- colnames(data$dat$Y)
  rownames(postCorr) <- colnames(data$dat$Y)
  
  if (is.null(clust)){
    df <- reshape2::melt(postCorr) 
  } else {
    df <- reshape2::melt(postCorr[rev(clust), c(clust)]) 
  }
  
  plt <- ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    # geom_text(aes(label = round(value, 2)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-1, 1)) +
    theme_minimal() +
    coord_fixed() +
    labs(x = "", y = "", fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank())
    # theme(axis.text.y = element_text(size = 5)) +
    # theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1))
  return(plt)
}

samples5 <- stackedList(readRDS("../output/results/rank5_xxx"))
mydata5 <- loadData(noFactors = ncol(samples5$Gamma))

samples6 <- stackedList(readRDS("../output/results/rank6_xxx"))
mydata6 <- loadData(noFactors = ncol(samples6$Gamma))

samples7 <- stackedList(readRDS("../output/results/rank7_xxx"))
mydata7 <- loadData(noFactors = ncol(samples7$Gamma))

samples8 <- stackedList(readRDS("../output/results/rank8_xxx"))
mydata8 <- loadData(noFactors = ncol(samples8$Gamma))

maps5 <- plot_factor_maps(mydata5, samples5)
maps6 <- plot_factor_maps(mydata6, samples6)
maps7 <- plot_factor_maps(mydata7, samples7)
maps8 <- plot_factor_maps(mydata8, samples8)

wrap_plots(c(maps5[c(2,3,5,4)], maps6[c(2,6,4,5,3)], maps7[c(6,3,7,5,2,4)],
             maps8[c(4,3,2,7,5,6,8)])) + 
  plot_layout(design=c("ABCD###
                        EFGHI##
                        JKLMNO#
                        PQRSTUV"), guides = "collect") & 
  theme(legend.position.inside = c(0.5,0.5), 
        plot.margin=margin(0,0,0,0), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        legend.text = element_text(size=15))

plot_covariate_effects(mydata5, samples5)
plot_covariate_effects(mydata6, samples6)
cov7 <- plot_covariate_effects(mydata7, samples7)
time7 <- plot_factor_proportions(mydata7, samples7)

pcl <- plot_clustered_loadings(mydata7, samples7, T)
plot_obs_correlations(mydata7, pcl[[3]]$rowInd)
plot_species_correlations(mydata7, samples7, pcl[[3]]$rowInd)

cov7 & theme(legend.position = "right",
             axis.text=element_text(size=15),
             legend.text=element_text(size=12))

time7 & theme(legend.position = "right",
              text=element_text(size=15, family="serif"),
              legend.text=element_text(size=12))

wrap_plots(c(maps7[c(2:7)])) + 
  plot_layout(design=c("123456"), guides = "collect", 
              heights = c(1,1)) & 
  theme(legend.position = "right", axis.text=element_text(size=12), plot.margin = margin(2,0,0,0)) 

# plot_clustered_loadings(mydata6, samples6)
load7 <- plot_clustered_loadings(mydata7, samples7)
load7[[1]] + guides(fill = guide_legend(nrow = 1)) & theme(legend.text = element_text(family="serif"),
                                                           legend.title = element_text(family="serif"))
#--------------------------#
# Re-ordered factors
#--------------------------#

samples <- stackedList(readRDS("../output/results/rank7_xxx"))
data <- loadData(noFactors = ncol(samples$Gamma))

ordered7factorMaps <- function(data, samples){
  dims <- data$dims
  
  D <- plgp::distance(data$coords)
  
  world <- st_as_sf(maps::map("world", fill=TRUE, plot =FALSE))
  FinSf <- world[6,]
  
  obsCoords <- data.frame(long=data$coords[,1], 
                          lat=data$coords[,2])
  
  c01Phi <- array(1, dim=c(dims$K, dims$L))
  c01Phi[,1] <- (t(data$dat$des01) %*% apply(samples$Phi[,1,], 1, mean)) / colSums(data$dat$des01)
  
  for (l in 2:dims$L){
    Xil <- samples$Xi[,l,]
    XBetal <- t(data$dat$des01) %*%  data$dat$X %*% samples$Beta[,l,]
    cltmp <- (rowMeans(XBetal + Xil) > 0.5) * 1
    c01Phi[,l] <- (t(data$dat$des01) %*% apply(samples$Phi[,l,], 1, mean)) / colSums(data$dat$des01) * cltmp
  }
  
  factorNames <- c("(1) Ref.", "(2) Urban", "(6) N. Old Growth", "(7) Fjell + wetland", 
                   "(3) Agriculture", "(4) S. Mixed Forest", "(5) Pine")
  
  mapList <- list()
  dat <- data.frame(X=obsCoords$long, Y=obsCoords$lat, Z=(c01Phi[,1]))
  mapList[[1]] <- ggplot(dat) +
    geom_point(aes(x=X, y=Y, fill=log(Z)), size=2., shape=21, show.legend = T ) +
    geom_sf(data=world, fill="transparent") +
    coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
    scale_fill_gradient2(low="white", mid="pink", high="darkred", midpoint=-8, limits=c(-9.25,-4.9)) +
    # scale_color_manual(values=c("black", "black")) +
    theme_bw() +
    theme(axis.text.x=element_text(colour = "transparent"), legend.title = element_blank()) +
    ylab("") + xlab("") + ggtitle(label=factorNames[1])
  
  for (l in 2:dims$L){
    dat <- data.frame(X=obsCoords$long, Y=obsCoords$lat, Z=(c01Phi[,l]))
    
    mapList[[l]] <- ggplot() + 
      geom_point(data=dat[dat$Z==0,], mapping=aes(x=X, y=Y), color="lightgray", fill="white", size=1.75, shape=22, show.legend = F) + 
      geom_point(data=dat[dat$Z>0,], mapping=aes(x=X, y=Y, fill=log(Z)), color="black", size=1.75, shape=22, show.legend = T) + 
      geom_sf(data=world, fill="transparent") +
      coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
      scale_fill_gradient2(low="white", high="darkred", midpoint=-8, limits=c(-9.25,-4.9)) +
      theme_bw() + 
      labs(title = factorNames[l]) + 
      theme(legend.title = element_blank(), plot.title = element_text(family = "serif"),
            legend.text = element_text(family = "serif"), 
            axis.text.y=element_text(colour = "transparent"),
            axis.text.x=element_text(colour = "transparent")) +
      ylab("") + xlab("")
    
    # if (l != 2){
    #   mapList[[l]] <- mapList[[l]] + theme(axis.text.y=element_text(colour = "transparent"))
    # }
    
    # if (l < 4){ 
    # mapList[[l]] <- mapList[[l]] + theme(axis.text.x=element_text(colour = "transparent"))  
    # }
  }
  return(mapList)
}
ordered7covariates <- function(data, samples){
  EB <- apply(samples$Beta, c(1,2), mean)
  EL <- apply(samples$Beta, c(1,2), quantile, probs=0.025)
  EU <- apply(samples$Beta, c(1,2), quantile, probs=0.975)
  
  ES <- sign(EL) == sign(EU)
  
  dirMat <- sign(EB[,2:data$dims$L])*ES[,2:data$dims$L]
  magMat <- sign(EB[,2:data$dims$L])*((abs(EB[,2:data$dims$L])^0.5))*ES[,2:data$dims$L]
  
  xlabs <- paste(colnames(data$dat$X))
  xlabs <- stringr::str_replace_all(xlabs, "_", " ")
  xlabs <- stringr::str_replace_all(xlabs, "poly", "")
  xlabs <- stringr::str_replace_all(xlabs, "1", "")
  rownames(magMat) <- xlabs
  magMat <- magMat[,c(1,4,5,6,2,3)]
  colnames(magMat) <- c("(2) Urban", "(3) Agriculture", "(4) S. Mixed Forest", "(5) Pine", "(6) N. Old Growth", "(7) Fjell + wetland")
  
  plt <- reshape::melt(magMat) %>% 
    ggplot() + 
    geom_tile(aes(y=as.factor(X2), x=X1, fill=value), color="black", size=0.1, show.legend = T) + 
    # scale_fill_viridis(option="Magma") + 
    scale_fill_gradient2(low=viridis::plasma(10)[2], mid="white", high=viridis::plasma(10)[8], 
                         midpoint=0, name="", limits=c(-1.5,2)) + 
    geom_point(data=data.frame(x=1:(data$dims$L-1), y=apply(magMat[2:21,], 2, which.max)+1),
               mapping=aes(x=y, y=x), color="black", shape="+", size=5) +
    geom_point(data=data.frame(x=1:(data$dims$L-1), y=apply(magMat[2:21,], 2, which.min)+1),
               mapping=aes(x=y, y=x), color="black", shape="-", size=5) +
    # scale_color_manual(values=c("transparent", "black"), guide="none") +
    theme_classic() + ylab("") + xlab("") +
    theme(axis.text.x=element_text(angle = 70, hjust=0.95, vjust=0.95, family="serif"),
          legend.text=element_text(family="serif"),
          axis.text.y=element_text(family="serif"))+
    coord_fixed()
  return(plt)
}
ordered7factorProportions <- function(data, samples){
  year <- data$dat$year
  
  year01 <- model.matrix(~factor(year)-1)
  CPhiMean <- apply(samples$C * samples$Phi, c(1,2), mean)
  
  yCPhi <- t(year01) %*% CPhiMean[,2:data$dims$L]
  traj <- yCPhi / matrix(rep(rowSums(yCPhi),data$dims$L-1), nrow=11)
  rownames(traj) <- paste(2006:2016)
  # colnames(traj) <- paste(2:data$dims$L)
  traj <- traj[,c(3,2,6,5,4,1)]
  colnames(traj) <- c("(2) Urban", "(3) Agriculture", "(4) S. Mixed Forest", "(5) Pine", "(6) N. Old Growth", "(7) Fjell + wetland")[6:1]
  
  plt <- reshape::melt(traj) %>% 
    ggplot() + 
    geom_bar(aes(x=as.factor(X1), y=value, fill=as.factor(X2)), stat = "identity", color="black") + 
    scale_fill_okabeito(name="Factor", order=data$dims$L:2) + 
    theme_minimal() + ylab("Factor relevance") + xlab("") +
    theme(axis.text = element_text(family="serif"),
          axis.title = element_text(family="serif"),
          legend.title = element_blank(),
          legend.text = element_text(family="serif"))
  return(plt)
}

ordered7clusterLoadings <- function(data, samples, retClust=F){
  Smed <- apply(samples$S, c(1,2), median)
  SGmed <- apply(samples$S * samples$Gamma, c(1,2), median)
  
  # row.names(SGmed) <- colnames(data$dat$Y)
  row.names(Smed) <- colnames(data$dat$Y)
  
  pmv <- matrix(NA, data$dims$P, 1000)
  pme <- matrix(NA, data$dims$P, 1000)
  for (it in 1:1000){
    E <- (samples$Phi[,,it*10] * samples$C[,,it*10]) %*% t(samples$Gamma[,,it*10] * samples$S[,,it*10])
    Ypred <- rpois(data$dims$N * data$dims$P, c(E)) %>% matrix(.,nrow=data$dims$N, ncol=data$dims$P)  
    pmv[,it] <- apply(Ypred, 2, var)
    pme[,it] <- apply(Ypred, 2, mean)
  } 
  
  empMarVar <- apply(data$dat$Y, 2, var)
  empMarE <- apply(data$dat$Y, 2, mean)
  
  height <- apply(pmv, 1, mean) / empMarVar
  heightE <- apply(pme, 1, mean) / empMarE
  
  hm2 <- heatmap.2(Smed, Colv = NA,  dendrogram = "row", key=F, margins=c(7.5,10),
                   keysize = 0.1, trace = "none", colsep=0:11, rowsep=1:10, sepcolor = "black", 
                   sepwidth = 0.01, cexRow=0.5, lmat = rbind(4:3,2:1), lwid = c(.25,1), cexCol = 1.25)
  
  meltThis <- (SGmed / matrix(rep(rowSums(SGmed), data$dims$L), nrow=data$dims$P, 
                              ncol=data$dims$L, byrow = F))[hm2$rowInd,c(1,2,5,6,7,3,4)]
  colnames(meltThis) <- c("(1) Ref.", "(2) Urban", "(3) Agriculture", "(4) S. Mixed Forest", 
                          "(5) Pine", "(6) N. Old Growth", "(7) Fjell + wetland")
  gooeySG <- reshape::melt(meltThis)
  gooeySG$segment_height <- gooeySG$value * height
  
  splitNames <- str_split(colnames(data$dat$Y)[hm2$rowInd], " ", simplify = T)
  G <- str_sub(splitNames[,1], 1, 1)
  namesShort <- paste0("italic('", str_glue("{G}. {splitNames[,2]}"), "')")
  
  plt <- ggplot() +
    geom_segment(aes(x=1:data$dims$P, xend=1:data$dims$P, y=rep(0, data$dims$P), 
                     yend=c(rep(c(1.02, 1.145), 66))), color="gray") +
    geom_col(data=gooeySG, mapping=aes(x = X1, y = segment_height, fill = as.factor(X2)), 
             color="black", size=0.2) + 
    geom_segment(aes(x=1, xend=data$dims$P, y=1, yend=1), color="black") +
    theme_void() + ylab("") + xlab("") +
    theme(axis.text.x = element_blank(), legend.text = element_text(family="serif"),
          legend.position = "bottom", legend.direction = "vertical", legend.byrow = T) +
    scale_fill_okabeito(name="", order=1:data$dims$L) + 
    labs(fill = "") +
    ylim(0,1.3) + xlim(-10,138) +
    annotate("text", x=1:data$dims$P, y=c(rep(c(1.025, 1.15), 66)), label=namesShort, 
             parse = TRUE, hjust = 0., size=1.7, angle=c(rep(c(0, 0), 66))) +
    coord_flip()
  
  gooeySG$segment_height <- gooeySG$value * heightE
  ggplot() +
    geom_segment(aes(x=1:data$dims$P, xend=1:data$dims$P, y=rep(0, data$dims$P), 
                     yend=c(rep(c(1.04, 1.15), 66))), color="gray") +
    geom_col(data=gooeySG, mapping=aes(x = X1, y = segment_height, fill = as.factor(X2)), 
             color="black", size=0.2, show.legend = F) + 
    geom_segment(aes(x=1, xend=data$dims$P, y=1, yend=1), color="black") +
    theme_void() + ylab("") + xlab("") +
    theme(axis.text.x = element_blank(), legend.text = element_text(family="serif"),
          legend.position = "bottom", legend.direction = "vertical", legend.byrow = T) +
    scale_fill_okabeito(name="", order=1:data$dims$L) + 
    labs(fill = "") +
    ylim(0,1.3) + xlim(-10,138) +
    annotate("text", x=1:data$dims$P, y=c(rep(c(1.045, 1.15), 66)), label=namesShort, 
             parse = TRUE, hjust = 0., size=1.7, angle=c(rep(c(0, 0), 66))) +
    coord_flip()
  if (retClust) {
    return(list(plt, pltE, hm2))
  } else{
    return(list(plt, pltE))
  }
}

dominantFactor <- function(data, samples){
  dims <- data$dims
  
  world <- st_as_sf(maps::map("world", fill=TRUE, plot =FALSE))
  FinSf <- world[6,]
  
  obsCoords <- data.frame(long=data$coords[,1], 
                          lat=data$coords[,2])
  
  domFactor <- apply(t(data$dat$des01) %*% apply(samples$C[,2:7,]*samples$Phi[,2:7,], c(1,2), mean), 
                c(1), which.max)
  newlevels <- c("(2) Urban", "(6) N. Old Growth", "(7) Fjell + wetland", 
                 "(3) Agriculture", "(4) S. Mixed Forest", "(5) Pine")
  
  
  dat <- data.frame(X=obsCoords$long, Y=obsCoords$lat, Z=newlevels[as.factor(domFactor)])
  plt <- ggplot(dat) +
    geom_sf(data=world, fill="lightgray", alpha=0.5) +
    geom_point(aes(x=X, y=Y, fill=(Z)), color="transparent", size=4, shape=21, show.legend = T) +
    coord_sf(xlim = c(18.5, 31.5), ylim = c(59.5, 70.5), expand = FALSE) +
    # scale_fill_gradient2(low="white", mid="pink", high="darkred", midpoint=-8, limits=c(-9.25,-4.9)) +
    # scale_color_manual(values=c("black", "black")) +
    scale_fill_okabeito(name="", order=2:data$dims$L) + 
    theme_bw() +
    theme(axis.text.x=element_text(colour = "transparent"), legend.title = element_blank(),
          legend.text = element_text(family = "serif"), axis.text.y = element_blank()) +
    ylab("") + xlab("") 
  return(plt)
}

maps7.2 <- ordered7factorMaps(data, samples)
covar7.2 <- ordered7covariates(data, samples)
props7.2 <- ordered7factorProportions(data, samples)
clust7.2 <- ordered7clusterLoadings(data, samples)[[1]]
  
wrap_plots(c(maps7.2[c(2:7)])) + 
  plot_layout(design=c("145623"), guides = "collect", 
              heights = c(1,1)) & 
  theme(legend.position = "right", axis.text=element_text(size=12), plot.margin = margin(2,0,0,0)) 

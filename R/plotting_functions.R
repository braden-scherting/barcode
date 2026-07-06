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


# Function for producing maps displayed in Figure 1
plot_factor_maps <- function(data, samples){
  dims <- data$dims
  
  D <- plgp::distance(data$coords)
  
  # Access the appropriate map outline
  world <- st_as_sf(maps::map("world", fill=TRUE, plot =FALSE))
  FinSf <- world[6,]
  
  obsCoords <- data.frame(long=data$coords[,1], 
                          lat=data$coords[,2])
  
  # Compute tile fill/outline based on posterior summaries
  c01Phi <- array(1, dim=c(dims$K, dims$L))
  c01Phi[,1] <- (t(data$dat$des01) %*% apply(samples$Phi[,1,], 1, mean)) / colSums(data$dat$des01)
  
  for (l in 2:dims$L){
    Xil <- samples$Xi[,l,]
    XBetal <- t(data$dat$des01) %*%  data$dat$X %*% samples$Beta[,l,]
    cltmp <- (rowMeans(XBetal + Xil) > 0.5) * 1
    c01Phi[,l] <- (t(data$dat$des01) %*% apply(samples$Phi[,l,], 1, mean)) / colSums(data$dat$des01) * cltmp
  }
  
  mapList <- list()
  
  # Because c_{i1}=1 for all, i, all tiles are outlined in the first map
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
  
  # Create and store the remaining maps, coloring the tile outline by sample factor presence 0/1
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

  }
  return(mapList)
}

# Code for reproducing Figure 3, a stacked and normalized bar chart
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

# Code for reproducing Figure 2, latent covariate effects
plot_covariate_effects <- function(data, samples){
  # Posterior summaries for evaluating statistical support
  EB <- apply(samples$Beta, c(1,2), mean)
  EL <- apply(samples$Beta, c(1,2), quantile, probs=0.025)
  EU <- apply(samples$Beta, c(1,2), quantile, probs=0.975)
  
  ES <- sign(EL) == sign(EU)
  
  dirMat <- sign(EB[,2:data$dims$L])*ES[,2:data$dims$L]
  magMat <- sign(EB[,2:data$dims$L])*((abs(EB[,2:data$dims$L])^0.5))*ES[,2:data$dims$L]
  
  # Clean up variable names
  xlabs <- paste(colnames(data$dat$X))
  xlabs <- stringr::str_replace_all(xlabs, "_", " ")
  xlabs <- stringr::str_replace_all(xlabs, "poly", "")
  xlabs <- stringr::str_replace_all(xlabs, "1", "")
  rownames(magMat) <- xlabs
  
  plt <- reshape::melt(magMat) %>% 
          ggplot() + 
          geom_tile(aes(y=as.factor(X2+1), x=X1, fill=value), color="black", size=0.1, show.legend = T) + 
          scale_fill_gradient2(low=viridis::plasma(10)[2], mid="white", high=viridis::plasma(10)[8], 
                               midpoint=0, name="", limits=c(-1.5,2)) + 
          # Markers for strongest positive effects
          geom_point(data=data.frame(x=1:(data$dims$L-1), y=apply(magMat[2:21,], 2, which.max)+1),
                     mapping=aes(x=y, y=x), color="black", shape="+", size=5) +
          # Markers for strongest negative effects
          geom_point(data=data.frame(x=1:(data$dims$L-1), y=apply(magMat[2:21,], 2, which.min)+1),
                     mapping=aes(x=y, y=x), color="black", shape="-", size=5) +
          theme_classic() + ylab("") + xlab("") +
          theme(axis.text.x=element_text(angle = 70, hjust=0.95, vjust=0.95, family="serif"),
                legend.text=element_text(family="serif"),
                axis.text.y=element_text(family="serif"))+
          coord_fixed()
  return(plt)
}

# Code for reproducing Figures 4 and A.2
plot_clustered_loadings <- function(data, samples, retClust=F){
  # Posterior summaries of species-specific factors
  Smed <- apply(samples$S, c(1,2), median)
  SGmed <- apply(samples$S * samples$Gamma, c(1,2), median)
  
  row.names(Smed) <- colnames(data$dat$Y)
  
  # Marginal posterior summaries
  pmv <- matrix(NA, data$dims$P, 1000)
  pme <- matrix(NA, data$dims$P, 1000)
  for (it in 1:1000){
    E <- (samples$Phi[,,it*10] * samples$C[,,it*10]) %*% t(samples$Gamma[,,it*10] * samples$S[,,it*10])
    Ypred <- rpois(data$dims$N * data$dims$P, c(E)) %>% matrix(.,nrow=data$dims$N, ncol=data$dims$P)  
    pmv[,it] <- apply(Ypred, 2, var)
    pme[,it] <- apply(Ypred, 2, mean)
  } 
  
  # Marginal empirical summaries
  empMarVar <- apply(data$dat$Y, 2, var)
  empMarE <- apply(data$dat$Y, 2, mean)
  
  # Ratios of posterior to empirical marginal mean and variance
  height <- apply(pmv, 1, mean) / empMarVar
  heightE <- apply(pme, 1, mean) / empMarE
  
  # Create heatmap to extract dendrogram & hierarchical clustering ordering based on 0/1 in S
  hm2 <- heatmap.2(Smed, Colv = NA,  dendrogram = "row", key=F, margins=c(7.5,10),
                   keysize = 0.1, trace = "none", colsep=0:11, rowsep=1:10, sepcolor = "black", 
                   sepwidth = 0.01, cexRow=0.5, lmat = rbind(4:3,2:1), lwid = c(.25,1), cexCol = 1.25)
  # Scale bar width by factor contribution and reorder per dendrogram
  gooeySG <- reshape::melt(
    (SGmed / matrix(rep(rowSums(SGmed), data$dims$L), nrow=data$dims$P, ncol=data$dims$L, byrow = F))[hm2$rowInd,]
  )
  # Scale bar width by posterior/empirical ratio
  gooeySG$segment_height <- gooeySG$value * height
  
  # Clean up species names
  splitNames <- str_split(colnames(data$dat$Y)[hm2$rowInd], " ", simplify = T)
  G <- str_sub(splitNames[,1], 1, 1)
  namesShort <- paste0("italic('", str_glue("{G}. {splitNames[,2]}"), "')")
  
  # Marginal variance: Figure 4
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
  # Marginal expectation: Figure A.2
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

# Code for reproducing Figure A.3 (left)
plot_species_correlations <- function(data, samples, clust=NULL){
  GamS <- samples$Gamma*samples$S
  PhiC <- samples$Phi*samples$C
  niter <- dim(samples$Gamma)[3]
  
  # Compute posterior marginal species correlation
  postCov <- matrix(0, dim(samples$Gamma)[1], dim(samples$Gamma)[1])
  for (it in 1:(niter)){
    EV <- diag(c((samples$Gamma[,,it] * samples$S[,,it]) %*% colMeans(samples$Phi[,,it] * samples$C[,,it])))
    VE <- (samples$Gamma[,,it] * samples$S[,,it]) %*% var(samples$Phi[,,it] * samples$C[,,it])%*% t(samples$Gamma[,,it] * samples$S[,,it])
    postCov <- postCov + (EV + VE)
  }
  postCorr <- cov2cor(postCov / (niter))
  colnames(postCorr) <- colnames(data$dat$Y)
  rownames(postCorr) <- colnames(data$dat$Y)
  
  # Optionally rearrange rows and columns based on species clusters
  if (is.null(clust)){
    df <- reshape2::melt(postCorr) 
  } else {
    df <- reshape2::melt(postCorr[rev(clust), c(clust)]) 
  }
  
  plt <- ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-1, 1)) +
    theme_minimal() +
    coord_fixed() +
    labs(x = "", y = "", fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank())
  return(plt)
}

# Code for reproducing Figure A.3 (right)
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
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-1, 1)) +
    theme_minimal() +
    coord_fixed() +
    labs(x = "", y = "", fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.y = element_blank()) +
    theme(axis.text.x = element_blank())
  return(plt)
}

# Read output based on model fits in `output/results/` and data to match
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

# Arranging plots to create Figure A.1
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

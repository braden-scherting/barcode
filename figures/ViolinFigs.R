library(patchwork)
library(ggplot2)
library(ggpattern)

RGC <- readRDS("results/RoundGaussConcurrent.rds")
PCU <- readRDS("results/PoissonCopulaUnconstrained.rds")


data.frame(RMSE=c(PCU$pCopRMSE), method=c(rep("Barcode", 100),
                                              rep("Poisson GLLVM", 100))) %>% 
ggplot() + 
  geom_boxplot(aes(y=log(RMSE), x=method)) + theme(aspect.ratio=1)

data.frame(MAE=c(PCU$pCopMAE), method=c(rep("Barcode", 100),
                                              rep("Poisson GLLVM", 100))) %>% 
  ggplot() + 
  geom_boxplot(aes(y=log(MAE), x=method)) + theme(aspect.ratio=1)


data.frame(RMSE=c(RGC$RoGauRMSE[,1:2]), method=c(rep("Barcode", 100),
                                                 rep("Poisson GLLVM", 100))) %>% 
  ggplot() + 
  geom_boxplot(aes(y=log(RMSE), x=method)) + theme(aspect.ratio=1)

data.frame(MAE=c(RGC$RoGauMAE[,1:2]), method=c(rep("Barcode", 100),
                                        rep("Poisson GLLVM", 100))) %>% 
  ggplot() + 
  geom_boxplot(aes(y=log(MAE), x=method)) + 
  theme_bw() + theme(aspect.ratio=1)


data.frame(reshape::melt(PCU)) |>
  ggplot() + 
  geom_boxplot(aes(x=L1, y=log(value), fill=as.factor(X2)), width=0.75) + 
  theme_bw() + theme(aspect.ratio=1)

data.frame(reshape::melt(RGC)) |>
  filter(X2<3) |> filter(L1!="RoGauProc") |> filter(value<exp(100)) |>
  ggplot() + 
  geom_boxplot(aes(x=L1, y=log(value), fill=as.factor(X2)), width=0.75, outliers = F, position=position_dodge(1)) + 
  geom_point(aes(x=L1, y=log(value), fill=as.factor(X2)), 
             position=position_jitterdodge(dodge.width = 1, jitter.width = 0.5), shape="x", alpha=0.25, size=3) + 
  theme(aspect.ratio=1) + theme_bw()
  
data.frame(reshape::melt(RGC)) |>
  filter(X2<3) |> filter(L1!="RoGauProc") |> filter(value<exp(100)) |>
  ggplot() + 
  geom_point(aes(x=L1, y=log(value), fill=as.factor(X2)), 
             position=position_jitterdodge(dodge.width = 1, jitter.width = 0.25), 
             shape="*", alpha=0.75, size=4, color="red", show.legend = F) + 
  geom_violin_pattern(aes(x=L1, y=log(value), pattern=as.factor(X2)), width=1, 
                       position=position_dodge(1), pattern_density=0.25, pattern_spacing=0.01, fill="white") + 
  scale_pattern_manual(values=c('wave', 'stripe'), labels=c("Barcode", "Pois. GLLVM")) + 
  ylab("log error") + xlab("") + ggtitle("Concurrent: rounded Gaussian") + 
  scale_x_discrete(labels= c("MAE", "RMSE")) + 
  labs(pattern=NULL) + theme_bw() 


data.frame(reshape::melt(PCU)) |>
  ggplot() + 
  geom_point(aes(x=L1, y=log(value), fill=as.factor(X2)), 
             position=position_jitterdodge(dodge.width = 1, jitter.width = 0.25), 
             shape="*", alpha=0.75, size=4, color="red", show.legend = F) + 
  geom_violin_pattern(aes(x=L1, y=log(value), pattern=as.factor(X2)), width=1, 
                      position=position_dodge(1), pattern_density=0.25, pattern_spacing=0.01, fill="white") + 
  scale_pattern_manual(values=c('wave', 'stripe'), labels=c("Barcode", "Pois. GLLVM")) + 
  ylab("log error") + xlab("") + ggtitle("Unconstrained: Poisson copula") + 
  scale_x_discrete(labels= c("MAE", "RMSE")) + 
  labs(pattern=NULL) + theme_bw()

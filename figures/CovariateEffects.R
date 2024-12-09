library(ggpattern)
library(patchwork)

source("SetupFunctions.R")
outStack <- stackedList(readRDS("./results/birds6factors05Dec2024.rds"))

EB <- apply(outStack$Beta, c(1,2), mean)
EL <- apply(outStack$Beta, c(1,2), quantile, probs=0.025)
EU <- apply(outStack$Beta, c(1,2), quantile, probs=0.975)

ES <- sign(EL) == sign(EU)

Cfrac <- colMeans(apply(outStack$C, c(1,2), median))
Sfrac <- colMeans(apply(outStack$S, c(1,2), median))

xlabs <- paste(colnames(tmp$dat$X)[2:tmp$dims$Q])
xlabs <- stringr::str_replace_all(xlabs, "_", " ")
xlabs <- stringr::str_replace_all(xlabs, "poly", "")
xlabs <- stringr::str_replace_all(xlabs, "1", "")

myl <- 4
for (myl in 2:6){
a <- data.frame(q=xlabs, e=EB[2:tmp$dims$Q,myl],
           l=EL[2:tmp$dims$Q,myl], 
           u=EU[2:tmp$dims$Q,myl],
           s=ES[2:tmp$dims$Q,myl]) %>% 
  ggplot() + 
  geom_segment(aes(y=q, x=l, xend=u, color=s)) + 
  geom_vline(aes(xintercept=0)) + 
  geom_point(aes(y=q, x=e, color=s), shape=16, size=2) + 
  geom_point(aes(y=q, x=l, color=s), shape="(", size=3) + 
  geom_point(aes(y=q, x=u, color=s), shape=")", size=3) + 
  scale_color_manual(values=c("gray", "red")) + 
  theme_minimal() + 
  theme(legend.position="none",
        axis.text.y = element_text(angle=25, vjust=0.2)) + 
  xlab("") + ylab("") + ggtitle(paste("(",myl-1,")", sep=""))

b <- ggplot() + 
  geom_rect(aes(xmin=0,xmax=4,ymin=0,ymax=1), color="black", fill="white") + 
  geom_rect(aes(xmin=0,xmax=4,ymin=0,ymax=0.5), color="black", fill="white") + 
  geom_rect_pattern(aes(xmin=0,xmax=4*Cfrac[myl],ymin=0.5,ymax=1), color="black", fill="white", pattern="stripe") +
  geom_rect_pattern(aes(xmin=0,xmax=4*Sfrac[myl],ymin=0,ymax=0.5), color="black", fill="white",pattern="stripe") +
  theme_void() + coord_fixed() +
  xlim(c(-1,4.05)) + 
  annotate("text", x=-0.65, y=0.25, label="Species") +
  annotate("text", x=-0.65, y=0.75, label="Samples")

print((a / b))
}


op <- par(); plot.new()
rect(0,0,1,0.5)
rect(0,0,1,0.25)
rect(0,0.25,Cfrac[2],0.5,density = 30)
rect(0,0,Sfrac[2],0.25,density = 30)



plot(c(0, 10), c(0, 10), type = "n", xlab = "", ylab = "",
     main = "2 x 11 rectangles; 'rect(100+i,300+i,  150+i,380+i)'")
i <- 4*(0:10)
j <- 10*(0:5)
rect(125+j, 360+j,   141+j, 405+j/2, col = c(NA,0),
     border = "gold", lwd = 2)

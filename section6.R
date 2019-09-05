# Set working directory to source destination
setwd("~/Documents/Immune_variation/3.Manuscripts/MS2/V5/analysis_V5")
library(readxl)
library(reshape2)
library(ggplot2)
library(cowplot)
library(stats)
library(dplyr)
library(lme4)
library(gtools)
library(Hmisc)
library(factoextra)
library(ggfortify)
library(scales)
library(dplyr)
rm(list = ls())

########## SECTION 6 - correlation heatmaps and networks ##############

rm(list = ls())

SBF <- read.csv("input/SB_frequency.csv")
gates <- read.csv("input/Subsets_gates.csv")

# remove neutrophils and lymphomonocytes
SBF <- SBF[, !names(SBF) %in% c("TLC","Neutrophil", "Lymphomonocyte")] 
gates <- dplyr::filter(gates, !Name %in% c("TLC", "Neutrophil", "Lymphomonocyte"))

SBF_Means <- aggregate(cbind(`B_cell`, `B1`, `B_immature`,`B_memory`,
                             `B_naive`, Plasmablast, Plasmablast_IgA, T_cell, CD4, CD8, Gamma_delta_T,
                             NKT, iNKT, CD4_naive, CD4_memory, CD4_EMRA, CD8_naive, CD8_memory, CD8_EMRA, 
                             Treg, nTreg, iTreg, CD39Treg, CD39nTreg, CD39iTreg, 
                             Monocyte, Monocyte_Classical, Monocyte_Inflammatory, Monocyte_Patrolling,
                             DC, mDC, pDC, NK_cell) ~ ID, data = SBF, FUN = mean, na.rm = TRUE)
SBF_Means <- melt(SBF_Means)
names <- gates
names <- names[,-3]
SBF_Means <- merge(SBF_Means, names, by.x = "variable", by.y = "Name", all.x = T)
SBF_Means <- dcast(SBF_Means, ID~Subset)
SBF_VAR <- aggregate(cbind(`B_cell`, `B1`, `B_immature`,`B_memory`,
                           `B_naive`, Plasmablast, Plasmablast_IgA, T_cell, CD4, CD8, Gamma_delta_T,
                           NKT, iNKT, CD4_naive, CD4_memory, CD4_EMRA, CD8_naive, CD8_memory, CD8_EMRA, 
                           Treg, nTreg, iTreg, CD39Treg, CD39nTreg, CD39iTreg, 
                           Monocyte, Monocyte_Classical, Monocyte_Inflammatory, Monocyte_Patrolling,
                           DC, mDC, pDC, NK_cell) ~ ID, data = SBF, FUN = sd, na.rm = TRUE)
SBF_VAR <- data.frame(ID=SBF_VAR[,1],(SBF_VAR[,-1])^2)
SBF_VAR <- melt(SBF_VAR)
names <- gates
names <- names[,-3]
SBF_VAR <- merge(SBF_VAR, names, by.x = "variable", by.y = "Name", all.x = T)
SBF_VAR <- dcast(SBF_VAR, ID~Subset)
identical(colnames(SBF_Means), colnames(SBF_VAR))

# Correlation heatmap of SBF_Means 
# Reference : https://bit.ly/1vpTlg2
cormat <- round(cor(SBF_Means[,-1], use="pairwise.complete.obs", method = "spearman"),2)
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
pdf("./output/Fig6.pdf", width = 7, height = 7)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  coord_fixed()+ 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 1) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
dev.off()

## Get the p values 
library("Hmisc")
cor2 <- rcorr(as.matrix(SBF_Means[,-1]), type = "spearman")
# Extract the correlation coefficients
COR_MEANS_R <- get_upper_tri(cor2$r)
COR_MEANS_R <- melt(COR_MEANS_R)
colnames(COR_MEANS_R) <- c("var1", "var2", "Rho")
# Extract p-values
COR_MEANS_P <- get_upper_tri(cor2$P)
COR_MEANS_P <- melt(COR_MEANS_P)
colnames(COR_MEANS_P) <- c("var1", "var2", "P")
identical(COR_MEANS_R[,1:2], COR_MEANS_P[,1:2])
COR_MEANS <- cbind(COR_MEANS_R, COR_MEANS_P$P)
colnames(COR_MEANS) <- c("var1", "var2","Rho", "P")
COR_MEANS <- COR_MEANS[which(!is.na(COR_MEANS$Rho)),]
COR_MEANS <- COR_MEANS[which(COR_MEANS$Rho !=1),]
COR_MEANS <- COR_MEANS[which(COR_MEANS$Rho !=-1.0),]
plot(COR_MEANS$Rho, COR_MEANS$P)
COR_MEANS$AdjP <- p.adjust(COR_MEANS$P, "fdr", nrow(COR_MEANS))
plot(COR_MEANS$Rho, COR_MEANS$AdjP)
colnames(COR_MEANS) <- c("Variable 1", "Variable 2", "Correlation coefficient", "P-value", "Adjusted P (FDR)")
write.csv(COR_MEANS, "./output/TableS8.csv", row.names = F)

## Correlation between subsets not falling into the same gate 
## Positive correlations 
##Patrolling monocyte vs DC 
p1 <- ggplot(SBF_Means, aes(`Monocyte`, `Dendritic cell`))+
  geom_point()+
  theme_gray()
p1
## CD4 mem vs CD8 mem 
p2 <- ggplot(SBF_Means, aes(`CD4 memory`, `CD8 memory`))+
  geom_point()+
  theme_gray()
p2
## GDT vs NKT 
p3 <- ggplot(SBF_Means, aes(`Gamma delta T`, `NKT`))+
  geom_point()+
  theme_gray()
p3
## CD4EMRA vs PB 
p4 <- ggplot(SBF_Means, aes(`CD4 EMRA`, `Plasmablast`))+
  geom_point()+
  theme_gray()
p4

library(cowplot)
pdf("./output/FigS7.pdf", height = 7, width = 7)
plot_grid(p1, p2, p3, p4, nrow = 2)
dev.off()

##Negative correlations 

##CD4 vs GDT ##
p1 <- ggplot(SBF_Means, aes(`Gamma delta T`, `CD4`))+
  geom_point()+
  theme_gray()
p1
## Mono vs iNKT 
p2 <- ggplot(SBF_Means, aes(Monocyte, `iNKT`))+
  geom_point()+
  theme_gray()+
  scale_y_log10()
p2
## DC vs iNKT 
p3 <- ggplot(SBF_Means, aes(`Dendritic cell`, `iNKT`))+
  geom_point()+
  theme_gray()+
  scale_y_log10()
p3
## naive B vs DC+ 
p4 <- ggplot(SBF_Means, aes(`Dendritic cell`, `B naive`))+
  geom_point()+
  theme_gray()
p4
## B1 vs iNKT 
p5<- ggplot(SBF_Means, aes(`B1`, `iNKT`))+
  geom_point()+
  theme_gray()+
  scale_y_log10()
p5
## Monocyte vs NKT 
p6 <- ggplot(SBF_Means, aes(`Monocyte`, `NKT`))+
  geom_point()+
  theme_gray()
p6


library(cowplot)
pdf("./output/FigS8.pdf", height = 9, width = 6)
plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3)
dev.off()

# qgraph
library(qgraph)
library(RColorBrewer)


# Group names given manually

Names <- read.csv("temp/temp_names_section6.csv")
Names <- Names[order(Names$Group),]
rownames(Names) <- 1: nrow(Names)
order_i_want <- factor(Names$Name)

dat <- SBF_Means
cordat <- cor(dat[,-c(1)])
cordat <- melt(cordat)
cordat <- dcast(cordat, Var1~Var2)
cordat$Var1 <- factor(cordat$Var1, levels = order_i_want)

pdf("./output/Fig7.pdf", width = 7, height = 8)
qgraph(cordat[,-1], 
       vsize=5,
       labels=TRUE, 
       groups = Names$Group,
       theme ="classic",
       legend.cex=0.5,
       legend=TRUE,
       layout="groups", 
       #layout="spring",
       #layout="circle",
       legend.mode="nodeNames",
       nodeNames=T,
       minimum=0.2,
       cut=0.4,
       posCol= "darkgreen", 
       negCol="red")

dev.off()

## Correlation heatmap of SBF_VAR 
# Reference : https://bit.ly/1vpTlg2
cormat <- round(cor(SBF_VAR[,-1], use="pairwise.complete.obs", method = "spearman"),2)
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
pdf("./output/FigS9.pdf", width = 7, height = 7)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  coord_fixed()+ 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 1) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
dev.off()
## Get the p values 
library("Hmisc")
cor2 <- rcorr(as.matrix(SBF_VAR[,-1]), type = "spearman")
# Extract the correlation coefficients
COR_VAR_R <- get_upper_tri(cor2$r)
COR_VAR_R <- melt(COR_VAR_R)
colnames(COR_VAR_R) <- c("var1", "var2", "Rho")
# Extract p-values
COR_VAR_P <- get_upper_tri(cor2$P)
COR_VAR_P <- melt(COR_VAR_P)
colnames(COR_VAR_P) <- c("var1", "var2", "P")
identical(COR_VAR_R[,1:2], COR_VAR_P[,1:2])
COR_VAR <- cbind(COR_VAR_R, COR_VAR_P$P)
colnames(COR_VAR) <- c("var1", "var2","Rho", "P")
COR_VAR <- COR_VAR[which(!is.na(COR_VAR$Rho)),]
COR_VAR <- COR_VAR[which(COR_VAR$Rho !=1),]
COR_VAR <- COR_VAR[which(COR_VAR$Rho !=-1.0),]
plot(COR_VAR$Rho, COR_VAR$P)
COR_VAR$AdjP <- p.adjust(COR_VAR$P, "fdr", nrow(COR_VAR))
plot(COR_VAR$Rho, COR_VAR$AdjP)
colnames(COR_VAR) <- c("Variable 1", "Variable 2", "Correlation coefficient", "P-value", "Adjusted P (FDR)")
write.csv(COR_VAR, "output/TableS9.csv", row.names = F)



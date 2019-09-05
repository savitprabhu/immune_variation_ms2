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


SBF <- read.csv("input/SB_frequency.csv")
gates <- read.csv("input/Subsets_gates.csv")

# remove neutrophils and lymphomonocytes
SBF <- SBF[, !names(SBF) %in% c("TLC","Neutrophil", "Lymphomonocyte")] 
gates <- dplyr::filter(gates, !Name %in% c("TLC", "Neutrophil", "Lymphomonocyte"))

SBF[,c(4:36)] = scale(SBF[,c(4:36)])
str(SBF)
SBF$Parameter = NULL
SBF_Means <- aggregate(SBF[,c(3:35)],by = list(SBF$ID), data = SBF, FUN = mean, na.rm = TRUE)
colnames(SBF_Means)[1] = "ID"

SBF_Means <- melt(SBF_Means)
names <- gates
names <- names[,-3]
SBF_Means <- merge(SBF_Means, names, by.x = "variable", by.y = "Name", all.x = T)
SBF_Means <- dcast(SBF_Means, ID~Subset)

SBF_VAR <- aggregate(SBF[,c(3:35)], by = list(SBF$ID), FUN = sd, na.rm = TRUE)

SBF_VAR <- data.frame(ID=SBF_VAR[,1],(SBF_VAR[,-1])^2)
SBF_VAR <- melt(SBF_VAR)
names <- gates
names <- names[,-3]
SBF_VAR <- merge(SBF_VAR, names, by.x = "variable", by.y = "Name", all.x = T)
SBF_VAR <- dcast(SBF_VAR, ID~Subset)
identical(colnames(SBF_Means), colnames(SBF_VAR))

########## Fig 3: Correltions between Mean and Var ###########
SBF_Means <- melt(SBF_Means)
SBF_Means$Parameter <- rep("Mean")
SBF_Means$Parameter <- paste(SBF_Means$Parameter, SBF_Means$variable, sep = "_")
SBF_Means <- SBF_Means[,-2]
SBF_VAR <- melt(SBF_VAR)
SBF_VAR$Parameter <- rep("VAR")
SBF_VAR$Parameter <- paste(SBF_VAR$Parameter, SBF_VAR$variable, sep = "_")
SBF_VAR <- SBF_VAR[,-2]
SBF_MV <- rbind(SBF_Means, SBF_VAR)
SBF_MV <- dcast(SBF_MV, ID~Parameter)

# Figure 3 
cor.test(SBF_MV$`Mean_CD4 memory`, SBF_MV$`VAR_CD4 memory`, method = 'spearman')
p1<-ggplot(SBF_MV, aes(`Mean_CD4 memory`, `VAR_CD4 memory`))+
  theme_bw()+
  geom_point(size=2, colour="black", shape = 1)+
  labs(title ="CD4 memory \n(R=-0.1, p=0.5)",
       x = "Mean of z-scores \n(baseline levels)",
       y = "Variance of z-scores (log scale) \n(corrected for subset mean)")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title = element_text(size = 10, face = "bold"))+
  theme(axis.text = element_text(size = 9, face = "bold"))+
  theme(axis.line = element_line(colour = "black"))+
  #scale_x_log10()+
  scale_y_log10()
p1

cor.test(SBF_MV$`Mean_Monocyte Classical`, SBF_MV$`VAR_Monocyte Classical`, method = 'spearman')
p2<-ggplot(SBF_MV, aes(`Mean_Monocyte Classical`, `VAR_Monocyte Classical`))+
  theme_bw()+
  geom_point(size=2, colour="black", shape = 1)+
  labs(title ="Classical monocyte \n(R=-0.02, p=0.8)",
       x = "Mean of z-scores \n(baseline levels)",
       y = "Variance of z-scores (log scale) \n(corrected for subset mean)")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title = element_text(size = 10, face = "bold"))+
  theme(axis.text = element_text(size = 9, face = "bold"))+
  theme(axis.line = element_line(colour = "black"))+
  #scale_x_log10()+
  scale_y_log10()
p2

cor.test(SBF_MV$`Mean_Treg`, SBF_MV$`VAR_Treg`, method = 'spearman')
p3<-ggplot(SBF_MV, aes(`Mean_Treg`, `VAR_Treg`))+
  theme_bw()+
  geom_point(size=2, colour="black", shape = 1)+
  labs(title ="Treg \n(R=0.27, p=0.07)",
       x = "Mean of z-scores \n(baseline levels)",
       y = "Variance of z-scores (log scale) \n(corrected for subset mean)")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title = element_text(size = 10, face = "bold"))+
  theme(axis.text = element_text(size = 9, face = "bold"))+
  theme(axis.line = element_line(colour = "black"))+
  #scale_x_log10()+
  scale_y_log10()
p3

cor.test(SBF_MV$`Mean_Plasmablast`, SBF_MV$`VAR_Plasmablast`, method = 'spearman')
p4<-ggplot(SBF_MV, aes(`Mean_Plasmablast`, `VAR_Plasmablast`))+
  theme_bw()+
  geom_point(size=2, colour="black", shape = 1)+
  labs(title ="Plasmablast \n(R=0.92, p< 2.2e-16)",
       x = "Mean of z-scores \n(baseline levels)",
       y = "Variance of z-scores (log scale) \n(corrected for subset mean)")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title = element_text(size = 10, face = "bold"))+
  theme(axis.text = element_text(size = 9, face = "bold"))+
  theme(axis.line = element_line(colour = "black"))+
  #scale_x_log10()+
  scale_y_log10()
p4

cor.test(SBF_MV$`Mean_CD4 EMRA`, SBF_MV$`VAR_CD4 EMRA`, method = 'spearman')
p5<-ggplot(SBF_MV, aes(`Mean_CD4 EMRA`, `VAR_CD4 EMRA`))+
  theme_bw()+
  geom_point(size=2, colour="black", shape = 1)+
  labs(title ="CD4 EMRA \n(R=0.77, p = 1.1e-10)",
       x = "Mean of z-scores \n(baseline levels)",
       y = "Variance of z-scores (log scale) \n(corrected for subset mean)")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title = element_text(size = 10, face = "bold"))+
  theme(axis.text = element_text(size = 9, face = "bold"))+
  theme(axis.line = element_line(colour = "black"))+
  #scale_x_log10()+
  scale_y_log10()
p5

cor.test(SBF_MV$`Mean_iTreg`, SBF_MV$`VAR_iTreg`, method = 'spearman')
p6<-ggplot(SBF_MV, aes(`Mean_iTreg`, `VAR_iTreg`))+
  theme_bw()+
  geom_point(size=2, colour="black", shape = 1)+
  labs(title ="iTreg \n(R=0.7, p=3.2e-07)",
       x = "Mean of z-scores \n(baseline levels)",
       y = "Variance of z-scores (log scale) \n(corrected for subset mean)")+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title = element_text(size = 10, face = "bold"))+
  theme(axis.text = element_text(size = 9, face = "bold"))+
  theme(axis.line = element_line(colour = "black"))+
  #scale_x_log10()+
  scale_y_log10()
p6

library(cowplot)
pdf("output/Fig3.pdf", height = 7, width = 10)
plot_grid(p1,p2, p3, p4, p5, p6, nrow = 2)
dev.off()

library(Hmisc)
cor2 <- rcorr(as.matrix(SBF_MV[,-1]), type = "spearman")
# Extract the correlation coefficients
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
COR_R <- get_upper_tri(cor2$r)
COR_R <- melt(COR_R)
colnames(COR_R) <- c("var1", "var2", "Rho")
# Extract p-values
COR_P <- get_upper_tri(cor2$P)
COR_P <- melt(COR_P)
colnames(COR_P) <- c("var1", "var2", "P")
identical(COR_R[,1:2], COR_P[,1:2])
COR <- cbind(COR_R, COR_P$P)
colnames(COR) <- c("var1", "var2","Rho", "P")
COR <- COR[which(!is.na(COR$Rho)),]
COR <- COR[which(COR$Rho !=1),]
library(stringr)
COLUMN1 <- str_split_fixed(COR$var1, "_", 2)
COLUMN2 <- str_split_fixed(COR$var2, "_", 2)
COR_NEW <- cbind(COLUMN1, COLUMN2, COR)
colnames(COR_NEW) <- c("A", "B", "C", "D", "var1", "var2", "Rho", "P")
COR_NEW <- COR_NEW[which(COR_NEW$B==COR_NEW$D),]
COR_NEW$AdjP <- p.adjust(COR_NEW$P, "fdr", nrow(COR_NEW))
COR_NEW$Significance <- COR_NEW$AdjP
COR_NEW[which(COR_NEW$Significance<0.05),10] <- rep("Adj. p < 0.05")
COR_NEW[which(COR_NEW$Significance != "Adj. p < 0.05"),10] <- rep("Adj. p > 0.05")
COR_NEW <- COR_NEW[order(COR_NEW$Rho),]
COR_NEW$B <- factor(COR_NEW$B, levels = COR_NEW$B)
COR_NEW <- COR_NEW[,c(2, 7:10)]
colnames(COR_NEW) <- c("Subset", "Rho", "P", "FDR", "Significance")
write.csv(COR_NEW, "output/TableS7.csv", row.names = F)


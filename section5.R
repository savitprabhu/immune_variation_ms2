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

########## SECTION 5 PCA analysis ##############
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
SBF_VAR[,-1] = log(SBF_VAR[,-1]) # to avoid extreme outliers in PCA plots
identical(colnames(SBF_Means), colnames(SBF_VAR))

# PCA analysis 
# PCA Based on means alone 
forpca1 <- SBF_Means[,-1]
str(forpca1)
summary(forpca1)
forplot1 <- prcomp(forpca1)
p1<- fviz_eig(forplot1, title="", ylab = "% of explained variance")
p2<- fviz_pca_ind(forplot1, geom = c("point"), title= "")
p1
p2

# PCA Based on variance alone 
forpca2 <- SBF_VAR[,-1]
forplot2 <- prcomp(forpca2)
p3<- fviz_eig(forplot2, title="", ylab = "% of explained variance")
p4<- fviz_pca_ind(forplot2, geom = c("point"), title= "")
p3
p4

# PCA Based on means + variance 
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

forpca3 <- SBF_MV[,-1]
forplot3 <- prcomp(forpca3)
p5<- fviz_eig(forplot3, title="", ylab = "% of explained variance")
p6<- fviz_pca_ind(forplot3, geom = c("point"), title= "")
p5
p6 # 2 outliers are causing the problem

pdf("./output/Fig5.pdf", height = 7.5, width = 5.5)
plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3, labels = "AUTO")
dev.off()

############# Reviewer's question about variable individuals ############
rm(list = ls())

demog = read.csv("input/SB_demographic.csv")
demog = demog[,c(1:3)]

SBF <- read.csv("input/SB_frequency.csv")
gates <- read.csv("input/Subsets_gates.csv")

SBF <- SBF[, !names(SBF) %in% c("TLC","Neutrophil", "Lymphomonocyte")] 
gates <- dplyr::filter(gates, !Name %in% c("TLC", "Neutrophil", "Lymphomonocyte"))

SBF[,c(4:36)] = scale(SBF[,c(4:36)])
str(SBF)
SBF$Parameter = NULL

SBF_VAR <- aggregate(SBF[,c(3:35)], by = list(SBF$ID), FUN = sd, na.rm = TRUE)
SBF_VAR <- data.frame(ID=SBF_VAR[,1],(SBF_VAR[,-1])^2)
SBF_VAR <- melt(SBF_VAR)
names <- gates
names <- names[,-3]
SBF_VAR <- merge(SBF_VAR, names, by.x = "variable", by.y = "Name", all.x = T)
SBF_VAR <- dcast(SBF_VAR, ID~Subset)
SBF_VAR[,-1] = log(SBF_VAR[,-1]) # to normalize and avoid extreme outliers 

SBF_VAR = merge(demog, SBF_VAR)

pcaplot <- prcomp(SBF_VAR[,-c(1:3)])

p1 = fviz_pca_ind(pcaplot, geom = c("point"), title= "Gender distribution", habillage=SBF_VAR$Sex)
p2 = fviz_pca_ind(pcaplot, geom = c("point"), title= "Age distribution", habillage=SBF_VAR$Age)

pdf("output/FigS6.pdf", width = 5, height = 9)
plot_grid(p1, p2, labels = "AUTO", nrow = 2)
dev.off()

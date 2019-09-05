# Set working directory to source destination
setwd("~/Documents/Immune_variation/3.Manuscripts/MS2/V5/analysis_V5")
library(readxl)
library(reshape2)
library(ggplot2)
library(cowplot)
library(stats)
library(dplyr)
library(scales)
rm(list = ls())

########## Boxplots of all subsets ##############

rm(list = ls())
SBF <- read.csv("input/SB_frequency.csv")
SBFm <- melt(SBF, id.vars = c("ID", "Bleed", "Parameter"))
pdf("output/FigS4.pdf", height = 5, width = 8)
# Neutrophils 
CELL <- dplyr::filter(SBFm, variable == "Neutrophil")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")
p2<- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("Neutrophils",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## LM ##
CELL <- dplyr::filter(SBFm, variable == "Lymphomonocyte")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))

title <- ggdraw() + draw_label("Lymphomonocyte",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## B cell ##
CELL <- dplyr::filter(SBFm, variable == "B_cell")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")
p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("B cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## B1 B cell ##
CELL <- dplyr::filter(SBFm, variable == "B1")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("B1 B cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## Immature B cell ##
CELL <- dplyr::filter(SBFm, variable == "B_immature")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("Immature B cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## Memory B cell ##
CELL <- dplyr::filter(SBFm, variable == "B_memory")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("Memory B cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## Naive B cell ##
CELL <- dplyr::filter(SBFm, variable == "B_naive")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("Naive B cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## PB ##
CELL <- dplyr::filter(SBFm, variable == "Plasmablast")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("Plasmablast",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## PB IgA+ ##
CELL <- dplyr::filter(SBFm, variable == "Plasmablast_IgA")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("IgA+ Plasmablast",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## T ##
CELL <- dplyr::filter(SBFm, variable == "T_cell")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("T cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## CD4 T ##
CELL <- dplyr::filter(SBFm, variable == "CD4")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("CD4 T cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## CD8 T ##
CELL <- dplyr::filter(SBFm, variable == "CD8")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("CD8 T cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## GD T ##
CELL <- dplyr::filter(SBFm, variable == "Gamma_delta_T")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("Gamma delta T cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## NKT ##
CELL <- dplyr::filter(SBFm, variable == "NKT")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("NKT cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## iNKT ##
CELL <- dplyr::filter(SBFm, variable == "iNKT")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("iNKT cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## CD4 naive ##
CELL <- dplyr::filter(SBFm, variable == "CD4_naive")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("naive CD4 T cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## CD4 memory ##
CELL <- dplyr::filter(SBFm, variable == "CD4_memory")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("memory CD4 T cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## CD4 EMRA ##
CELL <- dplyr::filter(SBFm, variable == "CD4_EMRA")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label(" CD4 T EMRA cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## CD8 naive ##
CELL <- dplyr::filter(SBFm, variable == "CD8_naive")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("naive CD8 T cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## CD8 memory ##
CELL <- dplyr::filter(SBFm, variable == "CD8_memory")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("memory CD8 T cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## CD8 EMRA ##
CELL <- dplyr::filter(SBFm, variable == "CD8_EMRA")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label(" CD8 T EMRA cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## Treg ##
CELL <- dplyr::filter(SBFm, variable == "Treg")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label(" CD4 Tregs ",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## nTreg ##
CELL <- dplyr::filter(SBFm, variable == "nTreg")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("nTregs ",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## iTreg ##
CELL <- dplyr::filter(SBFm, variable == "iTreg")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("iTregs ",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## CD39+ Treg ##
CELL <- dplyr::filter(SBFm, variable == "CD39Treg")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label(" CD39+ Tregs ",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## CD39+ nTreg ##
CELL <- dplyr::filter(SBFm, variable == "CD39nTreg")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label(" CD39+ nTregs ",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## CD39+ iTreg ##
CELL <- dplyr::filter(SBFm, variable == "CD39iTreg")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label(" CD39+ iTregs ",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## Mono ##
CELL <- dplyr::filter(SBFm, variable == "Monocyte")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("Monocyte",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## Mono_classical ##
CELL <- dplyr::filter(SBFm, variable == "Monocyte_Classical")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("Monocyte - Classical",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## Mono_infl ##
CELL <- dplyr::filter(SBFm, variable == "Monocyte_Inflammatory")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("Monocyte - Inflammatory",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## Mono_patrolling ##
CELL <- dplyr::filter(SBFm, variable == "Monocyte_Patrolling")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("Monocyte - Patrolling",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## DC ##
CELL <- dplyr::filter(SBFm, variable == "DC")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("Dendritic cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## mDC ##
CELL <- dplyr::filter(SBFm, variable == "mDC")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("mDC",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## pDC ##
CELL <- dplyr::filter(SBFm, variable == "pDC")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("pDC",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

## NK cell ##
CELL <- dplyr::filter(SBFm, variable == "NK_cell")
p1 <-ggplot(CELL, aes(variable, value))+
  geom_boxplot(fill="green", alpha=0.7)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", colour = "black"))+
  labs(x="Total", y="Frequency (%)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
CELLw <- dcast(CELL, ID+variable~Bleed, value.var = "value")
CELLw$Median <- apply(CELLw[,-c(1:2)], 1, median, na.rm=T)
CELLw <- CELLw[order(CELLw$Median),]
CELLw$ID <- factor(CELLw$ID, levels=unique(CELLw$ID))
CELLw$Median <- NULL
CELLw <- melt(CELLw, id.vars = c("ID", "variable"))
colnames(CELLw) <- c("ID", "subset", "variable", "value")

p2 <- ggplot(CELLw, aes(ID, value))+
  geom_boxplot(fill="skyblue", alpha=0.6, outlier.size = 0.3)+
  theme_light()+
  theme(axis.text.x = element_blank())+
  labs(x="Individuals", y=NULL)+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
title <- ggdraw() + draw_label("NK cell",fontface = 'bold')
p <- plot_grid(p1, p2, align = "h", nrow = 1, rel_widths = c(1/4, 3/4))
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins

dev.off()

############ Variance components model #############
library(lme4)
rm(list = ls())

SBF <- read.csv("input/SB_frequency.csv")
gates <- read.csv("input/Subsets_gates.csv")

# remove neutrophils and lymphomonocytes
SBF <- SBF[, !names(SBF) %in% c("TLC","Neutrophil", "Lymphomonocyte")] 
gates <- dplyr::filter(gates, !Name %in% c("TLC", "Neutrophil", "Lymphomonocyte"))


# Frequencies 
rm(list=setdiff(ls(), c("SBF", "gates")))

dat<- SBF
dat$ID<- as.factor(dat$ID)
dat$Bleed <- as.factor((dat$Bleed))
dat[,-c(1:3)]<-log(dat[,-c(1:3)])
dat[dat == -Inf] <- NA
datm<-melt(dat)

L<-lmer(B_cell~ (1|ID),data=dat)
B <- as.data.frame(VarCorr(L))
L<-lmer(B1~ (1|ID),data=dat)
B1 <- as.data.frame(VarCorr(L))
L<-lmer(B_immature~ (1|ID),data=dat)
B_immature <- as.data.frame(VarCorr(L))
L<-lmer(B_memory~ (1|ID),data=dat)
B_memory <- as.data.frame(VarCorr(L))
L<-lmer(B_naive~ (1|ID),data=dat)
B_naive <- as.data.frame(VarCorr(L))
L<-lmer(Plasmablast~ (1|ID),data=dat)
PB <- as.data.frame(VarCorr(L))
L<-lmer(Plasmablast_IgA~ (1|ID),data=dat)
PB_IgA <- as.data.frame(VarCorr(L))
L<-lmer(T_cell~ (1|ID),data=dat)
T_cell <- as.data.frame(VarCorr(L))
L<-lmer(CD4~ (1|ID),data=dat)
CD4 <- as.data.frame(VarCorr(L))
L<-lmer(CD8~ (1|ID),data=dat)
CD8 <- as.data.frame(VarCorr(L))
L<-lmer(Gamma_delta_T~ (1|ID),data=dat)
GDT <- as.data.frame(VarCorr(L))
L<-lmer(NKT~ (1|ID),data=dat)
NKT <- as.data.frame(VarCorr(L))
L<-lmer(iNKT~ (1|ID),data=dat)
iNKT <- as.data.frame(VarCorr(L))
L<-lmer(CD4_naive~ (1|ID),data=dat)
CD4n <- as.data.frame(VarCorr(L))
L<-lmer(CD4_memory~ (1|ID),data=dat)
CD4m <- as.data.frame(VarCorr(L))
L<-lmer(CD4_EMRA~ (1|ID),data=dat)
CD4EMRA <- as.data.frame(VarCorr(L))
L<-lmer(CD8_naive~ (1|ID),data=dat)
CD8n <- as.data.frame(VarCorr(L))
L<-lmer(CD8_memory~ (1|ID),data=dat)
CD8m <- as.data.frame(VarCorr(L))
L<-lmer(CD8_EMRA~ (1|ID),data=dat)
CD8EMRA <- as.data.frame(VarCorr(L))
L<-lmer(Treg~ (1|ID),data=dat)
Treg <- as.data.frame(VarCorr(L))
L<-lmer(nTreg~ (1|ID),data=dat)
nTreg <- as.data.frame(VarCorr(L))
L<-lmer(iTreg~ (1|ID),data=dat)
iTreg <- as.data.frame(VarCorr(L))
L<-lmer(CD39Treg~ (1|ID),data=dat)
CD39Treg <- as.data.frame(VarCorr(L))
L<-lmer(CD39nTreg~ (1|ID),data=dat)
CD39nTreg <- as.data.frame(VarCorr(L))
L<-lmer(CD39iTreg~ (1|ID),data=dat)
CD39iTreg <- as.data.frame(VarCorr(L))
L<-lmer(Monocyte~ (1|ID),data=dat)
Mono <- as.data.frame(VarCorr(L))
L<-lmer(Monocyte_Classical~ (1|ID),data=dat)
CM <- as.data.frame(VarCorr(L))
L<-lmer(Monocyte_Inflammatory~ (1|ID),data=dat)
IM <- as.data.frame(VarCorr(L))
L<-lmer(Monocyte_Patrolling~ (1|ID),data=dat)
PM <- as.data.frame(VarCorr(L))
L<-lmer(DC~ (1|ID),data=dat)
DC <- as.data.frame(VarCorr(L))
L<-lmer(mDC~ (1|ID),data=dat)
mDC <- as.data.frame(VarCorr(L))
L<-lmer(pDC~ (1|ID),data=dat)
pDC <- as.data.frame(VarCorr(L))
L<-lmer(NK_cell~ (1|ID),data=dat)
NK <- as.data.frame(VarCorr(L))

freq <- data.frame(B$vcov, B1$vcov, B_immature$vcov, B_memory$vcov, B_naive$vcov,
                   PB$vcov, PB_IgA$vcov, T_cell$vcov, CD4$vcov, CD8$vcov, GDT$vcov, NKT$vcov, iNKT$vcov,
                   CD4n$vcov, CD4m$vcov, CD4EMRA$vcov, CD8n$vcov, CD8m$vcov, CD8EMRA$vcov,
                   Treg$vcov, nTreg$vcov, iTreg$vcov, CD39Treg$vcov, CD39iTreg$vcov, CD39nTreg$vcov,
                   Mono$vcov, CM$vcov, IM$vcov, PM$vcov, DC$vcov, mDC$vcov, pDC$vcov, NK$vcov)
rownames(freq)<- c("Sample.ID", "Residual")
freq <- data.frame(t(freq))
freq$Subset <- gsub("\\..*","",rownames(freq))
freq$ICC <- freq$Residual/(freq$Residual+freq$Sample.ID)
freq$Subset <- factor(freq$Subset, levels = freq$Subset)
toplot <- freq[,c(3,4)]
rownames(toplot) <- 1:nrow(toplot)
toplot$Between <- 1-toplot$ICC
colnames(toplot) <- c("Name", "Within-individual", "Between-individual")
identical(as.character(toplot$Name), colnames(dat[,-c(1:3)]))
toplot$Name
colnames(dat[,-c(1:3)])

toplot$Name <- colnames(dat[,-c(1:3)]) # verified manually
gates <- gates[,-3]
toplot <- merge(toplot, gates, all.x = T)
toplot <- toplot[,c(4, 2, 3)]
toplot <- toplot[order(toplot$`Within-individual`),]
toplot$Subset <- factor(toplot$Subset, levels = unique(toplot$Subset))
#write.csv(toplot, "Var_comp_freq.csv", row.names = F)
toplot <- melt(toplot)
toplot$variable <- factor(toplot$variable, levels = c("Between-individual","Within-individual"))

ggplot(toplot, aes(Subset, value, fill=variable))+
  geom_bar(stat = "identity", width = 0.5, colour="black")+
  labs(x=NULL, y="Variance", fill="Partition")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  theme(text = element_text(face = "bold", colour = "black"))+
  scale_y_continuous(expand = c(0,0))

CV = read.csv("output/TableS4.csv")
CV = merge(CV, gates)
CV = CV[, c("Subset", "CV_betw_corrected")]
toplot2 = merge(CV, toplot)
toplot2 = toplot2[order(toplot2$CV_betw_corrected),]
toplot2$Subset = factor(toplot2$Subset, levels = unique(toplot2$Subset))


toplot2 = dcast(toplot2, Subset+CV_betw_corrected ~ variable)
str(toplot2)

toplot2$`Within-individual` = toplot2$CV_betw_corrected*toplot2$`Within-individual`
toplot2$`Between-individual` = toplot2$CV_betw_corrected*toplot2$`Between-individual`
toplot2$CV_betw_corrected = NULL
toplot2 = melt(toplot2)
str(toplot2)

# reviewer wanted the plot in this fashion.
pdf("output/Fig2.pdf", width = 7, height = 4)
ggplot(toplot2, aes(Subset, value, fill = variable))+
  geom_bar(stat = "identity" , width = 0.5, colour="black")+
  labs(x=NULL, y="Between-individual CV \n(corrected for technical CV)", 
       fill="Partition")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust = 1))+
  theme(text = element_text(face = "bold", colour = "black"))+
  scale_y_continuous(expand = c(0,0))
dev.off()

########## Comparison of within- vs between- variances (hypothesis testing) ##########

# Frequency 
rm(list = ls())
DATA <- read.csv("input/SB_frequency.csv")
DATA$Parameter <- NULL
DAT <- melt(DATA)
str(DAT)

pdf("output/FigS5.pdf", height = 10, width = 8)
par(mfrow=c(4,3))

# Neutrophil
dat <- dplyr::filter(DAT, variable =="Neutrophil")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Neutrophil<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Neutrophil
set.seed(12345)
avgdiff_Neutrophil <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Neutrophil, xlim = c(-abs(obsdiff_Neutrophil)-20, abs(obsdiff_Neutrophil)+20 ),
     main = "Neutrophil \n (% of total WBCs)", xlab = "Null distribution")
abline(v=obsdiff_Neutrophil, col="red",lwd=2)
Pval_Neutrophil<-(sum(abs(avgdiff_Neutrophil) > abs(obsdiff_Neutrophil))+1) / (length(avgdiff_Neutrophil)+1)

# Lymphomonocyte
dat <- dplyr::filter(DAT, variable =="Lymphomonocyte")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Lymphomonocyte<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Lymphomonocyte
set.seed(12345)
avgdiff_Lymphomonocyte <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Lymphomonocyte, xlim = c(-abs(obsdiff_Lymphomonocyte)-20, abs(obsdiff_Lymphomonocyte)+20 ),
     main = "Lymphomonocyte \n(% of total WBCs)", xlab = "Null distribution")
abline(v=obsdiff_Lymphomonocyte, col="red",lwd=2)

Pval_Lymphomonocyte<-(sum(abs(avgdiff_Lymphomonocyte) > abs(obsdiff_Lymphomonocyte))+1) / (length(avgdiff_Lymphomonocyte)+1)

# B_cell
dat <- dplyr::filter(DAT, variable =="B_cell")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_B_cell<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_B_cell
set.seed(12345)
avgdiff_B_cell <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_B_cell, xlim = c(-abs(obsdiff_B_cell)-20, abs(obsdiff_B_cell)+20 ),
     main = "B cell \n(% of total) ", xlab = "Null distribution")
abline(v=obsdiff_B_cell, col="red",lwd=2)

Pval_B_cell<-(sum(abs(avgdiff_B_cell) > abs(obsdiff_B_cell))+1) / (length(avgdiff_B_cell)+1)

# B1
dat <- dplyr::filter(DAT, variable =="B1")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_B1<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_B1
set.seed(12345)
avgdiff_B1 <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_B1, xlim = c(-abs(obsdiff_B1)-5, abs(obsdiff_B1)+5 ),
     main = "B1 B cell \n(% of B cells)", xlab = "Null distribution")
abline(v=obsdiff_B1, col="red",lwd=2)
Pval_B1<-(sum(abs(avgdiff_B1) > abs(obsdiff_B1))+1) / (length(avgdiff_B1)+1)

# B_immature
dat <- dplyr::filter(DAT, variable =="B_immature")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_B_immature<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_B_immature
set.seed(12345)
avgdiff_B_immature <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_B_immature, xlim = c(-abs(obsdiff_B_immature)-20, abs(obsdiff_B_immature)+20 ),
     main = "Immature B cells \n(% of B cells)", xlab = "Null distribution")
abline(v=obsdiff_B_immature, col="red",lwd=2)
Pval_B_immature<-(sum(abs(avgdiff_B_immature) > abs(obsdiff_B_immature))+1) / (length(avgdiff_B_immature)+1)

# B memory #
dat <- dplyr::filter(DAT, variable =="B_memory")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Bmemory<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Bmemory
set.seed(12345)
avgdiff_Bmemory <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Bmemory, xlim = c(-abs(obsdiff_Bmemory)-50, abs(obsdiff_Bmemory)+50 ),
     main = "B memory \n(% of B)", xlab = "Null distribution")
abline(v=obsdiff_Bmemory, col="red",lwd=2)
Pval_B_memory<-(sum(abs(avgdiff_Bmemory) > abs(obsdiff_Bmemory))+1) / (length(avgdiff_Bmemory)+1)


# B naive #
dat <- dplyr::filter(DAT, variable =="B_naive")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Bnaive<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Bnaive
set.seed(12345)
avgdiff_Bnaive <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Bnaive, xlim = c(-abs(obsdiff_Bnaive)-50, abs(obsdiff_Bnaive)+50 ),
     main = "B naive \n(% of B cells)", xlab = "Null distribution")
abline(v=obsdiff_Bnaive, col="red",lwd=2)

Pval_B_naive<-(sum(abs(avgdiff_Bnaive) > abs(obsdiff_Bnaive))+1) / (length(avgdiff_Bnaive)+1)

# Plasmablast
dat <- dplyr::filter(DAT, variable =="Plasmablast")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Plasmablast<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Plasmablast
set.seed(12345)
avgdiff_Plasmablast <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Plasmablast, xlim = c(-abs(obsdiff_Plasmablast)-10, abs(obsdiff_Plasmablast)+10 ),
     main = "Plasmablast \n(% of B cells)", xlab = "Null distribution")
abline(v=obsdiff_Plasmablast, col="red",lwd=2)

Pval_Plasmablast<-(sum(abs(avgdiff_Plasmablast) > abs(obsdiff_Plasmablast))+1) / (length(avgdiff_Plasmablast)+1)

# Plasmablast_IgA
dat <- dplyr::filter(DAT, variable =="Plasmablast_IgA")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Plasmablast_IgA<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Plasmablast_IgA
set.seed(12345)
avgdiff_Plasmablast_IgA <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Plasmablast_IgA, xlim = c(-abs(obsdiff_Plasmablast_IgA)-5, abs(obsdiff_Plasmablast_IgA)+5 ),
     main = "Plasmablast_IgA \n(% of B cells)", xlab = "Null distribution")
abline(v=obsdiff_Plasmablast_IgA, col="red",lwd=2)

Pval_Plasmablast_IgA<-(sum(abs(avgdiff_Plasmablast_IgA) > abs(obsdiff_Plasmablast_IgA))+1) / (length(avgdiff_Plasmablast_IgA)+1)

# T_cell
dat <- dplyr::filter(DAT, variable =="T_cell")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_T_cell<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_T_cell
set.seed(12345)
avgdiff_T_cell <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_T_cell, xlim = c(-abs(obsdiff_T_cell)-50, abs(obsdiff_T_cell)+50),
     main = "T_cell \n(% of total)", xlab = "Null distribution")
abline(v=obsdiff_T_cell, col="red",lwd=2)

Pval_T_cell<-(sum(abs(avgdiff_T_cell) > abs(obsdiff_T_cell))+1) / (length(avgdiff_T_cell)+1)


# CD4
dat <- dplyr::filter(DAT, variable =="CD4")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_CD4<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_CD4
set.seed(12345)
avgdiff_CD4 <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_CD4, xlim = c(-abs(obsdiff_CD4)-5, abs(obsdiff_CD4)+5 ),
     main = "CD4 \n(% of T cells)", xlab = "Null distribution")
abline(v=obsdiff_CD4, col="red",lwd=2)

Pval_CD4<-(sum(abs(avgdiff_CD4) > abs(obsdiff_CD4))+1) / (length(avgdiff_CD4)+1)

# CD8
dat <- dplyr::filter(DAT, variable =="CD8")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_CD8<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_CD8
set.seed(12345)
avgdiff_CD8 <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_CD8, xlim = c(-abs(obsdiff_CD8)-5, abs(obsdiff_CD8)+5 ),
     main = "CD8 \n(% of T cells)", xlab = "Null distribution")
abline(v=obsdiff_CD8, col="red",lwd=2)
Pval_CD8<-(sum(abs(avgdiff_CD8) > abs(obsdiff_CD8))+1) / (length(avgdiff_CD8)+1)

# Gamma_delta_T
dat <- dplyr::filter(DAT, variable =="Gamma_delta_T")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Gamma_delta_T<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Gamma_delta_T
set.seed(12345)
avgdiff_Gamma_delta_T <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Gamma_delta_T, xlim = c(-abs(obsdiff_Gamma_delta_T)-5, abs(obsdiff_Gamma_delta_T)+5 ),
     main = "Gamma_delta_T \n(% of T cells)", xlab = "Null distribution")
abline(v=obsdiff_Gamma_delta_T, col="red",lwd=2)

Pval_Gamma_delta_T<-(sum(abs(avgdiff_Gamma_delta_T) > abs(obsdiff_Gamma_delta_T))+1) / (length(avgdiff_Gamma_delta_T)+1)

# NKT
dat <- dplyr::filter(DAT, variable =="NKT")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_NKT<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_NKT
set.seed(12345)
avgdiff_NKT <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_NKT, xlim = c(-abs(obsdiff_NKT)-1, abs(obsdiff_NKT)+1 ),
     main = "NKT \n(% of T cells)", xlab = "Null distribution")
abline(v=obsdiff_NKT, col="red",lwd=2)

Pval_NKT<-(sum(abs(avgdiff_NKT) > abs(obsdiff_NKT))+1) / (length(avgdiff_NKT)+1)

# iNKT
dat <- dplyr::filter(DAT, variable =="iNKT")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_iNKT<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_iNKT
set.seed(12345)
avgdiff_iNKT <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_iNKT, xlim = c(-abs(obsdiff_iNKT)-0.25, abs(obsdiff_iNKT)+0.25 ),
     main = "iNKT \n(% of T cells)", xlab = "Null distribution")
abline(v=obsdiff_iNKT, col="red",lwd=2)

Pval_iNKT<-(sum(abs(avgdiff_iNKT) > abs(obsdiff_iNKT))+1) / (length(avgdiff_iNKT)+1)


# CD4 memory #
dat <- dplyr::filter(DAT, variable =="CD4_memory")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_CD4memory<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_CD4memory
set.seed(12345)
avgdiff_CD4memory <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_CD4memory, xlim = c(-abs(obsdiff_CD4memory)-50, abs(obsdiff_CD4memory)+50 ),
     main = "CD4 memory \n(% of CD4)", xlab = "Null distribution")
abline(v=obsdiff_CD4memory, col="red",lwd=2)

Pval_CD4_memory<-(sum(abs(avgdiff_CD4memory) > abs(obsdiff_CD4memory))+1) / (length(avgdiff_CD4memory)+1)


# CD4 naive #
dat <- dplyr::filter(DAT, variable =="CD4_naive")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_CD4naive<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_CD4naive
set.seed(12345)
avgdiff_CD4naive <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_CD4naive, xlim = c(-abs(obsdiff_CD4naive)-50, abs(obsdiff_CD4naive)+50 ),
     main = "CD4 naive \n(% of CD4)", xlab = "Null distribution")
abline(v=obsdiff_CD4naive, col="red",lwd=2)

Pval_CD4_naive<-(sum(abs(avgdiff_CD4naive) > abs(obsdiff_CD4naive))+1) / (length(avgdiff_CD4naive)+1)

# CD4 EMRA

# CD4 EMRA #
dat <- dplyr::filter(DAT, variable =="CD4_EMRA")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_CD4EMRA<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_CD4EMRA
set.seed(12345)
avgdiff_CD4EMRA <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_CD4EMRA, xlim = c(-abs(obsdiff_CD4EMRA)-10, abs(obsdiff_CD4EMRA)+10 ),
     main = "CD4 EMRA \n(% of CD4)", xlab = "Null distribution")
abline(v=obsdiff_CD4EMRA, col="red",lwd=2)
Pval_CD4_EMRA<-(sum(abs(avgdiff_CD4EMRA) > abs(obsdiff_CD4EMRA))+1) / (length(avgdiff_CD4EMRA)+1)

# CD8 memory #
dat <- dplyr::filter(DAT, variable =="CD8_memory")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_CD8memory<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_CD8memory
set.seed(12345)
avgdiff_CD8memory <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_CD8memory, xlim = c(-abs(obsdiff_CD8memory)-50, abs(obsdiff_CD8memory)+50 ),
     main = "CD8 memory \n(% of CD8)", xlab = "Null distribution")
abline(v=obsdiff_CD8memory, col="red",lwd=2)

Pval_CD8_memory<-(sum(abs(avgdiff_CD8memory) > abs(obsdiff_CD8memory))+1) / (length(avgdiff_CD8memory)+1)


# CD8 naive #
dat <- dplyr::filter(DAT, variable =="CD8_naive")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_CD8naive<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_CD8naive
set.seed(12345)
avgdiff_CD8naive <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_CD8naive, xlim = c(-abs(obsdiff_CD8naive)-50, abs(obsdiff_CD8naive)+50 ),
     main = "CD8 naive \n(% of CD8)", xlab = "Null distribution")
abline(v=obsdiff_CD8naive, col="red",lwd=2)

Pval_CD8_naive<-(sum(abs(avgdiff_CD8naive) > abs(obsdiff_CD8naive))+1) / (length(avgdiff_CD8naive)+1)

# CD8 EMRA

# CD8 EMRA #
dat <- dplyr::filter(DAT, variable =="CD8_EMRA")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_CD8EMRA<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_CD8EMRA
set.seed(12345)
avgdiff_CD8EMRA <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_CD8EMRA, xlim = c(-abs(obsdiff_CD8EMRA)-20, abs(obsdiff_CD8EMRA)+20 ),
     main = "CD8 EMRA \n(% of CD8)", xlab = "Null distribution")
abline(v=obsdiff_CD8EMRA, col="red",lwd=2)

Pval_CD8_EMRA<-(sum(abs(avgdiff_CD8EMRA) > abs(obsdiff_CD8EMRA))+1) / (length(avgdiff_CD8EMRA)+1)

# Treg
dat <- dplyr::filter(DAT, variable =="Treg")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Treg<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Treg
set.seed(12345)
avgdiff_Treg <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Treg, xlim = c(-abs(obsdiff_Treg)-2, abs(obsdiff_Treg)+2 ),
     main = "Treg \n(% of CD4)", xlab = "Null distribution")
abline(v=obsdiff_Treg, col="red",lwd=2)

Pval_Treg<-(sum(abs(avgdiff_Treg) > abs(obsdiff_Treg))+1) / (length(avgdiff_Treg)+1)

# iTreg
dat <- dplyr::filter(DAT, variable =="iTreg")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_iTreg<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_iTreg
set.seed(12345)
avgdiff_iTreg <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_iTreg, xlim = c(-abs(obsdiff_iTreg)-0.2, abs(obsdiff_iTreg)+0.2 ),
     main = "iTreg \n(% of CD4)", xlab = "Null distribution")
abline(v=obsdiff_iTreg, col="red",lwd=2)

Pval_iTreg<-(sum(abs(avgdiff_iTreg) > abs(obsdiff_iTreg))+1) / (length(avgdiff_iTreg)+1)

# nTreg
dat <- dplyr::filter(DAT, variable =="nTreg")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_nTreg<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_nTreg
set.seed(12345)
avgdiff_nTreg <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_nTreg, xlim = c(-abs(obsdiff_nTreg)-0.2, abs(obsdiff_nTreg)+0.2 ),
     main = "nTreg \n(% of CD4)", xlab = "Null distribution")
abline(v=obsdiff_nTreg, col="red",lwd=2)

Pval_nTreg<-(sum(abs(avgdiff_nTreg) > abs(obsdiff_nTreg))+1) / (length(avgdiff_nTreg)+1)

# CD39Treg
dat <- dplyr::filter(DAT, variable =="CD39Treg")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_CD39Treg<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_CD39Treg
set.seed(12345)
avgdiff_CD39Treg <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_CD39Treg, xlim = c(-abs(obsdiff_CD39Treg)-0.5, abs(obsdiff_CD39Treg)+0.5 ),
     main = "CD39Treg \n(% of CD4)", xlab = "Null distribution")
abline(v=obsdiff_CD39Treg, col="red",lwd=2)

Pval_CD39Treg<-(sum(abs(avgdiff_CD39Treg) > abs(obsdiff_CD39Treg))+1) / (length(avgdiff_CD39Treg)+1)

# CD39nTreg
dat <- dplyr::filter(DAT, variable =="CD39nTreg")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_CD39nTreg<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_CD39nTreg
set.seed(12345)
avgdiff_CD39nTreg <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_CD39nTreg, xlim = c(-abs(obsdiff_CD39nTreg)-0.2, abs(obsdiff_CD39nTreg)+0.2 ),
     main = "CD39nTreg \n(% of CD4)", xlab = "Null distribution")
abline(v=obsdiff_CD39nTreg, col="red",lwd=2)

Pval_CD39nTreg<-(sum(abs(avgdiff_CD39nTreg) > abs(obsdiff_CD39nTreg))+1) / (length(avgdiff_CD39nTreg)+1)

# CD39iTreg
dat <- dplyr::filter(DAT, variable =="CD39iTreg")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_CD39iTreg<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_CD39iTreg
set.seed(12345)
avgdiff_CD39iTreg <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_CD39iTreg, xlim = c(-abs(obsdiff_CD39iTreg)-5, abs(obsdiff_CD39iTreg)+5 ),
     main = "CD39iTreg \n(% of CD4)", xlab = "Null distribution")
abline(v=obsdiff_CD39iTreg, col="red",lwd=2)

Pval_CD39iTreg<-(sum(abs(avgdiff_CD39iTreg) > abs(obsdiff_CD39iTreg))+1) / (length(avgdiff_CD39iTreg)+1)


# Monocyte
dat <- dplyr::filter(DAT, variable =="Monocyte")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Monocyte<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Monocyte
set.seed(12345)
avgdiff_Monocyte <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Monocyte, xlim = c(-abs(obsdiff_Monocyte)-20, abs(obsdiff_Monocyte)+20 ),
     main = "Monocyte \n(% of total)", xlab = "Null distribution")
abline(v=obsdiff_Monocyte, col="red",lwd=2)

Pval_Monocyte<-(sum(abs(avgdiff_Monocyte) > abs(obsdiff_Monocyte))+1) / (length(avgdiff_Monocyte)+1)


# Monocyte_Classical
dat <- dplyr::filter(DAT, variable =="Monocyte_Classical")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Monocyte_Classical<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Monocyte_Classical
set.seed(12345)
avgdiff_Monocyte_Classical <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Monocyte_Classical, xlim = c(-abs(obsdiff_Monocyte_Classical)-50, abs(obsdiff_Monocyte_Classical)+50 ),
     main = "Monocyte_Classical \n(% of monocytes)", xlab = "Null distribution")
abline(v=obsdiff_Monocyte_Classical, col="red",lwd=2)

Pval_Monocyte_Classical<-(sum(abs(avgdiff_Monocyte_Classical) > abs(obsdiff_Monocyte_Classical))+1) / (length(avgdiff_Monocyte_Classical)+1)

# Monocyte_Inflammatory
dat <- dplyr::filter(DAT, variable =="Monocyte_Inflammatory")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Monocyte_Inflammatory<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Monocyte_Inflammatory
set.seed(12345)
avgdiff_Monocyte_Inflammatory <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Monocyte_Inflammatory, xlim = c(-abs(obsdiff_Monocyte_Inflammatory)-20, abs(obsdiff_Monocyte_Inflammatory)+20 ),
     main = "Monocyte_Inflammatory \n(% of monocytes)", xlab = "Null distribution")
abline(v=obsdiff_Monocyte_Inflammatory, col="red",lwd=2)

Pval_Monocyte_Inflammatory<-(sum(abs(avgdiff_Monocyte_Inflammatory) > abs(obsdiff_Monocyte_Inflammatory))+1) / (length(avgdiff_Monocyte_Inflammatory)+1)

# Monocyte_Patrolling
dat <- dplyr::filter(DAT, variable =="Monocyte_Patrolling")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_Monocyte_Patrolling<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_Monocyte_Patrolling
set.seed(12345)
avgdiff_Monocyte_Patrolling <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_Monocyte_Patrolling, xlim = c(-abs(obsdiff_Monocyte_Patrolling)-20, abs(obsdiff_Monocyte_Patrolling)+20 ),
     main = "Monocyte_Patrolling \n(% of monocytes)", xlab = "Null distribution")
abline(v=obsdiff_Monocyte_Patrolling, col="red",lwd=2)

Pval_Monocyte_Patrolling<-(sum(abs(avgdiff_Monocyte_Patrolling) > abs(obsdiff_Monocyte_Patrolling))+1) / (length(avgdiff_Monocyte_Patrolling)+1)


# DC
dat <- dplyr::filter(DAT, variable =="DC")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_DC<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_DC
set.seed(12345)
avgdiff_DC <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_DC, xlim = c(-abs(obsdiff_DC)-0.1, abs(obsdiff_DC)+0.1 ),
     main = "DC \n(% of total)", xlab = "Null distribution")
abline(v=obsdiff_DC, col="red",lwd=2)

Pval_DC<-(sum(abs(avgdiff_DC) > abs(obsdiff_DC))+1) / (length(avgdiff_DC)+1)

# pDC
dat <- dplyr::filter(DAT, variable =="pDC")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_pDC<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_pDC
set.seed(12345)
avgdiff_pDC <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_pDC, xlim = c(-abs(obsdiff_pDC)-30, abs(obsdiff_pDC)+30 ),
     main = "pDC \n(% of total)", xlab = "Null distribution")
abline(v=obsdiff_pDC, col="red",lwd=2)

Pval_pDC<-(sum(abs(avgdiff_pDC) > abs(obsdiff_pDC))+1) / (length(avgdiff_pDC)+1)

# mDC
dat <- dplyr::filter(DAT, variable =="mDC")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_mDC<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_mDC
set.seed(12345)
avgdiff_mDC <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_mDC, xlim = c(-abs(obsdiff_mDC)-20, abs(obsdiff_mDC)+20 ),
     main = "mDC \n(% of total)", xlab = "Null distribution")
abline(v=obsdiff_mDC, col="red",lwd=2)

Pval_mDC<-(sum(abs(avgdiff_mDC) > abs(obsdiff_mDC))+1) / (length(avgdiff_mDC)+1)


# NK_cell
dat <- dplyr::filter(DAT, variable =="NK_cell")
dat <- dcast(dat, ID~Bleed)
A<-dat$a
D<-dat$d
G<-dat$g
J<-dat$j
intra.var<-apply(dat[-c(1:2)],1,var, na.rm=T)
CA<-combinations(length(A),4,A,set=FALSE)
CD<-combinations(length(D),4,D,set=FALSE)
CG<-combinations(length(G),4,G,set = FALSE)
CJ<-combinations(length(J),4,J,set=FALSE)
C<-rbind(CA,CD,CG,CJ)
head(C)
dim(C)
set.seed(12345)
C <- C[sample(nrow(C), 1000, replace = F),]
N = nrow(dat)
pop.var<-apply(C,1,var, na.rm=T)
set.seed(12345)
inter.var<-sample(pop.var, N)
obsdiff_NK_cell<-mean(inter.var, na.rm = T)-mean(intra.var, na.rm = T)
obsdiff_NK_cell
set.seed(12345)
avgdiff_NK_cell <- replicate(500, {
  all <- sample(pop.var)
  controls1 <- all[1:N]
  controls2 <- all[(N+1):(2*N)]
  return(mean(controls1, na.rm = T) - mean(controls2, na.rm = T))
})
hist(avgdiff_NK_cell, xlim = c(-abs(obsdiff_NK_cell)-20, abs(obsdiff_NK_cell)+20 ),
     main = "NK_cell \n(% of total)", xlab = "Null distribution")
abline(v=obsdiff_NK_cell, col="red",lwd=2)

Pval_NK_cell<-(sum(abs(avgdiff_NK_cell) > abs(obsdiff_NK_cell))+1) / (length(avgdiff_NK_cell)+1)

pvalues <- data.frame(Pval_Neutrophil, Pval_Lymphomonocyte, 
                      Pval_B_cell, Pval_B1, Pval_B_immature, Pval_B_memory, Pval_B_naive, 
                      Pval_Plasmablast,Pval_Plasmablast_IgA, Pval_T_cell, Pval_CD4, Pval_CD8, 
                      Pval_Gamma_delta_T, Pval_NKT, Pval_iNKT, 
                      Pval_CD4_naive, Pval_CD4_memory, Pval_CD4_EMRA,
                      Pval_CD8_naive, Pval_CD8_memory, Pval_CD8_EMRA,
                      Pval_Treg, Pval_iTreg, Pval_nTreg,
                      Pval_CD39Treg, Pval_CD39iTreg, Pval_CD39nTreg,
                      Pval_Monocyte, Pval_Monocyte_Classical, Pval_Monocyte_Inflammatory, Pval_Monocyte_Patrolling,
                      Pval_DC, Pval_mDC, Pval_pDC, Pval_NK_cell  )
pvalues <- data.frame(t(pvalues))
colnames(pvalues) <- "Pval"
pvalues$FDR <- p.adjust(pvalues$Pval, nrow(pvalues), method = "fdr")
pvalues$Subset <- rownames(pvalues)
rownames(pvalues)<- 1:nrow(pvalues)
pvalues$Subset <- gsub("Pval_","", pvalues$Subset)
pvalues$Subset <- gsub("_"," ", pvalues$Subset)
pvalues <- dplyr::select(pvalues, 3, 1, 2)
dev.off()
write.csv(pvalues, "output/TableS6.csv", row.names = F)



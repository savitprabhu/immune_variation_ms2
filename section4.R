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

########## SECTION 4 ##############

rm(list = ls())
dat = read.csv("output/TableS7.csv")
dat = dat[order(dat$Rho),]
dat$Subset = factor(dat$Subset, levels = dat$Subset)

p1 = ggplot(dat, aes(Subset, Rho))+
  geom_bar(stat = "identity", aes(fill = Significance), colour = "black")+
  theme_bw()+
  labs(x = NULL,
       y = "Correlation coefficient \n(mean versus variance of z-scores)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_y_continuous(limits = c(-1, 1))+
  theme(legend.position = c(0.8, 0.2))

p1

CV = read.csv("output/TableS4.csv")
CV$Name = gsub("_", " ", CV$Name)

dat$Subset = as.character(dat$Subset)
dat$Subset = gsub("Plasmablast IgA+", "Plasmablast IgA", dat$Subset)
dat[which(dat$Subset == "Plasmablast IgA+" ), "Subset"] = "Plasmablast IgA"
dat$Subset = gsub("Dendritic cell", "DC", dat$Subset)

identical(c(sort(CV$Name)), c(sort(dat$Subset)))

sort(CV$Name)            
sort(dat$Subset)

toplot = merge(CV, dat, by.x = "Name", by.y = "Subset")
cor.test(toplot$Rho, toplot$CV_betw_corrected, method = "spearman")

p2 = ggplot(toplot, aes(Rho, CV_betw_corrected))+
  geom_point(shape = 1)+
  annotate("text", 0, 0.75, label = "Rho = 0.6 \n p value = 0.0002")+
  labs(x = "Correlation coefficient \n(from mean-variance plots Fig 4A)",
       y = "Between-individual CV \n (after correction for technical CV)")+
  theme_bw()+
  coord_fixed(ratio = 1)
p2

pdf("output/Fig4.pdf", width = 5, height = 8)
plot_grid(p1, p2, labels = "AUTO", nrow = 2, rel_heights = c(1.5,1))
dev.off()

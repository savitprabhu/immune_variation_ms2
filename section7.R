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



########## SECTION 7  Sibling analysis ##############

rm(list = ls())

## Serial bleed  = unrelated individuals
SB <- read.csv("input/SB_frequency.csv")
SB$CD4_CD8_ratio <- SB$CD4/SB$CD8
SB <- dplyr::filter(SB, Bleed == "a")
SB <- SB[complete.cases(SB),] # remove rows with NA
set.seed(321)
SBa <- SB[sample(nrow(SB), 20), ] 
SBa$Pair_Alphabet <- rep("a")
SB <- merge(SB, SBa, all.x = T)
SB[which(is.na(SB$Pair_Alphabet)),"Pair_Alphabet"] <- rep("b")
SBa <- dplyr::filter(SB, Pair_Alphabet=="a")
SBb <- dplyr::filter(SB, Pair_Alphabet=="b")
SBa$ID <- paste("Pair", 1:nrow(SBa), sep = "_")
SBb$ID <- paste("Pair", 1:nrow(SBb), sep = "_")
SB <- rbind(SBa, SBb)
SB$Relation <- rep("Unrelated")
SB <- melt(SB)
SB <- dcast(SB, ID+variable~Pair_Alphabet)
SB$Difference <- abs(SB$a-SB$b)
SB$Relation <- rep("Unrelated")
str(SB)

# Siblings = related individuals 
Sib <- read.csv("input/Sib_frequency.csv")
Sib$CD4_CD8_ratio <- Sib$CD4/Sib$CD8
Sib$Pair_Alphabet <- tolower(Sib$Pair_Alphabet)
Sib <- dplyr::filter(Sib, Pair_Alphabet %in% c("a", "b"))
Sib <- melt(Sib)
Sib <- dcast(Sib, ID+variable~Pair_Alphabet)
Sib$Difference <- abs(Sib$a-Sib$b)
Sib$Relation <- rep("Related")
str(Sib)
identical(order(unique(SB$variable)), order(unique(Sib$variable)))

pdf("./output/FigS10.pdf", height = 5, width = 5)
## variable == CD4_CD8_ratio 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD4_CD8_ratio")
sibling <- dplyr::filter(Sib, variable == "CD4_CD8_ratio")
obsdiff_CD4_CD8_ratio <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD4_CD8_ratio <- (sum(abs(avgdiff) > abs(obsdiff_CD4_CD8_ratio)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD4_CD8_ratio", xlab = "Null distribution")
abline(v=obsdiff_CD4_CD8_ratio, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_CD4_CD8_ratio, 3)))

## variable == Neutrophil 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "Neutrophil")
sibling <- dplyr::filter(Sib, variable == "Neutrophil")
obsdiff_Neutrophil <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_Neutrophil <- (sum(abs(avgdiff) > abs(obsdiff_Neutrophil)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "Neutrophil", xlab = "Null distribution")
abline(v=obsdiff_Neutrophil, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_Neutrophil, 3)))

## variable == Lymphomonocyte 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "Lymphomonocyte")
sibling <- dplyr::filter(Sib, variable == "Lymphomonocyte")
obsdiff_Lymphomonocyte <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_Lymphomonocyte <- (sum(abs(avgdiff) > abs(obsdiff_Lymphomonocyte)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "Lymphomonocyte", xlab = "Null distribution")
abline(v=obsdiff_Lymphomonocyte, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_Lymphomonocyte, 3)))

## variable == B_cell 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "B_cell")
sibling <- dplyr::filter(Sib, variable == "B_cell")
obsdiff_B_cell <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_B_cell <- (sum(abs(avgdiff) > abs(obsdiff_B_cell)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "B_cell", xlab = "Null distribution")
abline(v=obsdiff_B_cell, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_B_cell, 3)))

## variable == B1 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "B1")
sibling <- dplyr::filter(Sib, variable == "B1")
obsdiff_B1 <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_B1 <- (sum(abs(avgdiff) > abs(obsdiff_B1)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "B1", xlab = "Null distribution")
abline(v=obsdiff_B1, col="red", lwd=2)
text(mean(avgdiff)+2*sd(avgdiff),200, paste("p = ", round(pval_B1, 3)))

## variable == B_immature 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "B_immature")
sibling <- dplyr::filter(Sib, variable == "B_immature")
obsdiff_B_immature <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_B_immature <- (sum(abs(avgdiff) > abs(obsdiff_B_immature)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "B_immature", xlab = "Null distribution")
abline(v=obsdiff_B_immature, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_B_immature, 3)))

## variable == B_memory 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "B_memory")
sibling <- dplyr::filter(Sib, variable == "B_memory")
obsdiff_B_memory <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_B_memory <- (sum(abs(avgdiff) > abs(obsdiff_B_memory)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "B_memory", xlab = "Null distribution")
abline(v=obsdiff_B_memory, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_B_memory, 3)))

## variable == B_naive 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "B_naive")
sibling <- dplyr::filter(Sib, variable == "B_naive")
obsdiff_B_naive <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_B_naive <- (sum(abs(avgdiff) > abs(obsdiff_B_naive)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "B_naive", xlab = "Null distribution")
abline(v=obsdiff_B_naive, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),150, paste("p = ", round(pval_B_naive, 3)))

## variable == Plasmablast 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "Plasmablast")
sibling <- dplyr::filter(Sib, variable == "Plasmablast")
obsdiff_Plasmablast <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_Plasmablast <- (sum(abs(avgdiff) > abs(obsdiff_Plasmablast)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "Plasmablast", xlab = "Null distribution")
abline(v=obsdiff_Plasmablast, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),150, paste("p = ", round(pval_Plasmablast, 3)))

## variable == Plasmablast_IgA 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "Plasmablast_IgA")
sibling <- dplyr::filter(Sib, variable == "Plasmablast_IgA")
obsdiff_Plasmablast_IgA <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_Plasmablast_IgA <- (sum(abs(avgdiff) > abs(obsdiff_Plasmablast_IgA)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "Plasmablast_IgA", xlab = "Null distribution")
abline(v=obsdiff_Plasmablast_IgA, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),120, paste("p = ", round(pval_Plasmablast_IgA, 3)))

## variable == T_cell 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "T_cell")
sibling <- dplyr::filter(Sib, variable == "T_cell")
obsdiff_T_cell <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_T_cell <- (sum(abs(avgdiff) > abs(obsdiff_T_cell)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "T_cell", xlab = "Null distribution")
abline(v=obsdiff_T_cell, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_T_cell, 3)))

## variable == CD4 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD4")
sibling <- dplyr::filter(Sib, variable == "CD4")
obsdiff_CD4 <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD4 <- (sum(abs(avgdiff) > abs(obsdiff_CD4)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD4", xlab = "Null distribution")
abline(v=obsdiff_CD4, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_CD4, 3)))

## variable == CD8
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD8")
sibling <- dplyr::filter(Sib, variable == "CD8")
obsdiff_CD8 <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD8 <- (sum(abs(avgdiff) > abs(obsdiff_CD8)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD8", xlab = "Null distribution")
abline(v=obsdiff_CD8, col="red", lwd=2)
text(mean(avgdiff)-3*sd(avgdiff),200, paste("p = ", round(pval_CD8, 3)))


## variable == Gamma_delta_T
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "Gamma_delta_T")
sibling <- dplyr::filter(Sib, variable == "Gamma_delta_T")
obsdiff_Gamma_delta_T <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_Gamma_delta_T <- (sum(abs(avgdiff) > abs(obsdiff_Gamma_delta_T)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "Gamma_delta_T", xlab = "Null distribution")
abline(v=obsdiff_Gamma_delta_T, col="red", lwd=2)
text(mean(avgdiff)+2*sd(avgdiff),200, paste("p = ", round(pval_Gamma_delta_T, 3)))

## variable == NKT 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "NKT")
sibling <- dplyr::filter(Sib, variable == "NKT")
obsdiff_NKT <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_NKT <- (sum(abs(avgdiff) > abs(obsdiff_NKT)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "NKT", xlab = "Null distribution")
abline(v=obsdiff_NKT, col="red", lwd=2)
text(mean(avgdiff)+2*sd(avgdiff),200, paste("p = ", round(pval_NKT, 3)))


## variable == iNKT
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "iNKT")
sibling <- dplyr::filter(Sib, variable == "iNKT")
obsdiff_iNKT <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_iNKT <- (sum(abs(avgdiff) > abs(obsdiff_iNKT)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "iNKT", xlab = "Null distribution")
abline(v=obsdiff_iNKT, col="red", lwd=2)
text(mean(avgdiff)+2*sd(avgdiff),120, paste("p = ", round(pval_iNKT, 3)))


## variable == CD4_naive 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD4_naive")
sibling <- dplyr::filter(Sib, variable == "CD4_naive")
obsdiff_CD4_naive <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD4_naive <- (sum(abs(avgdiff) > abs(obsdiff_CD4_naive)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD4_naive", xlab = "Null distribution")
abline(v=obsdiff_CD4_naive, col="red", lwd=2)
text(mean(avgdiff)-3*sd(avgdiff),200, paste("p = ", round(pval_CD4_naive, 3)))


## variable == CD4_memory 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD4_memory")
sibling <- dplyr::filter(Sib, variable == "CD4_memory")
obsdiff_CD4_memory <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD4_memory <- (sum(abs(avgdiff) > abs(obsdiff_CD4_memory)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD4_memory", xlab = "Null distribution")
abline(v=obsdiff_CD4_memory, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),150, paste("p = ", round(pval_CD4_memory, 3)))


## variable == CD4_EMRA 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD4_EMRA")
sibling <- dplyr::filter(Sib, variable == "CD4_EMRA")
obsdiff_CD4_EMRA <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD4_EMRA <- (sum(abs(avgdiff) > abs(obsdiff_CD4_EMRA)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD4_EMRA", xlab = "Null distribution")
abline(v=obsdiff_CD4_EMRA, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),120, paste("p = ", round(pval_CD4_EMRA, 3)))


## variable == CD8_naive 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD8_naive")
sibling <- dplyr::filter(Sib, variable == "CD8_naive")
obsdiff_CD8_naive <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD8_naive <- (sum(abs(avgdiff) > abs(obsdiff_CD8_naive)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD8_naive", xlab = "Null distribution")
abline(v=obsdiff_CD8_naive, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_CD8_naive, 3)))


## variable == CD8_memory 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD8_memory")
sibling <- dplyr::filter(Sib, variable == "CD8_memory")
obsdiff_CD8_memory <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD8_memory <- (sum(abs(avgdiff) > abs(obsdiff_CD8_memory)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD8_memory", xlab = "Null distribution")
abline(v=obsdiff_CD8_memory, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_CD8_memory, 3)))


## variable == CD8_EMRA 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD8_EMRA")
sibling <- dplyr::filter(Sib, variable == "CD8_EMRA")
obsdiff_CD8_EMRA <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD8_EMRA <- (sum(abs(avgdiff) > abs(obsdiff_CD8_EMRA)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD8_EMRA", xlab = "Null distribution", xlim = c(-10,10))
abline(v=obsdiff_CD8_EMRA, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),150, paste("p = ", round(pval_CD8_EMRA, 3)))

## variable == Treg 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "Treg")
sibling <- dplyr::filter(Sib, variable == "Treg")
obsdiff_Treg <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_Treg <- (sum(abs(avgdiff) > abs(obsdiff_Treg)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "Treg", xlab = "Null distribution")
abline(v=obsdiff_Treg, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),120, paste("p = ", round(pval_Treg, 3)))

## variable == iTreg 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "iTreg")
sibling <- dplyr::filter(Sib, variable == "iTreg")
obsdiff_iTreg <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_iTreg <- (sum(abs(avgdiff) > abs(obsdiff_iTreg)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "iTreg", xlab = "Null distribution")
abline(v=obsdiff_iTreg, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_iTreg, 3)))

## variable == nTreg 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "nTreg")
sibling <- dplyr::filter(Sib, variable == "nTreg")
obsdiff_nTreg <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_nTreg <- (sum(abs(avgdiff) > abs(obsdiff_nTreg)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "nTreg", xlab = "Null distribution")
abline(v=obsdiff_nTreg, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),120, paste("p = ", round(pval_nTreg, 3)))

## variable == CD39Treg 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD39Treg")
sibling <- dplyr::filter(Sib, variable == "CD39Treg")
obsdiff_CD39Treg <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD39Treg <- (sum(abs(avgdiff) > abs(obsdiff_CD39Treg)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD39Treg", xlab = "Null distribution")
abline(v=obsdiff_CD39Treg, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),120, paste("p = ", round(pval_CD39Treg, 3)))

## variable == CD39nTreg 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD39nTreg")
sibling <- dplyr::filter(Sib, variable == "CD39nTreg")
obsdiff_CD39nTreg <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD39nTreg <- (sum(abs(avgdiff) > abs(obsdiff_CD39nTreg)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD39nTreg", xlab = "Null distribution")
abline(v=obsdiff_CD39nTreg, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_CD39nTreg, 3)))

## variable == CD39iTreg 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "CD39iTreg")
sibling <- dplyr::filter(Sib, variable == "CD39iTreg")
obsdiff_CD39iTreg <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_CD39iTreg <- (sum(abs(avgdiff) > abs(obsdiff_CD39iTreg)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "CD39iTreg", xlab = "Null distribution", xlim = c(-2.5,2.5))
abline(v=obsdiff_CD39iTreg, col="red", lwd=2)
text(mean(avgdiff)-3*sd(avgdiff),120, paste("p = ", round(pval_CD39iTreg, 3)))

## variable == Monocyte 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "Monocyte")
sibling <- dplyr::filter(Sib, variable == "Monocyte")
obsdiff_Monocyte <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_Monocyte <- (sum(abs(avgdiff) > abs(obsdiff_Monocyte)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "Monocyte", xlab = "Null distribution")
abline(v=obsdiff_Monocyte, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),120, paste("p = ", round(pval_Monocyte, 3)))

## variable == Monocyte_Classical 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "Monocyte_Classical")
sibling <- dplyr::filter(Sib, variable == "Monocyte_Classical")
obsdiff_Monocyte_Classical <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_Monocyte_Classical <- (sum(abs(avgdiff) > abs(obsdiff_Monocyte_Classical)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "Monocyte_Classical", xlab = "Null distribution", xlim = c(-10, 10))
abline(v=obsdiff_Monocyte_Classical, col="red", lwd=2)
text(mean(avgdiff)-3*sd(avgdiff),120, paste("p = ", round(pval_Monocyte_Classical, 3)))

## variable == Monocyte_Patrolling 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "Monocyte_Patrolling")
sibling <- dplyr::filter(Sib, variable == "Monocyte_Patrolling")
obsdiff_Monocyte_Patrolling <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_Monocyte_Patrolling <- (sum(abs(avgdiff) > abs(obsdiff_Monocyte_Patrolling)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "Monocyte_Patrolling", xlab = "Null distribution")
abline(v=obsdiff_Monocyte_Patrolling, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),120, paste("p = ", round(pval_Monocyte_Patrolling, 3)))

## variable == Monocyte_Inflammatory 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "Monocyte_Inflammatory")
sibling <- dplyr::filter(Sib, variable == "Monocyte_Inflammatory")
obsdiff_Monocyte_Inflammatory <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_Monocyte_Inflammatory <- (sum(abs(avgdiff) > abs(obsdiff_Monocyte_Inflammatory)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "Monocyte_Inflammatory", xlab = "Null distribution")
abline(v=obsdiff_Monocyte_Inflammatory, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),150, paste("p = ", round(pval_Monocyte_Inflammatory, 3)))

## variable == DC 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "DC")
sibling <- dplyr::filter(Sib, variable == "DC")
obsdiff_DC <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_DC <- (sum(abs(avgdiff) > abs(obsdiff_DC)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "DC", xlab = "Null distribution", xlim = c(-abs(obsdiff_DC),abs(obsdiff_DC)))
abline(v=obsdiff_DC, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),150, paste("p = ", round(pval_DC, 3)))

## variable == mDC 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "mDC")
sibling <- dplyr::filter(Sib, variable == "mDC")
obsdiff_mDC <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_mDC <- (sum(abs(avgdiff) > abs(obsdiff_mDC)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "mDC", xlab = "Null distribution")
abline(v=obsdiff_mDC, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_mDC, 3)))

## variable == pDC 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "pDC")
sibling <- dplyr::filter(Sib, variable == "pDC")
obsdiff_pDC <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_pDC <- (sum(abs(avgdiff) > abs(obsdiff_pDC)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "pDC", xlab = "Null distribution")
abline(v=obsdiff_pDC, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_pDC, 3)))

## variable == NK_cell 
set.seed(123)
unrelated <- dplyr::filter(SB, variable == "NK_cell")
sibling <- dplyr::filter(Sib, variable == "NK_cell")
obsdiff_NK_cell <- mean(unrelated$Difference, na.rm = T) - mean(sibling$Difference, na.rm = T)
n <- nrow(unrelated) # no. of unrelated obs
N <- nrow(unrelated)+nrow(sibling) # total obs
avgdiff <- replicate(1000, {
  ALL <- sample(c(unrelated$Difference, sibling$Difference))
  new_unrelated <- ALL[1:n]
  new_sibling <- ALL[(n+1):(N)]
  return(mean(new_unrelated, na.rm = T) - mean(new_sibling, na.rm = T))
})
pval_NK_cell <- (sum(abs(avgdiff) > abs(obsdiff_NK_cell)) + 1) / (length(avgdiff) + 1)
hist(avgdiff, main = "NK_cell", xlab = "Null distribution")
abline(v=obsdiff_NK_cell, col="red", lwd=2)
text(mean(avgdiff)-2*sd(avgdiff),200, paste("p = ", round(pval_NK_cell, 3)))

dev.off()
pattern1 <-grep("pval",names(.GlobalEnv),value=TRUE)
pattern1 <-do.call("list",mget(pattern1))
pvalues <- t(data.frame(pattern1))
pvalues <- data.frame(rownames(pvalues), pvalues[1:nrow(pvalues)])
colnames(pvalues) <- c("Subset", "p-value")
pvalues$Subset <- substring(pvalues$Subset, 6)

pattern2<-grep("obsdiff",names(.GlobalEnv),value=TRUE)
pattern2 <-do.call("list",mget(pattern2))
obs_diff <- t(data.frame(pattern2))
obs_diff <- data.frame(rownames(obs_diff), obs_diff[1:nrow(obs_diff)])
colnames(obs_diff) <- c("Subset", "Observed_difference")
obs_diff$Subset <- substring(obs_diff$Subset, 9)

df <- merge(pvalues, obs_diff)
write.csv(df, "./output/TableS11.csv", row.names = F)

########## Boxplots of significant differences #########
Sib_SB <- rbind(Sib, SB)


## plot Treg
toplot <- dplyr::filter(Sib_SB, variable == "Treg")
p1<-ggplot(toplot, aes(Relation, Difference, fill=Relation))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=1, width = 0.05)+
  labs(x=NULL, y="Difference between pair", title="Treg")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(colour = "black", face = "bold"))+
  theme(text = element_text(color = "black", face = "bold"))

## plot CD8_naive
toplot <- dplyr::filter(Sib_SB, variable == "CD8_naive")
p2<-ggplot(toplot, aes(Relation, Difference, fill=Relation))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=1, width = 0.05)+
  labs(x=NULL, y="Difference between pair", title="CD8 naive")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(colour = "black", face = "bold"))+
  theme(text = element_text(color = "black", face = "bold"))

## plot CD4_naive
toplot <- dplyr::filter(Sib_SB, variable == "CD4_naive")
p3<-ggplot(toplot, aes(Relation, Difference, fill=Relation))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=1, width = 0.05)+
  labs(x=NULL, y="Difference between pair", title="CD4 naive")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(colour = "black", face = "bold"))+
  theme(text = element_text(color = "black", face = "bold"))

## plot CD4
toplot <- dplyr::filter(Sib_SB, variable == "CD4")
p4<-ggplot(toplot, aes(Relation, Difference, fill=Relation))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=1, width = 0.05)+
  labs(x=NULL, y="Difference between pair", title="CD4 T cell")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(colour = "black", face = "bold"))+
  theme(text = element_text(color = "black", face = "bold"))

## plot Monocyte_Classical
toplot <- dplyr::filter(Sib_SB, variable == "Monocyte_Classical")
p5<-ggplot(toplot, aes(Relation, Difference, fill=Relation))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=1, width = 0.05)+
  labs(x=NULL, y="Difference between pair", title="Classical Monocyte")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(colour = "black", face = "bold"))+
  theme(text = element_text(color = "black", face = "bold"))

pdf("./output/Fig8.pdf", width = 7, height = 7)
plot_grid(p1, p2, p3, p4, p5, nrow = 2)
dev.off()

## Plots of variables which did not show a difference
## plot PB
toplot <- dplyr::filter(Sib_SB, variable == "Plasmablast")
p6<-ggplot(toplot, aes(Relation, Difference, fill=Relation))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=1, width = 0.05)+
  labs(x=NULL, y="Difference between pair", title="Plasmablast")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(colour = "black", face = "bold"))+
  theme(text = element_text(color = "black", face = "bold"))

## Plot CD4 EMRA
toplot <- dplyr::filter(Sib_SB, variable == "CD4_EMRA")
p7<-ggplot(toplot, aes(Relation, Difference, fill=Relation))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=1, width = 0.05)+
  labs(x=NULL, y="Difference between pair", title="CD4 EMRA")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(colour = "black", face = "bold"))+
  theme(text = element_text(color = "black", face = "bold"))

## Plot Inflammatory monocyte
toplot <- dplyr::filter(Sib_SB, variable == "Monocyte_Inflammatory")
p8<-ggplot(toplot, aes(Relation, Difference, fill=Relation))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=1, width = 0.05)+
  labs(x=NULL, y="Difference between pair", title="Inflammatory monocyte")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(colour = "black", face = "bold"))+
  theme(text = element_text(color = "black", face = "bold"))

## Plot iNKT
toplot <- dplyr::filter(Sib_SB, variable == "iNKT")
p9<-ggplot(toplot, aes(Relation, Difference, fill=Relation))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=1, width = 0.05)+
  labs(x=NULL, y="Difference between pair", title="iNKT cell")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(colour = "black", face = "bold"))+
  theme(text = element_text(color = "black", face = "bold"))

## Plot B1 B cell
toplot <- dplyr::filter(Sib_SB, variable == "B1")
p10<-ggplot(toplot, aes(Relation, Difference, fill=Relation))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=1, width = 0.05)+
  labs(x=NULL, y="Difference between pair", title="B1 B cell")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(colour = "black", face = "bold"))+
  theme(text = element_text(color = "black", face = "bold"))

pdf("./output/FigS11.pdf", width = 7, height = 7)
plot_grid(p6, p7, p8, p9, p10, nrow = 2)
dev.off()




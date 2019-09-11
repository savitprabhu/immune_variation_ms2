# Set working directory to source destination
setwd("~/Documents/Immune_variation/3.Manuscripts/MS2/V5/analysis_V5")
library(readxl)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(scales)

rm(list = ls())

# Demographics 
Demog <- read.csv("input/SB_demographic.csv")
table(Demog$Sex)
mean(Demog$Age)
sd(Demog$Age)

######### Fig S1 Range of cell frequencies and counts ########
SBF <- read.csv("input/SB_frequency.csv")
SBC <- read.csv("input/SB_count.csv")
gates <- read.csv("input/Subsets_gates.csv")

#Plot of normal range of frequencies 
rm(list=setdiff(ls(), c("SBF", "SBC", "gates")))

dat <- SBF
datm <- melt(dat)
datm <- merge(datm, gates, by.x = "variable", by.y = "Name", all = T)
datm$Parent.Gate <- factor(datm$Parent.Gate, levels = c("Total WBC", "Total PBMC" ,"Monocyte", "Dendritic cell", "B cell", "T cell", "CD4 T", "CD8 T"))
datm <- datm[order(datm$Parent.Gate),]
datm$Subset <- factor(datm$Subset, levels = unique(datm$Subset))
datm <- filter(datm, Subset !="TLC")
mypalette <- c("white", "blue", "red", "green", "brown", "yellow", "purple", "grey")

p1<-ggplot(datm, aes(Subset, value, fill = Parent.Gate))+
  geom_boxplot(outlier.size = 0.25)+
  scale_fill_manual(values = mypalette)+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(text = element_text(colour = "black", face = "bold"))+
  theme(axis.text.x = element_text(angle = 90, 
                                   size = 8, #face = "bold",
                                   colour = "black",
                                   vjust = 0.5, hjust = 1))+
  #ggtitle("Normal range of immune cell frequencies in adult blood")+
  labs(x=NULL, y="Frequency of subset \n(as % of parent gate)", fill="Parent gate")

p1
arrangement <- c(as.character(unique(datm$Subset)))

# Plot of normal range of counts 
dat <- SBC
datm <- melt(dat)
datm <- merge(datm, gates, by.x = "variable", by.y = "Name", all.x  = T)
str(datm)
TLC <- "TLC"
arrangement <- c(TLC, arrangement)
datm$Subset <- factor(datm$Subset, levels = arrangement)

p2<- ggplot(datm, aes(Subset, value))+
  geom_boxplot(outlier.size = 0.25)+
  theme_bw()+
  theme(text = element_text(colour = "black", face = "bold"))+
  theme(axis.text.x = element_text(angle = 90, 
                                   size = 8, #face = "bold",
                                   colour = "black",
                                   vjust = 0.5, hjust = 1))+
  #ggtitle("Normal range of immune cell frequencies in adult blood")+
  labs(x=NULL, y="Counts")+
  scale_y_log10()
p2

pdf("output/FigS1.pdf", height = 10, width = 7)
plot_grid(p1, p2, nrow = 2, align = "v", labels = "AUTO")
dev.off()

# Spreadsheet of range of cell frequencies and counts 
SBF <- read.csv("input/SB_frequency.csv")
SBC <- read.csv("input/SB_count.csv")

######### Table S2: Table of normal ranges ##########
library(fBasics)
SUMMARY_FREQ <- data.frame(basicStats(SBF[,-c(1:3)]))[c("Mean", "Stdev", "Median", "Minimum", "Maximum", "1. Quartile", "3. Quartile"),]
SUMMARY_FREQ <- data.frame(t(SUMMARY_FREQ))
SUMMARY_FREQ$name <- rownames(SUMMARY_FREQ)
SUMMARY_FREQ <- melt(SUMMARY_FREQ)
SUMMARY_FREQ$Parameter <- rep("Frequency")
SUMMARY_COUNT <- data.frame(basicStats(SBC[,-c(1:3)]))[c("Mean", "Stdev", "Median", "Minimum", "Maximum", "1. Quartile", "3. Quartile"),]
SUMMARY_COUNT <- data.frame(t(SUMMARY_COUNT))
SUMMARY_COUNT$name <- rownames(SUMMARY_COUNT)
SUMMARY_COUNT <- melt(SUMMARY_COUNT)
SUMMARY_COUNT$Parameter <- rep("Count")
SUMMARY <- rbind(SUMMARY_FREQ, SUMMARY_COUNT)
SUMMARY <- dcast(SUMMARY, name+Parameter~variable, value.var = "value")
colnames(SUMMARY) 
colnames(SUMMARY) <- c("Subset","Parameter","Mean","Stdev","Median",
                       "Minimum", "Maximum","Quartile (25)", "Quartile (75)")

write.csv(SUMMARY, "output/TableS2.csv", row.names = F)

######## comparison between males and females ########
rm(list = ls())
Demog <- read.csv("input/SB_demographic.csv")
Demog = Demog[,1:2]
SBF <- read.csv("input/SB_frequency.csv")
SBC <- read.csv("input/SB_count.csv")
# insert demographic details into SBF and SBC
SBF = merge(Demog, SBF)
SBC = merge(Demog, SBC)
aggregate(SBF[,-c(1:4)], by=list(SBF$Sex), mean, na.rm = T)

# SBF
SBF = dplyr::filter(SBF, Bleed == "a") # choose only 1 time point for comparison
SBF_M = dplyr::filter(SBF, Sex == "M")
SBF_F = dplyr::filter(SBF, Sex == "F")

SBF_M = data.frame(basicStats(SBF_M[,-c(1:4)]))[c("Mean", "Stdev", "Median", "Minimum", "Maximum", "1. Quartile", "3. Quartile"),]
SBF_F = data.frame(basicStats(SBF_F[,-c(1:4)]))[c("Mean", "Stdev", "Median", "Minimum", "Maximum", "1. Quartile", "3. Quartile"),]

SBF_TTESTS <- lapply(SBF[,-c(1:4)], function(i) t.test(i ~ SBF$Sex))
SBF_TTESTS$CD4$p.value

SBF_PVALS <- data.frame(sapply(SBF_TTESTS, '[[', 'p.value'))
colnames(SBF_PVALS) <- c("P_value")
SBF_PVALS$P_adjust <- p.adjust(SBF_PVALS$P_value, method = "bonferroni", n = nrow(SBF_PVALS))
SBF_PVALS = data.frame(t(SBF_PVALS))

SBF_F = round(SBF_F, 1)
SBF_M = round(SBF_M, 1)
SBF_PVALS = round(SBF_PVALS, 3)
SBF_F$Sex = rep("Female")
SBF_M$Sex = rep("Male")
SBF_PVALS$Sex = rep("P-values")
SBF_F$Parameter = rownames(SBF_F)
SBF_M$Parameter = rownames(SBF_M)
SBF_PVALS$Parameter = rownames(SBF_PVALS)
SBF_summary = rbind(SBF_F, SBF_M, SBF_PVALS)
SBF_summary$VARIABLE = rep("Frequency")
SBF_summary = SBF_summary[,c(38,36, 37, 1:35)]

# SBC
SBC = dplyr::filter(SBC, Bleed == "a") # choose only 1 time point for comparison
SBC_M = dplyr::filter(SBC, Sex == "M")
SBC_F = dplyr::filter(SBC, Sex == "F")

SBC_M = data.frame(basicStats(SBC_M[,-c(1:4)]))[c("Mean", "Stdev", "Median", "Minimum", "Maximum", "1. Quartile", "3. Quartile"),]
SBC_F = data.frame(basicStats(SBC_F[,-c(1:4)]))[c("Mean", "Stdev", "Median", "Minimum", "Maximum", "1. Quartile", "3. Quartile"),]

SBC_TTESTS <- lapply(SBC[,-c(1:4)], function(i) t.test(i ~ SBC$Sex))
SBC_TTESTS$CD4$p.value

SBC_PVALS <- data.frame(sapply(SBC_TTESTS, '[[', 'p.value'))
colnames(SBC_PVALS) <- c("P_value")
SBC_PVALS$P_adjust <- p.adjust(SBC_PVALS$P_value, method = "bonferroni", n = nrow(SBC_PVALS))
SBC_PVALS = data.frame(t(SBC_PVALS))

SBC_F = round(SBC_F, 1)
SBC_M = round(SBC_M, 1)
SBC_PVALS = round(SBC_PVALS, 3)
SBC_F$Sex = rep("Female")
SBC_M$Sex = rep("Male")
SBC_PVALS$Sex = rep("P-values")
SBC_F$Parameter = rownames(SBC_F)
SBC_M$Parameter = rownames(SBC_M)
SBC_PVALS$Parameter = rownames(SBC_PVALS)
SBC_summary = rbind(SBC_F, SBC_M, SBC_PVALS)
SBC_summary$VARIABLE = rep("Count")
SBC_summary = SBC_summary[,c(39,37, 38, 1:36)]

SBF_summary$TLC = rep(NA)
GENDER = rbind(SBF_summary, SBC_summary)
write.csv(GENDER, "output/TableS3.csv", row.names = F)

########## Seasonal variation ###########
rm(list = ls())
DATES = read_excel("input/Dates_of_bleed.xlsx")
SBF <- read.csv("input/SB_frequency.csv")
DATES = melt(DATES)
names(DATES) = c("ID", "Bleed", "Month")
SBF = merge(DATES, SBF)
str(SBF$Month)
SBF$Parameter = NULL
SBF[,-c(1:3)] = scale(SBF[,-c(1:3)])

SBFm = melt(SBF, id.vars = c("ID", "Bleed", "Month"))

pdf("output/FigS2.pdf", width = 12, height = 12)
ggplot(SBFm, aes(Month, value))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm")+
  facet_wrap(~variable, nrow = 5)+
  scale_x_datetime(labels = date_format("%b"))+
  labs(y = "Frequency of subsets \n expressed as % of parent gate \n (z-scores)")+
  theme_bw()
dev.off()



########### Variance calculation ###########
#as suggested by Satyajit 11-09-2019 #

setwd("~/Documents/Immune_variation/3.Manuscripts/MS2/V5/analysis_V5")
library(readxl)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(scales)

rm(list = ls())

SBF <- read.csv("input/SB_frequency.csv")
SBF <- SBF[, !names(SBF) %in% c("TLC","Neutrophil", "Lymphomonocyte")] 

# between-individual variation
BETW = SBF
BETW = dplyr::filter(SBF, Bleed == "a")
BETW = BETW[,-c(1:3)]
BETW = apply(BETW, 2, function(x){x/mean(x, na.rm = T)}) 
# normalizing each value by dividing it by mean of the subset. 
# Now mean of all subsets is centered around 1.
# but the variance is preserved
apply(BETW, 1, median, na.rm = T) # there is some skewing in a few subsets

BETWVAR = data.frame(apply(BETW, 2, function(x){sd(x, na.rm = T)^2})) # calculate variance
BETWVAR$Name = rownames(BETWVAR)
rownames(BETWVAR) = 1:nrow(BETWVAR)
names(BETWVAR)[1] = "VAR_betw"
BETWVAR = BETWVAR[, c("Name", "VAR_betw")]
BETWVAR$Name = gsub("_", " ",BETWVAR$Name)
BETWVAR = BETWVAR[order(BETWVAR$VAR_betw),]
BETWVAR$Name = factor(BETWVAR$Name, levels = BETWVAR$Name)

# within-individual variation
INTRA = SBF
INTRA$Parameter = NULL
INTRA[,-c(1:2)] = apply(INTRA[,-c(1:2)], 2, function(x){x/mean(x, na.rm = T)})
INTRA = aggregate(INTRA[,-c(1:2)], by=list(INTRA$ID), function(x){sd(x, na.rm = T)^2})
colnames(INTRA)[1] = "ID"
INTRA = data.frame(apply(INTRA[,-1], 2, mean, na.rm=T))
INTRA$Name = rownames(INTRA)
rownames(INTRA) = 1:nrow(INTRA)
colnames(INTRA) = c("VAR_INTRA", "Name")
INTRA$Name = gsub("_", " ", INTRA$Name)
INTRA = INTRA[order(INTRA$VAR_INTRA),]
INTRA$Name = factor(INTRA$Name, levels = BETWVAR$Name)

VAR = merge(INTRA, BETWVAR)
VAR = VAR[order(VAR$VAR_betw),]
VAR$Name = factor(VAR$Name)
toplot = melt(VAR)
toplot$variable = gsub("VAR_INTRA", "Within-individual",toplot$variable)
toplot$variable = gsub("VAR_betw", "Between-individual",toplot$variable)

pdf("output/Fig1.pdf", width = 8, height = 7)
ggplot(toplot, aes(Name, value))+
  geom_bar(stat = "identity", aes(fill = variable), position = "dodge")+
  theme_cowplot()+
  labs(x = NULL,
       y = "Variance \n(normalized for mean)",
       fill = NULL)+
  theme(axis.text.x = element_text(angle = 90, size = 9, vjust = 0.5, hjust = 1))+
  theme(text = element_text(size = 11))+
  theme(legend.position = c(0.3, 0.8))
dev.off()

# technical variation
BB40 = read_excel("input/Rawdata_BB40.xlsx")
BB43 = read_excel("input/Rawdata_BB43.xlsx")
BB58 = read_excel("input/Rawdata_BB58.xlsx")

tech = rbind(BB40, BB43, BB58)

tech[,-1] = apply(tech[,-1], 1, function(x){x/mean(x, na.rm = T)}) # row-wise mean standardize
tech$Name = gsub("_", " ", tech$Name)
sort(VAR$Name)
sort(tech$Name)
tech$Name = factor(tech$Name, levels = VAR$Name)

tech = melt(tech)
tech$variable = NULL
tech = aggregate(tech$value, by = list(tech$Name), function(x){sd(x, na.rm = T)^2})
colnames(tech) = c("Name", "VAR_tech")
identical(sort(as.character(VAR$Name)), sort(as.character(tech$Name))) # subset names are identical

VAR = merge(VAR, tech)
VAR$Name = factor(VAR$Name, levels = BETWVAR$Name)

toplot = melt(VAR)
toplot$variable = gsub("VAR_INTRA", "Within-individual",toplot$variable)
toplot$variable = gsub("VAR_betw", "Between-individual",toplot$variable)
toplot$variable = gsub("VAR_tech", "Technical replicates",toplot$variable)
toplot$variable = factor(toplot$variable, levels = c("Between-individual",
                                                     "Within-individual",
                                                     "Technical replicates"))

p1 = ggplot(toplot, aes(Name, value))+
  geom_bar(stat = "identity", aes(fill = variable), position = "dodge")+
  theme_cowplot()+
  labs(x = NULL,
       y = "Variance \n(normalized for mean)",
       fill = NULL,
       title = "Technical variation in immune subsets")+
  theme(axis.text.x = element_text(angle = 90, size = 9, vjust = 0.5, hjust = 1))+
  theme(text = element_text(size = 11))+
  theme(legend.position = c(0.3, 0.8))
p1

summary(lm(data = VAR, VAR_betw~VAR_INTRA+VAR_tech))


p2=ggplot(VAR, aes(VAR_INTRA, VAR_betw))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_y_log10()+
  scale_x_log10()+
  annotate("text", x= 0.2, y = 0.03, label = "Adjusted R-squared:  0.52 \nt-value = 5.98 \np-value = 1.46e-06",
           colour = "red")+
  theme(text = element_text(size = 11))+
  labs(x = "Within-individual variance \nof immune subsets", 
       y = "Between-individual variance \nof immune subsets",
       title = "After adjusting for technical variation")
p2

summary(lm(data = VAR, VAR_INTRA~VAR_tech))

p3=ggplot(VAR, aes(VAR_INTRA, VAR_tech))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_y_log10()+
  scale_x_log10()+
  annotate("text", x= 0.5, y = 0.03, label = "p-value = 0.39",
           colour = "red")+
  theme(text = element_text(size = 11))+
  labs(x = "Within-individual variance \nof immune subsets", 
       y = "Variance in technical replicates \nof immune subsets",
       title = "Within-individual variance vs \ntechnical variance")
p3

summary(lm(data = VAR, VAR_betw~VAR_tech))

p4 = ggplot(VAR, aes(VAR_betw, VAR_tech))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_y_log10()+
  scale_x_log10()+
  annotate("text", x= 1.2, y = 0.02, label = "p-value = 0.37",
           colour = "red")+
  theme(text = element_text(size = 11))+
  labs(x = "Between-individual variance \nof immune subsets", 
       y = "Variance in technical replicates \nof immune subsets",
       title = "Between-individual variance vs \ntechnical variance")

p4


# check for mean-variance plot for technical replicates
BB40 = read_excel("input/Rawdata_BB40.xlsx")
BB43 = read_excel("input/Rawdata_BB43.xlsx")
BB58 = read_excel("input/Rawdata_BB58.xlsx")

BB40$ID = rep("BB40")
BB43$ID = rep("BB43")
BB58$ID = rep("BB58")
tech = rbind(BB40, BB43, BB58)

tech = melt(tech)
tech = aggregate(tech$value, by = list(tech$Name), mean, na.rm = T)
colnames(tech) = c("Name", "Mean") # these are population sizes

tech$Name = gsub("_", " ", tech$Name)
tech = merge(tech, VAR)
summary(lm(data = tech, VAR_tech~Mean))

p5=ggplot(tech, aes(Mean, VAR_tech))+
  geom_point()+
  geom_smooth(method = "lm", 
              colour = "red", size = 0.5)+
  annotate("text", x= 60, y = 0.6, label = "Adjusted R-squared:  0.07 \nt-value = 1.84 \np-value = 0.08",
           colour = "red")+
  theme_cowplot()+
  #scale_x_log10()+
  scale_y_log10()+
  theme(text = element_text(size = 11))+
  labs(x = "Population size of each subset\n(frequency expressed as % of parent gate)",
       y = "Variance in technical replicates",
       title = "Variability in technical replicates vs \n population size of immune subsets")
p5

VAR = VAR[, c(1, 3, 2, 4)]
VAR_NEW = VAR
colnames(VAR_NEW) = c("Name", "Between-individual", "Within-individual", "Technical replicates")
write.csv(VAR_NEW, "output/TableS4.csv", row.names = F)

pdf("output/FigS3.pdf", width = 7, height = 7)
p1
p2
p3
p4
p5
dev.off()

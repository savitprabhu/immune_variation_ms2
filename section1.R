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
rm(list = ls())
### technical CV calculation
SBF <- read.csv("input/SB_frequency.csv")
tech <- read.csv("input/technical_control_CV.csv")

SBF <- SBF[, !names(SBF) %in% c("TLC","Neutrophil", "Lymphomonocyte")] 
tech <- dplyr::filter(tech, !Name %in% c("TLC","Neutrophil", "Lymphomonocyte") )

BETW = SBF

#BETW = dplyr::filter(SBF, Bleed == "a")
BETW = BETW[,-c(1:3)]
BETWCV = data.frame(apply(BETW, 2, sd, na.rm = T)/apply(BETW, 2, mean, na.rm=T))
BETWCV$Name = rownames(BETWCV)
rownames(BETWCV) = 1:nrow(BETWCV)
names(BETWCV)[1] = "CV_betw"
BETWCV = BETWCV[, c("Name", "CV_betw")]
tech$Name = as.character(tech$Name)
identical(sort(BETWCV$Name), sort(tech$Name))
tech$CV_tech = apply(tech[,-1], 1, mean, na.rm= T)

CV = merge(tech, BETWCV)
CV$CV_betw_corrected = CV$CV_betw - CV$CV_tech

INTRA = SBF
INTRA = melt(SBF)
INTRA = dcast(INTRA, ID+variable~Bleed, value.var = "value")
str(INTRA)
INTRA$CV = apply(INTRA[,c(3:6)], 1, sd, na.rm = T) / apply(INTRA[,-c(1:2)], 1, mean, na.rm = T)

TECHCV = CV[, c("Name", "CV_tech")]
names(TECHCV) = c("variable", "CV_tech")
INTRA = merge(INTRA, TECHCV)
INTRA$CV_intra_corrected = INTRA$CV - INTRA$CV_tech
INTRACV = data.frame(aggregate(INTRA$CV_intra_corrected, by = list(INTRA$variable), FUN = mean, na.rm = T))
names(INTRACV) = c("Name", "CV_intra_corrected")
CV = merge(CV, INTRACV)

CV[which(CV$CV_betw_corrected < 0), "CV_betw_corrected"] = rep(0) # because CVs can't be negative

CV[which(CV$CV_intra_corrected < 0), "CV_intra_corrected"] = rep(0)

CV$BB40 = NULL
CV$BB58 = NULL
CV$BB43 = NULL
CV[,-1] = round(CV[,-1], 4)
write.csv(CV, "output/TableS4.csv", row.names = F)

############# Plot of Figure 1 ############
toplot = CV
toplot$Name = gsub("_", " ", toplot$Name)
toplot = toplot[order(toplot$CV_betw_corrected),]
str(toplot$Name)
toplot$Name = factor(toplot$Name, levels = toplot$Name)

p1 = ggplot(toplot, aes(Name, CV_betw_corrected))+
  geom_bar(stat="identity", fill = "pink", colour = "black")+
  theme_bw()+
  labs(x = NULL,
       y = "Between-individual CV \n(corrected for technical CV)")+
  theme(axis.text.x = element_text(angle = 90, size = 9, vjust = 0.5, hjust = 1))
p1
#dev.off()
toplot = CV
toplot$Name = gsub("_", " ", toplot$Name)
toplot = toplot[order(toplot$CV_intra_corrected),]
str(toplot$Name)
toplot$Name = factor(toplot$Name, levels = toplot$Name)

p2 = ggplot(toplot, aes(Name, CV_intra_corrected))+
  geom_bar(stat="identity", fill = "pink", colour = "black")+
  theme_bw()+
  labs(x = NULL,
       y = "Within-individual CV \n(corrected for technical CV)")+
  theme(axis.text.x = element_text(angle = 90, size = 9, vjust = 0.5, hjust = 1))
p2  

upper = plot_grid(p1, p2, nrow = 1, labels = "AUTO")

plot(CV$CV_betw_corrected, CV$CV_tech)
cor.test(CV$CV_betw_corrected, CV$CV_tech)

plot(CV$CV_intra_corrected, CV$CV_tech)
cor.test(CV$CV_intra_corrected, CV$CV_tech)

p1=ggplot(CV, aes(CV_tech, CV_betw_corrected))+
  geom_point()+
  theme_bw()+
  labs(x = "CV of technical controls",
       y = "Between-individual CV \n(corrected)")+
  annotate("text", 0.9, 0.8, label = "R = 0.22 \n p = 0.2")
p1

p2=ggplot(CV, aes(CV_tech, CV_intra_corrected))+
  geom_point()+
  theme_bw()+
  labs(x = "CV of technical controls",
       y = "Within-individual CV \n(corrected)")+
  annotate("text", 0.9, 0.2, label = "R = -0.3 \n p = 0.09")
p2

cor.test(CV$CV_intra_corrected, CV$CV_betw_corrected)
p3=ggplot(CV, aes(CV_betw_corrected, CV_intra_corrected))+
  geom_point()+
  theme_bw()+
  labs(x = "Between-individual CV \n(corrected)",
       y = "Within-individual CV \n(corrected)")+
  annotate("text", 0.2, 0.2, label = "R = 0.7 \n p = 9.4e-06")
p3

lower = plot_grid(p1, p2, p3, nrow = 1, labels = c("C", "D", "E"))

pdf("output/Fig1.pdf", width = 8.5, height = 7)
plot_grid(upper, lower, nrow = 2, rel_heights = c(1.5,1))
dev.off()



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
SBF <- read.csv("input/SB_frequency.csv")
SBF <- SBF[, !names(SBF) %in% c("TLC","Neutrophil", "Lymphomonocyte")] 

# between-individual variation
BETW = SBF
BETW = dplyr::filter(SBF, Bleed == "a")
BETW = BETW[,-c(1:3)]
BETWCV = data.frame(apply(BETW, 2, function(x){sd(x, na.rm = T)/mean(x, na.rm = T)}))
BETWCV$Name = rownames(BETWCV)
rownames(BETWCV) = 1:nrow(BETWCV)
names(BETWCV)[1] = "CV_betw"
BETWCV = BETWCV[, c("Name", "CV_betw")]
BETWCV$Name = gsub("_", " ",BETWCV$Name)
BETWCV = BETWCV[order(BETWCV$CV_betw),]
BETWCV$Name = factor(BETWCV$Name, levels = BETWCV$Name)

# within-individual variation
INTRA = SBF
INTRA$Parameter = NULL
INTRA = aggregate(INTRA[,-c(1:2)], by=list(INTRA$ID), function(x){sd(x, na.rm = T)/mean(x, na.rm = T)})
colnames(INTRA)[1] = "ID"
INTRA = data.frame(apply(INTRA[,-1], 2, mean, na.rm=T))
INTRA$Name = rownames(INTRA)
rownames(INTRA) = 1:nrow(INTRA)
colnames(INTRA) = c("CV_INTRA", "Name")
INTRA$Name = gsub("_", " ", INTRA$Name)
INTRA = INTRA[order(INTRA$CV_INTRA),]
INTRA$Name = factor(INTRA$Name, levels = BETWCV$Name)

CV = merge(INTRA, BETWCV)

# technical variation
tech = read_excel("input/BB40_raw_values.xlsx")
tech$Name = gsub("_", " ", tech$Name)
sort(CV$Name)
sort(tech$Name)
tech$Name = factor(tech$Name, levels = CV$Name)
tech$CV_tech = apply(tech[,-1], 1, function(x){sd(x, na.rm = T)/mean(x, na.rm = T)})
tech = tech[, c("Name", "CV_tech")]
identical(sort(as.character(CV$Name)), sort(as.character(tech$Name))) # subset names are identical

CV = merge(CV, tech)

CV$Name = factor(CV$Name, levels = BETWCV$Name)

p1 = ggplot(CV, aes(Name, CV_betw))+
  geom_bar(stat = "identity", colour = "black", fill = "#CC79A7")+
  theme_bw()+
  labs(x = NULL,
       y = "Between-individual variation \n(CV)")+
  theme(axis.text.x = element_text(angle = 90, size = 9, vjust = 0.5, hjust = 1))

p2 = ggplot(CV, aes(Name, CV_INTRA))+
  geom_bar(stat = "identity", colour = "black", fill = "#CC79A7")+
  theme_bw()+
  labs(x = NULL,
       y = "Within-individual variation \n(CV)")+
  theme(axis.text.x = element_text(angle = 90, size = 9, vjust = 0.5, hjust = 1))

library(cowplot)
pdf("output/Fig1.pdf", width = 8, height = 8)
plot_grid(p1, p2, nrow = 2, labels = "AUTO")
dev.off()

write.csv(CV, "output/TableS4.csv", row.names = F)

#

fit1 = lm(data = CV, CV_INTRA~CV_betw)
summary(fit1)
fit2 = lm(data = CV, CV_INTRA~CV_betw+CV_tech)
summary(fit2)

pdf("output/FigS3_fit.pdf")
summary(fit1)
dev.off()

p1 = ggplot(CV, aes(CV_betw, CV_INTRA))+
  geom_point(shape = 1, colour = "blue")+
  geom_smooth(method = "lm", colour = "red", size = 0.5)+
  annotate("text", 1.5, 0.2, label = "Adj R2 = 0.56 \nIntercept = 0.13 \nSlope = 0.254 \nP = 2.5e-06", colour = "red")+
  labs(x = "Between-individual variation (CV)",
       y = "Within-individual variation (CV)",
       title = "After adjusting for technical variability in the model")+
  theme_bw()

ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    annotate("text", 1.5, 0.2, label = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "\nIntercept =",signif(fit$coef[[1]],5 ),
                       " \nSlope =",signif(fit$coef[[2]], 5),
                       " \nP =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegression(fit2)+
  theme_bw()
ggplotRegression(fit2)

summary(fit2)

# check for mean-variance plot for technical replicates
tech = read_excel("input/BB40_raw_values.xlsx")
tech$MEAN_tech = apply(tech[,-1], 1, function(x){mean(x, na.rm = T)})
tech$CV_tech = apply(tech[,-1], 1, function(x){sd(x, na.rm = T)/mean(x, na.rm = T)})
tech$VAR_tech = apply(tech[,-1], 1, function(x){sd(x, na.rm = T)/mean(x, na.rm = T)})


plot(tech$MEAN_tech, tech$VAR_tech, 
     xlab = "Population size of each subset",
     ylab = "Technical variability",
     col = "blue",
     xlim = c(0, 100), 
     ylim = c(0, 1.5))
text(40, 1.2, "rho = -0.2 \np = 0.23", col = "red")
cor.test(tech$MEAN_tech, tech$VAR_tech, method = "spearman")

p2 = ggplot(tech, aes(MEAN_tech, VAR_tech))+
  geom_point(colour = "blue")+
  geom_smooth(method = "lm", colour = "red")+
  annotate("text", x= 40, y = 1.2, label = "rho = - 0.2 \nP = 0.23",
           colour = "red")+
  theme_bw()+
  labs(x = "Population size of each subset",
       y = "Variability in technical replicates",
       title = "Degree of variability in technical replicates \nas a function of population size of immune subsets")



pdf("output/FigS3.pdf", width = 5, height = 10)
plot_grid(p1, p2, nrow = 2, labels = "AUTO")
dev.off()



# rough work : different way of plotting #
toplot = melt(CV)
toplot$variable = gsub("CV_INTRA", "Within-individual variation", toplot$variable)
toplot$variable = gsub("CV_betw", "Between-individual variation", toplot$variable)
toplot$variable = gsub("CV_tech", "Technical-replicate variation", toplot$variable)
toplot$variable = factor(toplot$variable, 
                         levels = c("Between-individual variation",
                                    "Within-individual variation",
                                    "Technical-replicate variation"))
toplot = filter(toplot, variable != "Technical-replicate variation")

ggplot(toplot, aes(Name, value))+
  geom_bar(stat = "identity", colour = "black", fill = "#CC79A7")+
  theme_bw()+
  labs(x = NULL,
       y = "Coefficient of variation (CV)")+
  theme(axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5, hjust = 1))+
  facet_wrap(~variable, nrow = 1)+
  coord_flip()





# trial 
# removing technical CV outliers
CV = CV[which(CV$CV_tech < 1),]
fit1 = lm(data = CV, CV_INTRA~CV_tech)
summary(fit1)
plot(CV$CV_INTRA, CV$CV_tech)

fit2 = lm(data = CV, CV_INTRA~CV_betw+CV_tech)
summary(fit2)
plot(CV$CV_betw, CV$CV_INTRA)


fit3 = lm(data = CV, CV_betw~CV_tech)
summary(fit3)
plot(CV$CV_tech, CV$CV_betw)

fit4 = lm(data = CV, CV_betw~CV_INTRA+CV_tech)
summary(fit4)
plot(CV$CV_tech, CV$CV_INTRA)




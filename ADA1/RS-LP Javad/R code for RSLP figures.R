#January 2019####
#Zachary Bendiks 
#RSLP Paper - Barouei

rm(list = ls())

#1) Beta diversity plot showing the four treatment groups#### 
library(vegan)
library(MASS)

ADA_map <- read.table("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/mapping_ADA_all_groups.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)
ADA_dm <- read.table("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/wUniFrac distance matrix for RSLP Jan 19.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)

#needed to change the group names to match Javad's notation
ADA_map_updated_names <- ADA_map
ADA_map_updated_names$Treatment <- gsub("HF$", "CL", ADA_map_updated_names$Treatment)
ADA_map_updated_names$Treatment <- gsub("HF\\+RS", "RS", ADA_map_updated_names$Treatment)
ADA_map_updated_names$Treatment <- gsub("HF\\+LP", "LP", ADA_map_updated_names$Treatment)
ADA_map_updated_names$Treatment <- gsub("RS\\+LP", "RS\\+LP", ADA_map_updated_names$Treatment)

ADAdist <- as.dist(ADA_dm)
PCoA <- cmdscale(ADAdist)
str(PCoA)
ordiplot(PCoA)
plot(PCoA)

#create vector of sample names in the order they appear in the DM
sample_names <- rownames(ADA_dm)
#subset mapping file by sample names in DM
nrow(ADA_map_updated_names) #40
nrow(ADA_dm) #39 
ADAmap_sub <- ADA_map_updated_names[ADA_map_updated_names$SampleID %in% sample_names,]
nrow(ADAmap_sub) #39 (sample ADA08 missing from the DM)  
#reorder map object based on sample names in PCoA using the match argument 
ADAmap_sub_ordered <- ADAmap_sub[match(sample_names, ADAmap_sub$SampleID),]
ADAmap_sub_ordered$Treatment <- factor(ADAmap_sub_ordered$Treatment)

str(PCoA)
colnames(PCoA) <- c("PC1 - 63.43%", "PC2 - 11.94%")

#generating plot (note that the ^ means beginning of a string, while $ means the end.  The \\ symbol means to ignore the regex and treat it literally)
setwd("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/")
tiff("ADA RS LP PCoA Jan19.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 300)

par(mar=c(4,4,0.25,0.25))
plot(PCoA)

ordispider(PCoA[grep("^CL$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "CL"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^CL$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "grey80", cex = 1.2))

ordispider(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "black", cex = 1.2))

ordispider(PCoA[grep("RS\\+LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS+LP"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS\\+LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "black", cex = 1.2))

ordispider(PCoA[grep("^LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "LP"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "grey80", cex = 1.2))
#must use pt.bg instead of bg, or legend symbols will not fill! 
legend(x=0.1, y = -0.055, legend = c("CL", "RS", "LP", "RS+LP"), pch = c(21,21,22,22), pt.bg = c("grey80", "black", "grey80", "black"), col = "black", bty = "n", cex = 0.8, horiz = FALSE, ncol = 1, text.width = 0.15, xjust = 0)

dev.off()
#generating PDF version as well (higher quality, can easily convert to EPS file)
pdf("ADA RS LP PCoA.pdf", width = 4, height = 3.5)

par(mar=c(4,4,0.25,0.25))
plot(PCoA)

ordispider(PCoA[grep("^CL$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "CL"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^CL$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "grey80", cex = 1.2))

ordispider(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "black", cex = 1.2))

ordispider(PCoA[grep("RS\\+LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS+LP"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS\\+LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "black", cex = 1.2))

ordispider(PCoA[grep("^LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "LP"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "grey80", cex = 1.2))
#must use pt.bg instead of bg, or legend symbols will not fill! 
legend(x=0.1, y = -0.055, legend = c("CL", "RS", "LP", "RS+LP"), pch = c(21,21,22,22), pt.bg = c("grey80", "black", "grey80", "black"), col = "black", bty = "n", cex = 0.8, horiz = FALSE, ncol = 1, text.width = 0.15, xjust = 0)

dev.off()

#2) Alpha diversity plot of all 4 groups#### 
ADA_alphadiv <- read.table("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/faith PD for RSLP Jan19.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)

#install.packages("colorspace")
#install.packages("ggplot2")
library(ggplot2)
#install.packages("PMCMR")
library(PMCMR)
#install.packages("ggsignif")
library(ggsignif)

mean_HF <- mean(ADA_alphadiv$faith_pd[1:9])
SD_HF <- sd(ADA_alphadiv$faith_pd[1:9])

mean_RS <- mean(ADA_alphadiv$faith_pd[10:19])
SD_RS <- sd(ADA_alphadiv$faith_pd[10:18])

mean_LP <- mean(ADA_alphadiv$faith_pd[20:29])
SD_LP <- sd(ADA_alphadiv$faith_pd[20:29])

mean_RSLP <- mean(ADA_alphadiv$faith_pd[30:39])
SD_RSLP <- sd(ADA_alphadiv$faith_pd[30:39])

ADA_alphadiv$Group[30:39] <- rep("RS+LP", 10)
ADA_alphadiv$Group[1:9] <- rep("HF", 9)
ADA_alphadiv$Group[10:19] <- rep("RS", 10)
ADA_alphadiv$Group[20:29] <- rep("LP", 10)

#nonparametric (must make Group a factor before this works)
kruskal.test(data = ADA_alphadiv, faith_pd ~ as.factor(Group)) #Chi-squared = 11.156, df = 3, p = 0.01091 
posthoc.kruskal.dunn.test(x=ADA_alphadiv$faith_pd, g = as.factor(ADA_alphadiv$Group)) #HF v. RS and RS v. LP are the only pairwise differences  
#parametric
ADA_lm <- lm(data = ADA_alphadiv, formula = faith_pd~as.factor(Group))
ADA_anova_results <- anova(ADA_lm) #F = 4.59, df = 3, p = 0.0082
ADA_anova_model <- aov(ADA_lm) #model fitting
TukeyHSD(x = ADA_anova_model, conf.level = 0.95) #HF v. RS and RS v. LP are the only pairwise differences  

#changing the names of the groups and setting as factor for order rearrangement 
ADA_alphadiv_new_names <- ADA_alphadiv
ADA_alphadiv_new_names$Group <- gsub("HF$", "CL", ADA_alphadiv_new_names$Group)
ADA_alphadiv_new_names$Group <- factor(ADA_alphadiv_new_names$Group, levels = c("CL", "RS", "LP", "RS+LP")) #rearranges the x-axis order so RS appears second and Lp third 

#generating the plot 
par(mar=c(4.5,4.5,0.25,0.25))

ADA1_alpha_plot <- ggplot(ADA_alphadiv_new_names, aes(x=Group, y=faith_pd))+
  geom_boxplot()+
  ylim(7,17)+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 10, color = "black"), 
    axis.text.x = element_text(size=15, color = "black"), 
    axis.title.y = element_text(size = 15))+
  labs(x = NULL, y = "Faith PD at 41141 seqs")+
  geom_signif(comparisons = list(c("LP", "RS"), c("CL", "RS")), 
              map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), 
              test = wilcox.test,
              y_position = 16.5)

setwd("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019")
pdf("ADA RS LP alpha div plot.pdf", width = 4, height = 3.5)
ADA1_alpha_plot
dev.off()

tiff("ADA RS LP alpha div plot.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 300)
ADA1_alpha_plot
dev.off()


#3) L. plantarum abundance across all 4 groups#### 
#I manually removed "Constructed from biom file" line, changed "#OTU ID" to "OTU_ID", and deleted the "taxonomy" header before reading table in
ADA_table <- read.table("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/unannotated_rarefied_41141_table_edited.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)

#subsetting the table so only the ASVs of interest are represented 
ADA_table_Lp_ASVs_only <- subset(ADA_table, 
                                 ADA_table$OTU_ID == "79e4ffe038c224862d47da04769daecf"|
                                   ADA_table$OTU_ID == "4ff06c5727dbcdc7ca727c36d57d6a39"|
                                   ADA_table$OTU_ID == "4ddcd585962b7b0132f00958ccb19812"|
                                   ADA_table$OTU_ID == "39577a6f14ba9d3059bb6ef1a18185e2")

#dividing each value by the rarefaction depth to get abundance
ADA_table_Lp_ASVs_only_percent_abundance <- (ADA_table_Lp_ASVs_only[,2:ncol(ADA_table_Lp_ASVs_only)]/colSums(ADA_table[,2:ncol(ADA_table)]))*100
ADA_table_Lp_ASVs_only_percent_abundance <- data.frame(t(ADA_table_Lp_ASVs_only_percent_abundance))
#adds column showing the total abundance of Lp ASVs 
ADA_table_Lp_ASVs_only_percent_abundance$Total <- rowSums(ADA_table_Lp_ASVs_only_percent_abundance)
ADA_table_Lp_ASVs_only_percent_abundance$Group[30:39] <- rep("RS+LP", 10)
ADA_table_Lp_ASVs_only_percent_abundance$Group[1:9] <- rep("CL", 9)
ADA_table_Lp_ASVs_only_percent_abundance$Group[10:19] <- rep("RS", 10)
ADA_table_Lp_ASVs_only_percent_abundance$Group[20:29] <- rep("LP", 10)

ADA_table_Lp_ASVs_only_percent_abundance$Group <- factor(ADA_table_Lp_ASVs_only_percent_abundance$Group, levels = c("CL", "RS", "LP", "RS+LP"))
#creating axis label expression that can be fed into "lab" argument in ggplot
y_axis <- expression(paste(italic("L. plantarum"), " abundance (%)")) 

#nonparametric test
kruskal.test(data = ADA_table_Lp_ASVs_only_percent_abundance, Total ~ as.factor(Group)) #chi-squared = 29.545, df = 3, p-value = 1.72e-06 
posthoc.kruskal.dunn.test(x=ADA_table_Lp_ASVs_only_percent_abundance$Total, g = as.factor(ADA_table_Lp_ASVs_only_percent_abundance$Group)) #CL-LP, CL-RSLP, LP-RS, RS-RSLP (Holm adjusted p values)
#parametric
ADA_Lacto_lm <- lm(data = ADA_table_Lp_ASVs_only_percent_abundance, formula = Total~as.factor(Group))
ADA__Lacto_anova_results <- anova(ADA_Lacto_lm) #F = 14.69, df = 3, p = 2.35e-6
ADA_Lacto_anova_model <- aov(ADA_Lacto_lm) #model fitting
TukeyHSD(x = ADA_Lacto_anova_model, conf.level = 0.95) #LP-CL, RS+LP-CL, LP-RS, RS+LP-RS

ADA_RSLP_Lacto_plot <- ggplot(ADA_table_Lp_ASVs_only_percent_abundance, aes(x=Group, y=Total))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 13))+
  labs(x = NULL, y = y_axis)+ #adjusting axis titles, using formatted text expression for Y axis 
  geom_signif(comparisons = list(c("CL", "LP")), map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), test = wilcox.test, y_position = 3.8, tip_length = 0)+ #adding significance notation
  geom_signif(comparisons = list(c("CL", "RS+LP")), map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), test = wilcox.test, y_position = 4.1, tip_length = 0)+
  geom_signif(comparisons = list(c("RS", "LP")), map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), test = wilcox.test, y_position = 3.5, tip_length = 0)+
  geom_signif(comparisons = list(c("RS", "RS+LP")), map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), test = wilcox.test, y_position = 3.2, tip_length = 0)+
  ylim(0, 4.2)


setwd("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/")
pdf("ADA RS LP abundance plot.pdf", width = 4, height = 3.5)
ADA_RSLP_Lacto_plot
dev.off()

tiff("ADA RS LP abundance plot.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 300)
ADA_RSLP_Lacto_plot
dev.off()

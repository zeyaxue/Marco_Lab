#May 2019####
#Zachary Bendiks 
#RSLP Paper (ADA 1)- Barouei

#generating figures to summarize microbiome results for manuscript.  

rm(list = ls())
set.seed(1)

#1) Beta diversity plot showing the four treatment groups#### 
library(vegan)
library(MASS)

ADA_map <- read.table("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/mapping_ADA_all_groups.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)
ADA_dm <- read.table("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/wUniFrac distance matrix for RSLP May 19.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)

#needed to change the group names to match Javad's notation
ADA_map_updated_names <- ADA_map
#ADA_map_updated_names$Treatment <- gsub("HF$", "CL", ADA_map_updated_names$Treatment)
ADA_map_updated_names$Treatment <- gsub("HF\\+RS", "RS", ADA_map_updated_names$Treatment)
ADA_map_updated_names$Treatment <- gsub("HF\\+LP", "LP", ADA_map_updated_names$Treatment)
ADA_map_updated_names$Treatment <- gsub("RS\\+LP", "RS\\+LP", ADA_map_updated_names$Treatment)
rm(ADA_map)

ADAdist <- as.dist(ADA_dm)
PCoA <- cmdscale(ADAdist)
str(PCoA)
ordiplot(PCoA)
plot(PCoA)

#create vector of sample names in the order they appear in the DM
sample_names <- rownames(ADA_dm)
#subset mapping file by sample names in DM
nrow(ADA_map_updated_names) #40
nrow(ADA_dm) #40 (39 if ADA08 filtered out, but wasn't in this analysis)
#ADAmap_sub <- ADA_map_updated_names[ADA_map_updated_names$SampleID %in% sample_names,]
  
#reorder map object based on sample names in PCoA using the match argument 
ADAmap_sub_ordered <- ADA_map_updated_names[match(sample_names, ADA_map_updated_names$SampleID),]
ADAmap_sub_ordered$Treatment <- factor(ADAmap_sub_ordered$Treatment)

str(PCoA)
colnames(PCoA) <- c("PC1 - 72.36%", "PC2 - 7.85%")

#generating plot (note that the ^ means beginning of a string, while $ means the end.  The \\ symbol means to ignore the regex and treat it literally)
setwd("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/")
tiff("ADA RS LP PCoA May19.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 300)

par(mar=c(4,4,0.25,0.25))
plot(PCoA)

ordispider(PCoA[grep("^HF$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "HF"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^HF$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "grey80", cex = 1.2))

ordispider(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "black", cex = 1.2))

ordispider(PCoA[grep("RS\\+LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS+LP"], col = "black", lwd = 1.75)
ordispider(PCoA[grep("^LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "LP"], col = "grey80", lwd = 1.75)

with(ADAmap_sub_ordered, points(PCoA[grep("RS\\+LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "black", cex = 1.2))
with(ADAmap_sub_ordered, points(PCoA[grep("^LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "grey80", cex = 1.2))
#must use pt.bg instead of bg, or legend symbols will not fill! 
legend(x=-0.245, 
       y = 0.15, 
       legend = c("HF", "RS", "LP", "RS+LP"), 
       pch = c(21,21,22,22), 
       pt.bg = c("grey80", "black", "grey80", "black"), 
       pt.cex = 1.25,
       col = "black", 
       bty = "n", 
       cex = 0.9, 
       horiz = FALSE, 
       ncol = 1, 
       text.width = 0.15, 
       xjust = 0)

dev.off()

#generating PDF version as well (higher quality, can easily convert to EPS file)
pdf("ADA RS LP PCoA.pdf", width = 4, height = 3.5)

par(mar=c(4,4,0.25,0.25))
plot(PCoA)

ordispider(PCoA[grep("^HF$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "HF"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^HF$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "grey80", cex = 1.2))

ordispider(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "black", cex = 1.2))

ordispider(PCoA[grep("RS\\+LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS+LP"], col = "black", lwd = 1.75)
ordispider(PCoA[grep("^LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "LP"], col = "grey80", lwd = 1.75)

with(ADAmap_sub_ordered, points(PCoA[grep("RS\\+LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "black", cex = 1.2))
with(ADAmap_sub_ordered, points(PCoA[grep("^LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "grey80", cex = 1.2))
#must use pt.bg instead of bg, or legend symbols will not fill! 
legend(x=-0.245, 
       y = 0.15, 
       legend = c("HF", "RS", "LP", "RS+LP"), 
       pch = c(21,21,22,22), 
       pt.bg = c("grey80", "black", "grey80", "black"), 
       pt.cex = 1.25,
       col = "black", 
       bty = "n", 
       cex = 0.9, 
       horiz = FALSE, 
       ncol = 1, 
       text.width = 0.15, 
       xjust = 0)

dev.off()

#2) Alpha diversity plot of all 4 groups#### 
#download from the rarefaction curve visualization from QIIME2 
ADA_alphadiv <- read.table("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/faith PD for RSLP May19.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)

#install.packages("colorspace")
#install.packages("ggplot2")
library(ggplot2)
#install.packages("PMCMR")
library(PMCMR)
#install.packages("ggsignif")
library(ggsignif)

mean_HF <- mean(ADA_alphadiv$faith_pd[1:10])
SD_HF <- sd(ADA_alphadiv$faith_pd[1:10])

mean_RS <- mean(ADA_alphadiv$faith_pd[11:20])
SD_RS <- sd(ADA_alphadiv$faith_pd[11:20])

mean_LP <- mean(ADA_alphadiv$faith_pd[21:30])
SD_LP <- sd(ADA_alphadiv$faith_pd[21:30])

mean_RSLP <- mean(ADA_alphadiv$faith_pd[31:40])
SD_RSLP <- sd(ADA_alphadiv$faith_pd[31:40])

ADA_alphadiv$Group[31:40] <- rep("RS+LP", 10)
ADA_alphadiv$Group[1:10] <- rep("HF", 10)
ADA_alphadiv$Group[11:20] <- rep("RS", 10)
ADA_alphadiv$Group[21:30] <- rep("LP", 10)

#nonparametric (must make Group a factor before this works)
kruskal.test(data = ADA_alphadiv, faith_pd ~ as.factor(Group)) #Chi-squared = 12.476, df = 3, p = 0.0059 
posthoc.kruskal.dunn.test(x=ADA_alphadiv$faith_pd, g = as.factor(ADA_alphadiv$Group)) #HF v. RS and RS v. LP are the only pairwise differences  
#parametric testing group only 
replications(data = ADA_alphadiv, formula = faith_pd~Group)
ADA_lm_group <- lm(data = ADA_alphadiv, formula = faith_pd~Group)
ADA_anova_results_group <- anova(ADA_lm_group) #F = 2.25, df = 3, p = 0.09959
ADA_anova_model_group <- aov(ADA_lm_group) #model fitting
TukeyHSD(x = ADA_anova_model_group, conf.level = 0.95) #no significant difference  


#2.1) Alpha diversity interaction effects####
#go back to group object and add columns for LP and RS factors.
ADA_alphadiv$Lp <- c(rep("no LP", 20), rep("LP", 20))
ADA_alphadiv$Lp <- as.factor(ADA_alphadiv$Lp)
ADA_alphadiv$Starch <- c(rep("no RS", 10), rep("RS", 10), rep("no RS", 10), rep("RS", 10))
ADA_alphadiv$Starch <- as.factor(ADA_alphadiv$Starch)

#parametric two-way ANOVA 
ADA_lm_2way_faith <- lm(data = ADA_alphadiv, formula = faith_pd ~ Starch * Lp)
ADA_anova_results_2way_faith <- anova(ADA_lm_2way_faith) #Starch was significant, but no interaction effects between the two 
ADA_anova_model_2way_faith <- aov(ADA_lm_2way_faith) #model fitting
TukeyHSD(x = ADA_anova_model_2way_faith, conf.level = 0.95) 

#changing the names of the groups and setting as factor for order rearrangement 
ADA_alphadiv_new_names <- ADA_alphadiv
#ADA_alphadiv_new_names$Group <- gsub("HF$", "CL", ADA_alphadiv_new_names$Group)
ADA_alphadiv_new_names$Group <- factor(ADA_alphadiv_new_names$Group, levels = c("HF", "RS", "LP", "RS+LP")) #rearranges the x-axis order so RS appears second and Lp third 

#2.2) Generating plot####
par(mar=c(4.5,4.5,0.25,0.25))

ADA1_alpha_plot <- ggplot(ADA_alphadiv_new_names, aes(x=Group, y=faith_pd))+
  geom_boxplot()+
  ylim(7,17)+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 9, color = "black"), 
    axis.text.x = element_text(size=12, color = "black"), 
    axis.title.y = element_text(size = 12))+
  labs(x = NULL, y = "Faith PD at 33,925 seqs")+
  geom_signif(comparisons = list(c("LP", "RS")), 
              map_signif_level = c("p = .078"=0.05), #based on Tukey's test above 
              test = t.test,
              y_position = 16.5)

setwd("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019")
pdf("ADA RS LP alpha div plot.pdf", width = 3, height = 3)
ADA1_alpha_plot
dev.off()

tiff("ADA RS LP alpha div plot.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 300)
ADA1_alpha_plot
dev.off()

#2.3) Shannon diversity index ####
ADA_alphadiv_shannon <- read.table("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/shannon index for RSLP Jan19.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)

#install.packages("colorspace")
#install.packages("ggplot2")
#library(ggplot2)
#install.packages("PMCMR")
#library(PMCMR)
#install.packages("ggsignif")
#library(ggsignif)

ADA_alphadiv_shannon$Group[31:40] <- rep("RS+LP", 10)
ADA_alphadiv_shannon$Group[1:10] <- rep("HF", 10)
ADA_alphadiv_shannon$Group[11:20] <- rep("RS", 10)
ADA_alphadiv_shannon$Group[21:30] <- rep("LP", 10)

ADA_alphadiv_shannon$Lp <- c(rep("no LP", 20), rep("LP", 20))
ADA_alphadiv_shannon$Lp <- as.factor(ADA_alphadiv_shannon$Lp)
ADA_alphadiv_shannon$Starch <- c(rep("no RS", 10), rep("RS", 10), rep("no RS", 10), rep("RS", 10))
ADA_alphadiv_shannon$Starch <- as.factor(ADA_alphadiv_shannon$Starch)

#nonparametric (must make Group a factor before this works)
kruskal.test(data = ADA_alphadiv_shannon, shannon ~ as.factor(Group)) #chi-squared = 22.658, df = 3, p-value = 4.758e-05
posthoc.kruskal.dunn.test(x=ADA_alphadiv_shannon$shannon, g = as.factor(ADA_alphadiv_shannon$Group)) 
#testing interaction effects 
ADA_lm_shannon <- lm(data = ADA_alphadiv_shannon, formula = shannon ~ Starch * Lp)
ADA_anova_results_shannon <- anova(ADA_lm_shannon) #Starch was significant, LP almost significant, but no interaction effects between the two 
ADA_anova_model_shannon <- aov(ADA_lm_shannon) #model fitting
TukeyHSD(x = ADA_anova_model_shannon, conf.level = 0.95) #HF v. RS and RS v. LP are the only pairwise differences  

#rearranging table 
ADA_alphadiv_shannon_new_names <- ADA_alphadiv_shannon
ADA_alphadiv_shannon_new_names$Group <- factor(ADA_alphadiv_shannon_new_names$Group, levels = c("HF", "RS", "LP", "RS+LP")) #rearranges the x-axis order so RS appears second and Lp third 

par(mar=c(4.5,4.5,0.25,0.25))

ADA1_alpha_plot_shannon <- ggplot(ADA_alphadiv_shannon_new_names, aes(x=Group, y=shannon))+
  geom_boxplot()+
  ylim(3,6.5)+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 9, color = "black"), 
    axis.text.x = element_text(size=12, color = "black"), 
    axis.title.y = element_text(size = 12))+
  labs(x = NULL, y = "Shannon Index at 33,925 seqs")+
  geom_signif(comparisons = list(c("HF", "RS"), c("RS", "LP"), c("LP", "RS+LP")), 
              map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), 
              test = wilcox.test,
              y_position = 6.1)+
  geom_signif(comparisons = list(c("HF", "RS+LP")), 
              map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), 
              test = wilcox.test,
              y_position = 6.4)

setwd("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019")
pdf("ADA RS LP alpha div plot shannon.pdf", width = 3, height = 3)
ADA1_alpha_plot_shannon
dev.off()

tiff("ADA RS LP alpha div plot shannon.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 300)
ADA1_alpha_plot_shannon
dev.off()

#3) L. plantarum abundance across all 4 groups#### 
#I manually removed "Constructed from biom file" line, changed "#OTU ID" to "OTU_ID", and deleted the "taxonomy" header before reading table in
ADA_table <- read.table("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/unannotated_rarefied_33925_table_for_plot.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)

#converting counts to relative abundance 
ADA_table_percent_abundance <- ADA_table[,2:ncol(ADA_table)]/colSums(ADA_table[,2:ncol(ADA_table)])*100
rownames(ADA_table_percent_abundance) <- ADA_table[,1]
#subsetting the table so only the ASVs of interest are represented 
ADA_table_Lp_ASVs_only_percent_abundance <- subset(ADA_table_percent_abundance, 
                                                   rownames(ADA_table_percent_abundance) == "4ff06c5727dbcdc7ca727c36d57d6a39"|
                                                  rownames(ADA_table_percent_abundance) == "39577a6f14ba9d3059bb6ef1a18185e2")

tADA_table_Lp_ASVs_only_percent_abundance <- data.frame(t(ADA_table_Lp_ASVs_only_percent_abundance))

#adds column showing the total abundance of Lp ASVs 
tADA_table_Lp_ASVs_only_percent_abundance$Total <- rowSums(tADA_table_Lp_ASVs_only_percent_abundance)
tADA_table_Lp_ASVs_only_percent_abundance$Group[31:40] <- rep("RS+LP", 10)
tADA_table_Lp_ASVs_only_percent_abundance$Group[1:10] <- rep("HF", 10)
tADA_table_Lp_ASVs_only_percent_abundance$Group[11:20] <- rep("RS", 10)
tADA_table_Lp_ASVs_only_percent_abundance$Group[21:30] <- rep("LP", 10)

tADA_table_Lp_ASVs_only_percent_abundance$Group <- factor(tADA_table_Lp_ASVs_only_percent_abundance$Group, levels = c("HF", "RS", "LP", "RS+LP"))
#creating axis label expression that can be fed into "lab" argument in ggplot
y_axis <- expression(paste(italic("L. plantarum"), " abundance (%)")) 


#nonparametric test
kruskal.test(data = tADA_table_Lp_ASVs_only_percent_abundance, Total ~ Group) #chi-squared = 30.973, df = 3, p-value = 8.611e-7 
posthoc.kruskal.dunn.test(x=tADA_table_Lp_ASVs_only_percent_abundance$Total, g = tADA_table_Lp_ASVs_only_percent_abundance$Group) 
#parametric
ADA_Lacto_lm <- lm(data = tADA_table_Lp_ASVs_only_percent_abundance, formula = Total~Group)
ADA__Lacto_anova_results <- anova(ADA_Lacto_lm) #F = 15.195, df = 3, p = 1.499e-6
ADA_Lacto_anova_model <- aov(ADA_Lacto_lm) #model fitting
TukeyHSD(x = ADA_Lacto_anova_model, conf.level = 0.95) #LP-CL, RS+LP-CL, LP-RS, RS+LP-RS

tADA_table_Lp_ASVs_only_percent_abundance[,3] <- signif(tADA_table_Lp_ASVs_only_percent_abundance[,3], 6)

ADA_RSLP_Lacto_plot <- ggplot(tADA_table_Lp_ASVs_only_percent_abundance, aes(x=Group, y=Total))+
  geom_boxplot()+
  geom_point(size = 2)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 9, color = "black"), 
        axis.text.x = element_text(size=12, color = "black"), 
        axis.title.y = element_text(size = 12))+
  labs(x = NULL, y = y_axis)+ #adjusting axis titles, using formatted text expression for Y axis 
  geom_signif(comparisons = list(c("HF", "LP")), map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), test = wilcox.test, y_position = 3.85, tip_length = 0)+ #adding significance notation
  geom_signif(comparisons = list(c("HF", "RS+LP")), map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), test = wilcox.test, y_position = 4.2, tip_length = 0)+
  geom_signif(comparisons = list(c("RS", "LP")), map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), test = wilcox.test, y_position = 3.15, tip_length = 0)+
  geom_signif(comparisons = list(c("RS", "RS+LP")), map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), test = wilcox.test, y_position = 3.5, tip_length = 0)+
  ylim(0, 4.3)


setwd("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/")
pdf("ADA RS LP abundance plot.pdf", width = 3, height = 3)
ADA_RSLP_Lacto_plot
dev.off()

tiff("ADA RS LP abundance plot.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 300)
ADA_RSLP_Lacto_plot
dev.off()


#3.1) Lp abundance interaction effects####
tADA_table_Lp_ASVs_only_percent_abundance$Lp <- c(rep("no LP", 20), rep("LP", 20))
tADA_table_Lp_ASVs_only_percent_abundance$Lp <- as.factor(tADA_table_Lp_ASVs_only_percent_abundance$Lp)
tADA_table_Lp_ASVs_only_percent_abundance$Starch <- c(rep("no RS", 10), rep("RS", 10), rep("no RS", 10), rep("RS", 10))
tADA_table_Lp_ASVs_only_percent_abundance$Starch <- as.factor(tADA_table_Lp_ASVs_only_percent_abundance$Starch)

#parametric ANOVA
ADA_lm_lacto <- lm(data = tADA_table_Lp_ASVs_only_percent_abundance, formula = Total ~ Starch * Lp)
ADA_anova_results_lacto <- anova(ADA_lm_lacto) #Lp was significant, RS was not, and no interaction effects between the two 
ADA_anova_model_lacto <- aov(ADA_lm_lacto) #model fitting
TukeyHSD(x = ADA_anova_model_lacto, conf.level = 0.95) #HF v. RS and RS v. LP are the only pairwise differences  


#4) Bar chart taxa summary ####
#need to convert the table to long format for ggplot2 
library(ggplot2)
#install.packages("reshape2", dependencies = TRUE)
#install.packages("stringi", dependencies = TRUE)
library(reshape2)
#remove #constructed from biom, replace #OTU ID with OTU_ID, and delete taxonomy column header before import 
#REMOVE ANY SPACES 
ADA_taxa_table <- read.table("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/reanalysis/Jan 2019/rarefied_table_collapsed_L6_for_taxa_summary_edited.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)
rownames(ADA_taxa_table) <- ADA_taxa_table[,1]
ADA_taxa_table <- ADA_taxa_table[,-1] #putting taxa strings as rownames then deleting OTU_ID column 

#retaining just the top 10 taxa based on total sums 
ADA_taxa_sums <- data.frame(rowSums(ADA_taxa_table))
ADA_taxa_sums_order <- order(ADA_taxa_sums, decreasing = TRUE)
ADA_taxa_sums[ADA_taxa_sums_order,]
ADA_taxa_table_order <- ADA_taxa_table[ADA_taxa_sums_order, ,drop = FALSE] #setting 100,000 TMM reads as cutoff will remove most phyla

top10_taxa <- ADA_taxa_table_order[1:10, , drop=FALSE]
top10_taxa <- row.names(top10_taxa)
low_taxa <- ADA_taxa_table_order[11:nrow(ADA_taxa_table_order), , drop=FALSE]
low_taxa <- row.names(low_taxa)

#aggregate low abundance phyla together and add back to the table 
low_taxa_table <- ADA_taxa_table_order[row.names(ADA_taxa_table_order) %in% low_taxa,]
low_taxa_table_sums <- colSums(low_taxa_table)

top10_taxa_table <- ADA_taxa_table_order[row.names(ADA_taxa_table_order) %in% top10_taxa,]
top10_taxa_table <- rbind(top10_taxa_table, low_taxa_table_sums)
row.names(top10_taxa_table)[11] <- "Other"

#converting coutns to relative abundance
colSums(top10_taxa_table) #still all at 33925 
str(top10_taxa_table)
top10_taxa_table_abun <- prop.table(as.matrix(top10_taxa_table), 2) #doesn't work on data frames for whatever reason...
colSums(top10_taxa_table_abun) #should all = 1 

#cleaning up taxa strings 
taxa_names <- rownames(top10_taxa_table_abun)
taxa_names <- t(data.frame(strsplit(taxa_names, ";"))) #split string by ;, convert to data frame, and transpose 
taxa_names <- gsub(".__", "", taxa_names) #removes all taxonomic level designations 
taxa_names <- paste(taxa_names[,4], taxa_names[,5], taxa_names[,6], sep = "|") #combine order through genus levels, separated by "|"
taxa_names <- gsub("__", "", taxa_names) #cleaning up some of the strings 
taxa_names[11] <- "Other" #couldn't figure out how to negate regex for | so I jsut changed the name manually 
#adding cleaned up taxa strings back to the original table 
rownames(top10_taxa_table_abun) <- taxa_names
#decided to not include the "Other" group and instead just leave the area blank 
top10_taxa_table_abun_no_other <- top10_taxa_table_abun[-11,]
#changing the sample names to 1-10 so ggplot doesn't add a bunch of unneeded columns when I plot with facet_grid
colnames(top10_taxa_table_abun_no_other) <- c(rep(1:10, 4))
#want to do PERCENT abundance and not use decimals, so multiply all values by 100 
top10_taxa_table_abun_no_other_percent <- top10_taxa_table_abun_no_other*100

#melting table to get it in long format 
top10_taxa_table_abun_melt <- melt(top10_taxa_table_abun_no_other_percent)
top10_taxa_table_abun_melt$group <- c(rep("HF", 100), rep("RS", 100), rep("LP", 100), rep("RS+LP",100))
#convert group column to factor and specify levels so the facet plots are made in that order 
top10_taxa_table_abun_melt$group <- factor(top10_taxa_table_abun_melt$group, levels = c("HF", "RS", "LP", "RS+LP"))
names(top10_taxa_table_abun_melt) <- c("Taxa","sample","abundance", "group")


#generating plot with a diverging, manual color palette (easier to see differences across groups)
taxa_abundance_barplot <- ggplot(data = top10_taxa_table_abun_melt,
       aes(as.factor(sample), y = abundance, fill = Taxa)) +
  geom_bar(stat = "identity", color = "black")+
  facet_grid(~group)+
  ylab("Relative Abundance (%)") + xlab("Samples")+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 15))+
  scale_fill_manual(values = c("forestgreen", "yellow", "blue", "red", "purple", "limegreen", "hotpink", "brown", "orange", "darkturquoise"))


#generate plot, save to file 
tiff("ADA_RSLP_Taxa_barchart.tiff", width = 12, height = 3.5, units = "in", compression = "lzw", res = 600)
taxa_abundance_barplot
dev.off()

pdf("ADA_RSLP_Taxa_barchart.pdf", width = 12, height = 3.5)
taxa_abundance_barplot
dev.off()


#5) Combining all figures into one multi-panel plot ####
#cannot use par functions or 'layout' to combine base and ggplot2 plots.  Need a way to capture ggplot2 objects  
#install.packages("cowplot")
#install.packages("gridGraphics") #required to use cowplot 
library(cowplot)
library(gridExtra)
#library(gridGraphics)
#capturing each plot that I need for multipanel with recordPlot()
#first plot (Beta diversity) 
par(mar=c(4,4,0.25,0.25))
plot(PCoA)

ordispider(PCoA[grep("^HF$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "HF"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^HF$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "grey80", cex = 1.2))

ordispider(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "black", cex = 1.2))

ordispider(PCoA[grep("RS\\+LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS+LP"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS\\+LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "black", cex = 1.2))

ordispider(PCoA[grep("^LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "LP"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^LP", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "grey80", cex = 1.2))
#must use pt.bg instead of bg, or legend symbols will not fill! 
legend(x=-0.245, 
       y = 0.15, 
       legend = c("HF", "RS", "LP", "RS+LP"), 
       pch = c(21,21,22,22), 
       pt.bg = c("grey80", "black", "grey80", "black"), 
       pt.cex = 1.25,
       col = "black", 
       bty = "n", 
       cex = 0.8, 
       horiz = FALSE, 
       ncol = 1, 
       text.width = 0.15, 
       xjust = 0)
plot(PCoA)
plot1 <- grid.grab() 

#second plot (alpha diversity) 
#already saved as ggplot object, so just have to recall and record it
ADA1_alpha_plot
plot2 <- grid.grab()

#third plot (LP abundance) 
#already saved as ggplot object, so just have to recall and record it
ADA_RSLP_Lacto_plot
plot3 <- grid.grab()

#fourth plot (taxa summary barchart) 
#already saved as ggplot object, so just have to recall and record it
taxa_abundance_barplot
plot4 <- grid.grab()

#can set up a 3x3 object with the plot positions specified like so 
panel_layout <- rbind(c(1,2,3), 4)
layout(1)
#layout(1) #use to make jsut one plot, not a multipanel
#par(mar=c(3,4,1.5,0.1)) #adjust plot margins 

grid.arrange(plot1, plot2, plot3, plot4)
dev.off()


#for some reason the vegan ordination plot won't save as a grob.  Fucking annoying
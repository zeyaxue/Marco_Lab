#ADA1 RS-LP QIIME2 analysis
#Zach Bendiks
#July 17, 2018 

rm(list = ls())

#1) Beta diversity plot showing the four treatment groups#### 
library(vegan)
library(MASS)

ADA_map <- read.table("E:/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/mapping_ADA_all_groups.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)
ADA_dm <- read.table("E:/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/distance-matrix.tsv", sep = "", stringsAsFactors = FALSE, header = TRUE)

#needed to change the group names to match Javad's notation
ADA_map_updated_names <- ADA_map
ADA_map_updated_names$Treatment <- gsub("HF$", "CL", ADA_map_updated_names$Treatment)
ADA_map_updated_names$Treatment <- gsub("HF\\+RS", "RS", ADA_map_updated_names$Treatment)
ADA_map_updated_names$Treatment <- gsub("HF\\+LP", "Lp", ADA_map_updated_names$Treatment)
ADA_map_updated_names$Treatment <- gsub("RS\\+LP", "RS\\+Lp", ADA_map_updated_names$Treatment)

ADAdist <- as.dist(ADA_dm)
PCoA <- cmdscale(ADAdist)
str(PCoA)
ordiplot(PCoA)
plot(PCoA)

#create vector of sample names in the order they appear in the DM
sample_names <- rownames(ADA_dm)
#subset mapping file by sample names in DM
nrow(ADA_map_updated_names) #40
nrow(ADA_dm) #38 
ADAmap_sub <- ADA_map_updated_names[ADA_map_updated_names$SampleID %in% sample_names,]
nrow(ADAmap_sub) #38 (samples ADA08 and ADA16 missing from the DM)  
#reorder map object based on sample names in PCoA using the match argument 
ADAmap_sub_ordered <- ADAmap_sub[match(sample_names, ADAmap_sub$SampleID),]
ADAmap_sub_ordered$Treatment <- factor(ADAmap_sub_ordered$Treatment)

str(PCoA)
colnames(PCoA) <- c("PC1 - 59.14%", "PC2 - 11.82%")

#generating plot (note that the ^ means beginning of a string, while $ means the end.  The \\ symbol means to ignore the regex and treat it literally)
setwd("E:/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/")
tiff("ADA RS LP PCoA.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 300)

plot(PCoA)

ordispider(PCoA[grep("^CL$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "CL"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^CL$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "grey80", cex = 1.2))

ordispider(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "black", cex = 1.2))

ordispider(PCoA[grep("RS\\+Lp", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS+Lp"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS\\+Lp", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "black", cex = 1.2))

ordispider(PCoA[grep("^Lp", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "Lp"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^Lp", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "grey80", cex = 1.2))

legend(x=0.055, y = 0.095, legend = c("CL", "RS", "Lp", "RS+Lp"), pch = c(19,19,15,15), col = c("grey80", "black", "grey80", "black"), bty = "n", cex = 1.0, horiz = FALSE, ncol = 1, text.width = 0.15, xjust = 0)

dev.off()
#generating PDF version as well (higher quality, can easily convert to EPS file)
pdf("ADA RS LP PCoA.pdf", width = 5, height = 5)

plot(PCoA)

ordispider(PCoA[grep("^CL$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "CL"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^CL$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "grey80", cex = 1.2))

ordispider(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS$", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 21, bg = "black", cex = 1.2))

ordispider(PCoA[grep("RS\\+Lp", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "RS+Lp"], col = "black", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("RS\\+Lp", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "black", cex = 1.2))

ordispider(PCoA[grep("^Lp", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], group = ADAmap_sub_ordered$Treatment[ADAmap_sub_ordered$Treatment == "Lp"], col = "grey80", lwd = 1.75)
with(ADAmap_sub_ordered, points(PCoA[grep("^Lp", ADAmap_sub_ordered$Treatment), 1:2, drop = FALSE], col = "black", pch = 22, bg = "grey80", cex = 1.2))

legend(x=0.055, y = 0.095, legend = c("CL", "RS", "Lp", "RS+Lp"), pch = c(19,19,15,15), col = c("grey80", "black", "grey80", "black"), bty = "n", cex = 1.0, horiz = FALSE, ncol = 1, text.width = 0.15, xjust = 0)

dev.off()




#2) Alpha diversity plot of all 4 groups#### 
ADA_alphadiv <- read.table("E:/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/alpha-diversity.tsv", sep = "", stringsAsFactors = FALSE, header = TRUE)

#install.packages("colorspace")
#install.packages("ggplot2")
library(ggplot2)
#install.packages("PMCMR")
library(PMCMR)
#install.packages("ggsignif")
library(ggsignif)

mean_HF <- mean(ADA_alphadiv$faith_pd[1:9])
SD_HF <- sd(ADA_alphadiv$faith_pd[1:9])

mean_RS <- mean(ADA_alphadiv$faith_pd[10:18])
SD_RS <- sd(ADA_alphadiv$faith_pd[10:18])

mean_LP <- mean(ADA_alphadiv$faith_pd[19:28])
SD_LP <- sd(ADA_alphadiv$faith_pd[19:28])

mean_RSLP <- mean(ADA_alphadiv$faith_pd[29:38])
SD_RSLP <- sd(ADA_alphadiv$faith_pd[29:38])

ADA_alphadiv$Group[29:38] <- rep("HF+RS+LP", 10)
ADA_alphadiv$Group[1:9] <- rep("HF", 9)
ADA_alphadiv$Group[10:18] <- rep("HF+RS", 9)
ADA_alphadiv$Group[19:28] <- rep("HF+LP", 10)

#nonparametric
kruskal.test(data = ADA_alphadiv, faith_pd ~ as.factor(Group)) #p = 0.04162, must make group a factor before KW test will work 
posthoc.kruskal.dunn.test(x=ADA_alphadiv$faith_pd, g = as.factor(ADA_alphadiv$Group)) #only sig. difference is between HF+LP and HF+RS 
#parametric
ADA_lm <- lm(data = ADA_alphadiv, formula = faith_pd~as.factor(Group))
ADA_anova_results <- anova(ADA_lm) #p=0.01886
ADA_anova_model <- aov(ADA_lm) #model fitting
TukeyHSD(x = ADA_anova_model, conf.level = 0.95) #also shows that the only sig difference is between HF+RS and HF+LP

#changing the names of the groups and setting as factor for order rearrangement 
ADA_alphadiv_new_names <- ADA_alphadiv
ADA_alphadiv_new_names$Group <- gsub("HF$", "CL", ADA_alphadiv_new_names$Group)
ADA_alphadiv_new_names$Group <- gsub("HF\\+RS", "RS", ADA_alphadiv_new_names$Group)
ADA_alphadiv_new_names$Group <- gsub("HF\\+LP", "Lp", ADA_alphadiv_new_names$Group)
ADA_alphadiv_new_names$Group <- gsub("RS\\+LP", "RS\\+Lp", ADA_alphadiv_new_names$Group)
ADA_alphadiv_new_names$Group <- factor(ADA_alphadiv_new_names$Group, levels = c("CL", "RS", "Lp", "RS+Lp")) #rearranges the x-axis order so RS appears second and Lp third 

ADA1_alpha_plot <- ggplot(ADA_alphadiv_new_names, aes(x=Group, y=faith_pd))+
  geom_boxplot()+
  ylim(10,16)+
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(size=15, color = "black"), 
    axis.title.y = element_text(size = 15))+
  labs(x = NULL, y = "Faith PD")+
  geom_signif(comparisons = list(c("Lp", "RS")), map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), test = wilcox.test)

setwd("G:/My Drive/TRANSCEND/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/")
pdf("ADA RS LP alpha div plot.pdf", width = 3, height = 4)
ADA1_alpha_plot
dev.off()




#3) Lactobacillus abundance across all 4 groups#### 
#had to remove the # symbols and delete all spaces.  I also moved the taxa strings over and manually added OTU numbers 
ADA_table <- read.table("E:/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/ADA_RSLP_feature-table_headers_edited.txt", sep = "", stringsAsFactors = FALSE, header = TRUE)

ADA_table_Lacto_only <- ADA_table[grep("Lactobacillus", ADA_table$taxonomy),]
ADA_table_Lacto_only_scaled <- ADA_table_Lacto_only[,2:39]/colSums(ADA_table[,2:39]) #converting values to relative abundance 
ADA_table_Lacto_only_scaled <- t(ADA_table_Lacto_only_scaled)

ADA_table_Lacto_only_scaled <- data.frame(ADA_table_Lacto_only_scaled)
ADA_table_Lacto_only_scaled$group[29:38] <- rep("RS+Lp", 10)
ADA_table_Lacto_only_scaled$group[1:9] <- rep("CL", 9)
ADA_table_Lacto_only_scaled$group[10:18] <- rep("RS", 9)
ADA_table_Lacto_only_scaled$group[19:28] <- rep("Lp", 10)
ADA_table_Lacto_only_scaled$group <- factor(ADA_table_Lacto_only_scaled$group, levels = c("CL", "RS", "Lp", "RS+Lp"))
#creating axis label expression that can be fed into "lab" argument in ggplot
y_axis <- expression(paste(italic("Lactobacillus"), " abundance (rarefied to 35,738 seqs)")) 

#nonparametric test
kruskal.test(data = ADA_table_Lacto_only_scaled, X29 ~ as.factor(group)) #p = 0.01165, must make group a factor before KW test will work 
posthoc.kruskal.dunn.test(x=ADA_table_Lacto_only_scaled$X29, g = as.factor(ADA_table_Lacto_only_scaled$group)) #no pairwise differences between the groups
#parametric
ADA_Lacto_lm <- lm(data = ADA_table_Lacto_only_scaled, formula = X29~as.factor(group))
ADA__Lacto_anova_results <- anova(ADA_Lacto_lm) #p=0.005042
ADA_Lacto_anova_model <- aov(ADA_Lacto_lm) #model fitting
TukeyHSD(x = ADA_Lacto_anova_model, conf.level = 0.95) #LP-RS and RSLP-LP are the only pairwise comparisons that differ.  RS seems to inhibit LActo expansion

ADA_RSLP_Lacto_plot <- ggplot(ADA_table_Lacto_only_scaled, aes(x=group, y=X29))+
     geom_boxplot()+
     theme_bw()+
     theme(axis.text.x = element_text(size = 12),
           axis.text.y = element_text(size = 12),
           axis.title = element_text(size = 13))+
     labs(x = "Group", y = y_axis)+ #adjusting axis titles, using formatted text expression for Y axis 
     geom_signif(comparisons = list(c("Lp", "RS")), map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), test = wilcox.test)+ #adding significance notation
     geom_signif(comparisons = list(c("Lp", "RS+Lp")), map_signif_level = c("*"=0.001, "*"=0.01, "*"=0.05), test = wilcox.test)+
     ylim(0.015, 0.0375)

setwd("E:/ADA1/RS+LP Manunscript (Javad)/RS_LP/QIIME2/")
pdf("ADA RS LP lactobacillus abundance plot.pdf", width = 4.5, height = 4.5)
ADA_RSLP_Lacto_plot
dev.off()




#4) Differential abundance analysis####
ap <- available.packages()
View(ap)
"yingtools2" %in% rownames(ap)
#not available in CRAN, had to download GitHub repository using devtools package 
install.packages("devtools")
library(devtools)
devtools::install_github("ying14/yingtools2")
library(yingtools2)
#loading additional packages 
library(dplyr)
library(phyloseq)
library(ggplot2)
library(data.table)

lefse2 <- function (phy, class, subclass = NA, subject = NA, anova.alpha = 0.05, 
                    wilcoxon.alpha = 0.05, lda.cutoff = 2, wilcoxon.within.subclass = FALSE, 
                    one.against.one = FALSE, mult.test.correction = 0, make.lefse.plots = FALSE, 
                    by_otus = FALSE, levels = rank_names(phy)) { 
     keepvars <- c(class, subclass, subject, "sample")
     keepvars <- unique(keepvars[!is.na(keepvars)])
     samp <- get.samp(phy)[, keepvars]
     if (by_otus) {
          otu <- get.otu.melt(phy, sample_data = FALSE)
          otu.levels <- otu %>% 
               mutate(taxon = otu) %>% 
               group_by(sample,taxon) %>% 
               summarize(pctseqs = sum(pctseqs)) %>% 
               mutate(taxon = gsub(" ", "_", taxon))
     }else{
          otu <- get.otu.melt(phy, sample_data = FALSE)
          otu.list <- lapply(1:length(levels), function(i) {
               lvls <- levels[1:i]
               lvl <- levels[i]
               otu.level <- otu
               otu.level$taxon <- do.call(paste, c(lapply(lvls, function(l) otu[[l]]), sep = "|"))
               otu.level$rank <- lvl
               otu.level2 <- otu.level %>% 
                    group_by(sample, taxon,  rank) %>% 
                    summarize(pctseqs = sum(pctseqs)) %>% 
                    ungroup()
               return(otu.level2)
          })
          otu.levels <- bind_rows(otu.list) %>% 
               mutate(taxon = gsub(" ","_", taxon))
     }
     otu.tbl <- otu.levels %>% 
          dcast(sample ~ taxon, value.var = "pctseqs",fill = 0) %>% 
          left_join(samp, by = "sample") %>% 
          select_(.dots = c(keepvars,lazyeval::interp(~everything())))
     if (is.na(subject) | subject != "sample") {
          otu.tbl <- otu.tbl %>% select(-sample)
     }
     tbl <- otu.tbl %>% t()
     write.table(tbl, "lefse.txt", quote = FALSE, sep = "\t", 
                 col.names = FALSE)
     opt.class <- paste("-c", which(keepvars %in% class))
     opt.subclass <- ifelse(is.na(subclass), "", paste("-s", which(keepvars %in% subclass)))
     opt.subject <- ifelse(is.na(subject), "", paste("-u", which(keepvars %in% subject)))
     format.command <- paste("python format_input.py lefse.txt lefse.in",  opt.class, opt.subclass, opt.subject, "-o 1000000")
     system(format.command)
     lefse.command <- paste("python run_lefse.py lefse.in lefse.res", 
                            "-a", anova.alpha, "-w", wilcoxon.alpha, "-l", lda.cutoff, 
                            "-e", as.numeric(wilcoxon.within.subclass), "-y", as.numeric(one.against.one), 
                            "-s", mult.test.correction)
     system(lefse.command)
     print("Wrote lefse.res")
     lefse.out <- read.table("lefse.res", header = FALSE, sep = "\t") %>% 
          rename(taxon = V1, log.max.pct = V2, direction = V3,lda = V4, p.value = V5)
     if (make.lefse.plots) {
          system("python plot_res.py lefse.res lefse_lda.png")
          print("Wrote lefse_lda.png")
          system("python plot_cladogram.py lefse.res lefse_clado.pdf --format pdf")
          print("Wrote lefse_clado.pdf")
     }
     return(lefse.out)
}

#create phyloseq object containing ADA1 OTU table and mapping file 
#NOTE that OTU table cannot have any spaces or R will act up 
setwd("/home/marcolabuser/Zach/ADA_RS_LP/diversity/lefse/phyloseq")
otu <- read.table(file = "/home/marcolabuser/Zach/ADA_RS_LP/diversity/lefse/phyloseq/feature-table.txt", header = TRUE) 
head(otu)
str(otu)
numeric(otu[2:ncol(otu)])

otu_table(otu, taxa_are_rows = TRUE)
tax <- read.csv(file = "/home/marcolabuser/Zach/ADA_RS_LP/diversity/lefse/phyloseq/taxonomy.csv", header = TRUE)
head(tax)
merged_file <- merge(otu, tax, by.x=c("OTUID"), by.y=c("OTUID"))
head(merged_file)
#load cid94 data
data("cid94")
#move to directory where you downloaded the lefse zip, in my case it was in downloads:
setwd("/home/marcolabuser/Zach/nsegata-lefse-82605a2ae7b7")

#run the new function on the phyloseq object, specifying class and subclass variable names that must be in your sample_data slot of phyloseq object
lefse.tbl <- lefse2(cid94,class="Consistency")


#make a lefse style barplot.. here i filtered for 2 consistencies, big LDA values:
lefse.tbl %>%
 filter(!is.na(lda),
        lda>4,
        direction %in% c("formed stool","liquid")) %>%
 mutate(lda=ifelse(direction=="liquid",-lda,lda)) %>%
 ggplot() +
 geom_bar(aes(x=reorder(taxon,lda),y=lda,fill=direction),color="black",stat="identity") +
 coord_flip()

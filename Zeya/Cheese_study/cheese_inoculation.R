library(ggplot2)
library(reshape2)

path <- "G:/My Drive/UC_Davis/Marco_lab/milk_microbiota/2017_Hilmar_sampling/isolate inoculation/"


tab <- read.csv(file.path(path,"sample information sheet.csv"))

tab.milk <- subset(tab, Species %in% c( "milk_no_slits", "milk_slits", "saline")) 
tab.milk <- subset(tab.milk, Aging_Time %in% c("0D", "11D"))

tab.iso <- subset(tab, Species %in% c("Leuconostoc lactis", "Leuconostoc mesenteroides", 
                                                "Lactobacillus fermentum", "saline"))
tab.iso <- subset(tab.iso, Aging_Time %in% c("Inoculation", "0D", "5D"))

tab.LABr <- subset(tab.iso, Rifampicin %in% c("R", "NC"))
tab.LAB <- subset(tab.iso, Rifampicin %in%  c("WT", "NC"))

# Define function for plot log CFU count over time
CFUoverT <- function(x, path.out, w, h) {
  tabm <- melt(x, id.vars = c("SampleID", "Inoculum", "Species", 
                              "Aging_Time", "Aging_Temperature", "Rifampicin", "pH"))
  tabm$Aging_Time <- factor(tabm$Aging_Time, levels = c("Inoculation", "0D", "5D", "11D"))
  tabm$Species <- factor(tabm$Species, levels = c("saline", "Leuconostoc lactis", 
                                                  "Leuconostoc mesenteroides",
                                                  "Lactobacillus fermentum", 
                                                  "milk_no_slits", "milk_slits"))
  
  pdf(path.out, w, h)
  print(ggplot(tabm, aes(x = Aging_Time, y = value)) +
          geom_boxplot() +
          ylab("Log10(CFU/mL)") +
          facet_grid(variable ~ Species) +
          theme_bw(base_size = 16)) 
  dev.off()
}

CFUoverT(tab.LAB, path.out = file.path(path, "fig/30C_aging_isolates.pdf"), w = 11, h = 7)
CFUoverT(tab.LABr, path.out = file.path(path, "fig/30C_aging_isolatesR.pdf"), w = 11, h = 7)
CFUoverT(tab.milk, path.out = file.path(path, "fig/30C_aging_milk.pdf"), w = 8.65, h = 7)

# KW for CFU count 
TotalCFUKW <- function(x, y, path.out){
  x <- subset(x, Species %in% y)
  fit <- aov(Total_CFU.g ~ Aging_Time, data = x)
  
  attach(x)
  # because Nemenyi is no appropriate for groups with unequal sample sizes
  Z <- posthoc.kruskal.dunn.test(Total_CFU.g, Aging_Time, p.adjust.method = "none")[[3]]
  detach()
  
  write.csv(Z, path.out)
}
LABCFUKW <- function(x, y, path.out){
  x <- subset(x, Species %in% y)
  fit <- aov(LAB_CFU.g ~ Aging_Time, data = x)
  
  attach(x)
  # because Nemenyi is no appropriate for groups with unequal sample sizes
  Z <- posthoc.kruskal.dunn.test(LAB_CFU.g, Aging_Time, p.adjust.method = "none")[[3]]
  detach()
  
  write.csv(Z, path.out)
}
LABRCFUKW <- function(x, y, path.out){
  x <- subset(x, Species %in% y)
  fit <- aov(Rif_LAB_CFU.g ~ Aging_Time, data = x)
  
  attach(x)
  # because Nemenyi is no appropriate for groups with unequal sample sizes
  Z <- posthoc.kruskal.dunn.test(Rif_LAB_CFU.g, Aging_Time, p.adjust.method = "none")[[3]]
  detach()
  
  write.csv(Z, path.out)
}

TotalCFUKW(tab.milk, "saline", file.path(path, "fig/Total_salineKW.csv"))
TotalCFUKW(tab.milk, "milk_no_slits", file.path(path, "fig/Total_milknoslitsKW.csv"))
TotalCFUKW(tab.milk, "milk_slits", file.path(path, "fig/Total_milkslitsKW.csv"))
LABCFUKW(tab.milk, "saline", file.path(path, "fig/LAB_salineKW.csv"))
LABCFUKW(tab.milk, "milk_no_slits", file.path(path, "fig/LAB_milknoslitsKW.csv"))
LABCFUKW(tab.milk, "milk_slits", file.path(path, "fig/LAB_milkslitsKW.csv"))


TotalCFUKW(tab.iso, "saline", file.path(path, "fig/Total_salineKWiso.csv"))
TotalCFUKW(tab.iso, "Leuconostoc lactis", file.path(path, "fig/Total_LLKW.csv"))
TotalCFUKW(tab.iso, "Leuconostoc mesenteroides", file.path(path, "fig/Total_LMKW.csv"))
TotalCFUKW(tab.iso, "Lactobacillus fermentum", file.path(path, "fig/Total_LBFKW.csv"))
LABCFUKW(tab.iso, "saline", file.path(path, "fig/LAB_salineKWiso.csv"))
LABCFUKW(tab.iso, "Leuconostoc lactis", file.path(path, "fig/LAB_LLKW.csv"))
LABCFUKW(tab.iso, "Leuconostoc mesenteroides", file.path(path, "fig/LAB_LMKW.csv"))
LABCFUKW(tab.iso, "Lactobacillus fermentum", file.path(path, "fig/LAB_LBFKW.csv"))
LABRCFUKW(tab.LABr, "Leuconostoc lactis", file.path(path, "fig/LABR_LLKW.csv"))
LABRCFUKW(tab.LABr, "Leuconostoc mesenteroides", file.path(path, "fig/LABR_LMKW.csv"))
LABRCFUKW(tab.LABr, "Lactobacillus fermentum", file.path(path, "fig/LABR_LBFKW.csv"))

# Slit area analysis
tab.area <- read.csv(file.path(path,"image analysis.csv"))
tab.area$Species <- factor(tab.area$Species, levels = c("saline", "Lactobacillus plantarum", 
                                          "milk_no_slits", "milk_slits",
                                          "Leuconostoc lactis ", 
                                          "Leuconostoc mesenteroides", 
                                          "Lactobacillus fermentum "))

# milk inoculated samples 
tab.arem <- subset(tab.area, Species %in% c("saline", "milk_no_slits", "milk_slits"))
tab.arem <- subset(tab.arem, Inoculum != "NC4")
tab.arem <- subset(tab.arem, Inoculum != "NC2")

# Bacterial isolate inoculated samples 
tab.areb <- subset(tab.area, Species %in% c("saline", "Leuconostoc lactis ", 
                                            "Leuconostoc mesenteroides", 
                                            "Lactobacillus fermentum ",
                                            "Lactobacillus plantarum"))
tab.areb <- subset(tab.areb, Inoculum != "NC3")
tab.areb <- subset(tab.areb, Inoculum != "NC1")

# Define function to plot slit area 
SlitArea <- function(x, path.out, w, h){
  pdf(path.out, width = w, height = h)
  print(ggplot(x, aes(x = Species, y = slit_area.)) +
          geom_boxplot(aes(x = Species, y = slit_area.)) +
          ylab("slit area%")+
          xlab("Inoculum")+
          theme_bw(base_size = 16) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

SlitArea(tab.arem, path.out = file.path(path, "fig/slit_area_milk.pdf"), w = 4, h =3.7)
SlitArea(tab.areb, path.out = file.path(path, "fig/slit_area_isolates.pdf"), w = 4, h = 4.63)

## Kruskal Wallis (KW) line figures for each taxa 
library(PMCMR)
SlitAreaKW <- function(x, path.out){
  fit <- aov(slit_area. ~ Species, data = x)
  
  attach(x)
  # because Nemenyi is no appropriate for groups with unequal sample sizes
  Z <- posthoc.kruskal.dunn.test(slit_area., Species, p.adjust.method = "none")[[3]]
  detach()
  
  write.csv(Z, path.out)
}

SlitAreaKW(tab.arem, path.out = file.path(path, "fig/slit_area_milk_KW.csv"))
SlitAreaKW(tab.areb, path.out = file.path(path, "fig/slit_area_isolates_KW.csv"))

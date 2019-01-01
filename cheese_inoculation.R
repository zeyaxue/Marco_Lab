library(ggplot2)
library(reshape2)

path <- "G:/My Drive/UC_Davis/Marco_lab/milk_microbiota/2017_Hilmar_sampling/isolate inoculation/"


tab <- read.csv(file.path(path,"sample information sheet.csv"))

tab.milk <- subset(tab, Inoculumn_species %in% c("Frozen milk", "NC")) 

tab.iso <- subset(tab, Inoculumn_species %in% c("Leuconostoc lactis", "Leuconostoc mesenteroides", 
                                                "Lactobacillus fermentum", "NC"))
tab.LABr <- subset(tab.iso, Rifampicin %in% c("R", "NC"))
tab.LAB <- subset(tab.iso, Rifampicin %in%  c("WT", "NC"))

# Define function for plot log CFU count over time
CFUoverT <- function(x, w, h, path.out) {
  tabm <- melt(x, id.vars = c("SampleID", "Inoculum", "Inoculumn_species", 
                              "Aging_Time", "Aging_Temperature", "Rifampicin"))
  tabm$Aging_Time <- factor(tabm$Aging_Time, levels = c("Inoculation", "0D", "5D"))
  tabm$Inoculumn_species <- factor(tabm$Inoculumn_species, levels = c("NC", 
                                                                      "Leuconostoc lactis", "Leuconostoc mesenteroides", 
                                                                      "Lactobacillus fermentum"))
  
  #pdf(path.out, w, h)
  ggplot(tabm, aes(x = Aging_Time, y = value, color = Inoculumn_species)) +
    geom_boxplot() +
    facet_wrap("variable") +
    theme_bw()
  #dev.off()
}


CFUoverT(tab.LABr)






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

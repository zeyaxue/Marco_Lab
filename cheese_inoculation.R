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


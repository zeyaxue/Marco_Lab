# 0 # Install DAtest 
## https://github.com/Russel88/DAtest

# 1 # Load necessary packages 
library(phyloseq);packageVersion("phyloseq")
library(DAtest)
# Define output/saving path for files
path.out <- "G:/My Drive/UC_Davis/Marco_lab/milk_microbiota/Hilmar_weekly_samples/QIIME2/DAtest"

# 2 # Load input data for DAtest & pre-processing to clean up the data 
## This should be a unrarefied dataset.
## In my case, this phyloseq object that has be removed of 
### (1)contanmination reads using decontam 
### (2) Unassigned phylum using phyloseq 
### (3) taxa < 0.00005 using phyloseq
path.in <- "G:/My Drive/UC_Davis/Marco_lab/milk_microbiota/Hilmar_weekly_samples/QIIME2/decontam/phyloseq/ps05.Phy05.rds"
ps05.phy05 <- readRDS(path.in) 
saveRDS(ps05.phy05, file.path(path.out, "ps05.phy05.rds"))
ps05.phy05
#  add tree 
path.tree <- "G:/My Drive/UC_Davis/Marco_lab/milk_microbiota/Hilmar_weekly_samples/QIIME2/decontam/phyloseq/tree.nwk"
tree <- read_tree(path.tree)
ps05.phy05.tree <- merge_phyloseq(ps05.phy05, phy_tree(tree))
saveRDS(ps05.phy05.tree, file.path(path.out, "ps05.phy05.tree.rds"))
ps05.phy05.tree
# remove negative control and seq control samples
ps <- subset_samples(ps05.phy05.tree, SampleOrCtrl == "Sample")
# remove sample with less than 100 total reads
ps <- prune_samples(sample_sums(ps) >= 100, ps)
saveRDS(ps, file.path(path.out, "ps.rds"))
ps #  USE THIS for DAtest


# 3 # DAtest using CheeseOutcome as predictor
# Use function for individial methods
source(file.path(path.out,"DAtest_ANCOM.R"))  # call DA.anc()
# plot results
plot(mytest)
# Print summary statistics (medians)
summary(mytest)
# Details from the run:
mytest$details
# Average run times of each method:
mytest$run.times
## Power analysis
po.anc <- powerDA(ps, predictor = "CheeseOutcome", covars = "SampleType", R=10,
                  test = "anc", k=c(5,10,15))
plot(po.anc)
summary(po.anc)


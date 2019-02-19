MC\_paper2\_phyloseq\_import
================
Zeya Zhengyao Xue
January 14, 2019

This file follows the analysis in QIIME2, which imports mapping file, feature table, taxonomy and tree as a 4-component phyloseq object.

Setting up packages and working directory
-----------------------------------------

``` r
library(phyloseq);packageVersion("phyloseq")
```

    ## [1] '1.22.3'

``` r
library(vegan);packageVersion("vegan")
```

    ## Warning: package 'vegan' was built under R version 3.4.4

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-1

    ## [1] '2.5.1'

``` r
library(ggplot2);packageVersion("ggplot2")
```

    ## Warning: package 'ggplot2' was built under R version 3.4.4

    ## [1] '3.1.0'

``` r
library(reshape2);packageVersion("reshape2")
```

    ## Warning: package 'reshape2' was built under R version 3.4.3

    ## [1] '1.4.3'

``` r
library(cowplot);packageVersion("cowplot")
```

    ## Warning: package 'cowplot' was built under R version 3.4.3

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    ## [1] '0.9.2'

``` r
library(plyr)
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 3.4.4

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:plyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
path <- "~/Google Drive File Stream/My Drive/UC_Davis/Marco_lab/milk_microbiota/Mock_community/MC_paper2_QIIME2/"
```

Import QIIME2 output and make phyloseq object
---------------------------------------------

2/5/19: add Phase 3 results and import as the same ps object

``` r
# Transpose the QIIME2 output in excel
# For the theoretical values, I can adjust its absolute reads to suit data set by
# changing values in the feature table
SeqTab <- read.table(file.path(path, "feature_table/5file-combined_table-dada2-mc2.txt"), header = TRUE, stringsAsFactors = FALSE)
colnames(SeqTab) <- gsub("X","",colnames(SeqTab))
row.names(SeqTab) <- SeqTab$OTUID
SeqTab <- SeqTab[,-1]
SeqTab <- as.matrix.data.frame(SeqTab)

# Read in metadata file
samdf <- read.csv(file.path(path,"mapping/5file-MC_paper2_samdf.csv"))
rownames(samdf) <- samdf$SampleID

# Read in exported rooted tree 
tree <- read_tree(file.path(path,"5file-tree.nwk"))

# Parse out taxonomic assignment in excel and remove confidence column
# keep only taxa included in the "feature-table-mc2.csv" file
TaxTab <- read.table(file.path(path,"5file-taxonomy-mc2.tsv"), header = TRUE, sep = '\t', na.strings = NA)
rownames(TaxTab) <- TaxTab$FeatureID
TaxTab <- TaxTab[,-1]
TaxTab <- as.matrix.data.frame(TaxTab)

# merge as the phyloseq ready to use
ps <- phyloseq(otu_table(SeqTab, taxa_are_rows = TRUE), tax_table(TaxTab),sample_data(samdf), phy_tree(tree))
ps # 1052 taxa and 212 samples 
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1052 taxa and 212 samples ]
    ## sample_data() Sample Data:       [ 212 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1052 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 1052 tips and 1050 internal nodes ]

``` r
# Filter out singleton
ps <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)
ps # 1049 taxa and 212 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1049 taxa and 212 samples ]
    ## sample_data() Sample Data:       [ 212 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1049 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 1049 tips and 1047 internal nodes ]

``` r
write.csv(ps@otu_table, file.path(path,"feature_table/5file-feature-table-mc2.csv"))
```

Inspect Library Sizes
---------------------

``` r
#add a colum for read counts/library size
samdf$LibrarySize <- sample_sums(ps)
#sort by the library size in increasing order
samdf <- samdf[order(samdf$LibrarySize),]
samdf$Index <- seq(nrow(samdf))
ggplot(samdf, aes(x=Index, y=LibrarySize, color=SampleOrCtrl))+geom_point()
```

![](figs/unnamed-chunk-3-1.png) As expected, some of the UHT\_NC had pretty high read coverage. Also, I decided to not do decontam with the MC data because I know already which taxa are contamination and which taxa are expected.

Filter ASVs based on taxonomy and abundance
-------------------------------------------

``` r
# Show available ranks in the dataset
rank_names(ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

``` r
# Create table, number of features for each phylum
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

    ## 
    ##                        [Thermi]  Actinobacteria   Bacteroidetes 
    ##              36               7             115             193 
    ##     Chloroflexi   Crenarchaeota   Cyanobacteria   Euryarchaeota 
    ##               5               1              11               4 
    ##      Firmicutes    Fusobacteria  Proteobacteria    Spirochaetes 
    ##             434               3             214               8 
    ##   Synergistetes     Tenericutes Verrucomicrobia 
    ##               1              15               2

``` r
# remove asv that was not assigned at phylum level
ps.Phy <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c(""))
ps.Phy # 1013 taxa and 212 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1013 taxa and 212 samples ]
    ## sample_data() Sample Data:       [ 212 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1013 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 1013 tips and 1011 internal nodes ]

``` r
# Explore abundance to decide filter threshold
## Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(ps.Phy),
                MARGIN = ifelse(taxa_are_rows(ps.Phy), yes = 1, no = 2),
                FUN = function(x)(sum(x>0)))
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(ps.Phy),
                     tax_table(ps.Phy))
# plot the abudance per phylum 
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(ps.Phy),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.01,  linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](figs/unnamed-chunk-4-1.png)

``` r
# Most bacteria are in the Actinobacteria, Bacteroidetes, Firmircutes and Proteobacteria
# But I don't actually have Bacteroidetes members...

# Remove very low abudance ASVs independent of sample prevalence
x = taxa_sums(ps.Phy)
keepTaxa = which((x / sum(x)) > 0.00005)
ps.Phy05 <- prune_taxa(names(keepTaxa), ps.Phy)
ps.Phy05 # 197 taxa and 212 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 197 taxa and 212 samples ]
    ## sample_data() Sample Data:       [ 212 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 197 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 197 tips and 196 internal nodes ]

Determine rarefaction level for Alpha and Beta analysis later on
----------------------------------------------------------------

``` r
set.seed(123)
# Define function to calculate alpha diversity and rarefaction levels
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c("SampleID", 'measures'), value.name = 'Alpha_diversity')
    molten_alpha_diversity
  }
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  #convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  rarefaction_curve_data
}
rarefaction_curve_data <- calculate_rarefaction_curves(ps.Phy05, c("Observed", "Shannon","Simpson"), 
                                                       rep(c(1, 2500, 5000, 7500, 10000), each = 10))
rarefaction_curve_data$SampleID <- gsub("X","",rarefaction_curve_data$SampleID)
# summarize alpha diversity
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'SampleID', 'measures'), 
                                        summarise, 
                                        Alpha_diversity_mean = mean(Alpha_diversity), 
                                        Alpha_diversity_sd = sd(Alpha_diversity))
# Add sample data
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, 
                                                data.frame(sample_data(ps.Phy)), 
                                                by.x = 'SampleID', by.y = 'row.names')
```

    ## Warning in merge.data.frame(rarefaction_curve_data_summary,
    ## data.frame(sample_data(ps.Phy)), : column name 'SampleID' is duplicated in
    ## the result

``` r
rarefaction_curve_data_summary_verbose <- rarefaction_curve_data_summary_verbose[, -1]

# plot out rarefaction curve
ggplot(rarefaction_curve_data_summary_verbose, aes(x = Depth, y = Alpha_diversity_mean, 
                                                   ymin = Alpha_diversity_mean - Alpha_diversity_sd, 
                                                   ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                                                   color = SampleID, group = SampleID)) +
  scale_fill_brewer(palette = "Set2") + 
  geom_line()+
  geom_pointrange(size=0.1)+
  theme(legend.position="none") + #  remove legend
  facet_wrap(facets = ~ measures, scales = "free")
```

![](figs/unnamed-chunk-5-1.png)

``` r
# export the library size to a csv file
write.csv(sample_sums(ps.Phy05), file.path(path, "mapping/library_size_ps.Phy05_5file.csv"))

# The diversity leveled out after 2500 seqs per sample, rarefy at 3990
# to keep the sample with the fewest reads that is no a NC (MC.50.2)
ps.Phy.rare <- rarefy_even_depth(ps.Phy05, sample.size = 3990, replace = TRUE, 
                                 rngseed = 123 , trimOTUs = TRUE, verbose = TRUE) 
```

    ## `set.seed(123)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(123); .Random.seed` for the full vector

    ## ...

    ## 9 samples removedbecause they contained fewer reads than `sample.size`.

    ## Up to first five removed samples are:

    ## milliQ.110CB.NC10TB.NC10TC.NC10UB.1  

    ## ...

``` r
ps.Phy.rare # 197 taxa and 203 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 197 taxa and 203 samples ]
    ## sample_data() Sample Data:       [ 203 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 197 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 197 tips and 196 internal nodes ]

Divide the ps object for different dataset
------------------------------------------

### For taxonomy analysis based on "ps.Phy.REC"

``` r
# Clean up the taxonomy level
#make the deepest taxonomic identification in TaxTab as the recommened taxonomy (REC)
RECps <- function(x, ps) {
  TaxTab2 <- as.data.frame(x)
  
  list.g <- as.character(TaxTab2$Genus)
  list.f <- as.character(TaxTab2$Family)
  list.o <- as.character(TaxTab2$Order)
  list.c <- as.character(TaxTab2$Class)
  list.p <- as.character(TaxTab2$Phylum)
  list.k <- as.character(TaxTab2$Kingdom)
  list.REC <- character(length(list.g))
  
  for(i in 1:dim(TaxTab2)[1]){
    G = which(TaxTab2$Genus[i] == "")
    F = which(TaxTab2$Family[i] == "")
    O = which(TaxTab2$Order[i] == "")
    C = which(TaxTab2$Class[i] == "")
    P = which(TaxTab2$Phylum[i] == "")
    K = which(TaxTab2$Kingdom[i] == "")
    if(length(G) == 0){
      list.REC[i] <- list.g[i]
    } else if(length(F) == 0){
      list.REC[i] <- list.f[i]
    } else if(length(O) == 0){
      list.REC[i] <- list.o[i]
    } else if(length(C) == 0){
      list.REC[i] <- list.c[i]
    } else if(length(P) == 0){
      list.REC[i] <- list.p[i]
    } else if(length(K) == 0){
      list.REC[i] <- list.k[i]
    } else
      list.REC[i] <- "meow"
  }
  
  TaxTab2$REC <- list.REC
  TaxTab2$REC <- factor(TaxTab2$REC)
  merge_phyloseq(ps, TaxTab2 %>% as.matrix() %>% tax_table())
}

ps.Phy.glom <- tax_glom(ps.Phy05, taxrank = "Genus", NArm = FALSE)
ps.Phy.REC <- RECps(TaxTab, ps.Phy.glom)
ps.Phy.REC # 100 taxa and 212 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 100 taxa and 212 samples ]
    ## sample_data() Sample Data:       [ 212 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 100 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 100 tips and 99 internal nodes ]

``` r
# Define function to get the average and raw value of the feature table 
AveFeatTab <- function(ps){
  df <- psmelt(ps)
  df.cast <- dcast(df, REC ~ For_AVE, mean, value.var = "Abundance")
  df.cast[, -1] <- lapply( df.cast[ , -1], function(x) x/sum(x, na.rm=TRUE) )
  df.cast
}
FeatTab <- function(ps){
  df <- psmelt(ps)
  df.cast <- dcast(df, For_AVE + SampleID ~ REC, value.var = "Abundance")
  df.cast
}


# 1. DNA extraction lysis method using the Total kit; BCMC1
samdf_BCMC1_total_lyse <- read.csv(file.path(path, "mapping/MC_paper2_samdf_BCMC1_total_lyse.csv"))
rownames(samdf_BCMC1_total_lyse) <- samdf_BCMC1_total_lyse$SampleID
ps1 <- phyloseq(sample_data(samdf_BCMC1_total_lyse), otu_table(ps.Phy.REC),
                tax_table(ps.Phy.REC), phy_tree(ps.Phy.REC))
ps1 # 91 taxa and 23 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 100 taxa and 23 samples ]
    ## sample_data() Sample Data:       [ 23 samples by 13 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 100 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 100 tips and 99 internal nodes ]

``` r
write.csv(AveFeatTab(ps1), file.path(path, "feature_table/FeatTab_BCMC1_total_lyseAVE.csv"))


# 2. DNA extraction kit with MoBio
samdf_BCMC1_MoBio <- read.csv(file.path(path, "mapping/MC_paper2_samdf_BCMC1_MoBio.csv"))
rownames(samdf_BCMC1_MoBio) <- samdf_BCMC1_MoBio$SampleID
ps2 <- phyloseq(sample_data(samdf_BCMC1_MoBio), otu_table(ps.Phy.REC),
                tax_table(ps.Phy.REC), phy_tree(ps.Phy.REC))
ps2 # 91 taxa and 6 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 100 taxa and 6 samples ]
    ## sample_data() Sample Data:       [ 6 samples by 13 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 100 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 100 tips and 99 internal nodes ]

``` r
write.csv(AveFeatTab(ps2), file.path(path, "feature_table/FeatTab_BCMC1_MoBioAVE.csv"))

# 3. DNA extraction kit choice on 10 ml UHT milk; BCMC2
samdf_BCMC2_kit_10ml <- read.csv(file.path(path, "mapping/MC_paper2_samdf_BCMC2_kit_10ml.csv"))
rownames(samdf_BCMC2_kit_10ml) <- samdf_BCMC2_kit_10ml$SampleID
ps3 <- phyloseq(sample_data(samdf_BCMC2_kit_10ml), otu_table(ps.Phy.REC),
                tax_table(ps.Phy.REC), phy_tree(ps.Phy.REC))
ps3 # 91 taxa and 30 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 100 taxa and 30 samples ]
    ## sample_data() Sample Data:       [ 30 samples by 13 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 100 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 100 tips and 99 internal nodes ]

``` r
write.csv(AveFeatTab(ps3), file.path(path, "feature_table/FeatTab_BCMC2_kit_10mlAVE.csv"))

# 4. PMA treatment
samdf_PMA <- read.csv(file.path(path, "mapping/MC_paper2_samdf_PMA.csv"))
rownames(samdf_PMA) <- samdf_PMA$SampleID
ps4 <- phyloseq(sample_data(samdf_PMA), otu_table(ps.Phy.REC),
                tax_table(ps.Phy.REC), phy_tree(ps.Phy.REC))
ps4 # 91 taxa and 15 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 100 taxa and 15 samples ]
    ## sample_data() Sample Data:       [ 15 samples by 13 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 100 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 100 tips and 99 internal nodes ]

``` r
write.csv(AveFeatTab(ps4), file.path(path, "feature_table/FeatTab_PMA.csv"))

# 5. Sample storage
samdf_sample_storage <- read.csv(file.path(path, "mapping/MC_paper2_samdf_sample_storage.csv"))
rownames(samdf_sample_storage) <- samdf_sample_storage$SampleID
ps5 <- phyloseq(sample_data(samdf_sample_storage), otu_table(ps.Phy.REC),
                tax_table(ps.Phy.REC), phy_tree(ps.Phy.REC))
ps5 # 91 taxa and 24 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 100 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 100 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 100 tips and 99 internal nodes ]

``` r
write.csv(AveFeatTab(ps5), file.path(path, "feature_table/FeatTab_sample_storage.csv"))

# 6. Sample Volume
samdf_total_vol <- read.csv(file.path(path, "mapping/MC_paper2_samdf_total_vol.csv"))
rownames(samdf_total_vol) <- samdf_total_vol$SampleID
ps6 <- phyloseq(sample_data(samdf_total_vol), otu_table(ps.Phy.REC),
                tax_table(ps.Phy.REC), phy_tree(ps.Phy.REC))
ps6 # 91 taxa and 18 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 100 taxa and 18 samples ]
    ## sample_data() Sample Data:       [ 18 samples by 12 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 100 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 100 tips and 99 internal nodes ]

``` r
write.csv(AveFeatTab(ps6), file.path(path, "feature_table/FeatTab_total_vol.csv"))

# 7. phase 3 kit, cleanup and lysing method 
samdf_phase3 <- read.csv(file.path(path, "mapping/phase3_samdf.csv"))
rownames(samdf_phase3) <- samdf_phase3$SampleID
ps7 <- phyloseq(sample_data(samdf_phase3), otu_table(ps.Phy.REC),
                tax_table(ps.Phy.REC), phy_tree(ps.Phy.REC))
ps7 # 100 taxa and 51 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 100 taxa and 51 samples ]
    ## sample_data() Sample Data:       [ 51 samples by 15 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 100 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 100 tips and 99 internal nodes ]

``` r
write.csv(AveFeatTab(ps7), file.path(path, ("feature_table/FeatTab_phase3.csv")))


# create ps list with all phyloseq objects
psList <- list(ps1, ps2, ps3, ps4, ps5, ps6, ps7)
names(psList) <- c(1:7)

# use the function to write feature table (not average) list and write out csv
FeatTab.list <- lapply(psList, FeatTab)
mapply(write.csv, FeatTab.list, file = paste0(path, "feature_table/FeatTab_", names(FeatTab.list), ".csv"))
```

    ## $`1`
    ## NULL
    ## 
    ## $`2`
    ## NULL
    ## 
    ## $`3`
    ## NULL
    ## 
    ## $`4`
    ## NULL
    ## 
    ## $`5`
    ## NULL
    ## 
    ## $`6`
    ## NULL
    ## 
    ## $`7`
    ## NULL

#### Taxonomy plot

### For alpha and beta div based on ps.Phy.rare

``` r
ps1.rare <- rarefy_even_depth(ps1, sample.size = 3990, replace = TRUE, rngseed = 123 , trimOTUs = TRUE, verbose = TRUE) 
```

    ## `set.seed(123)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(123); .Random.seed` for the full vector

    ## ...

    ## 75OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
ps2.rare <- rarefy_even_depth(ps2, sample.size = 3990, replace = TRUE, rngseed = 123 , trimOTUs = TRUE, verbose = TRUE) 
```

    ## `set.seed(123)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(123); .Random.seed` for the full vector

    ## ...

    ## 90OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
ps3.rare <- rarefy_even_depth(ps3, sample.size = 3990, replace = TRUE, rngseed = 123 , trimOTUs = TRUE, verbose = TRUE) 
```

    ## `set.seed(123)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(123); .Random.seed` for the full vector

    ## ...

    ## 1 samples removedbecause they contained fewer reads than `sample.size`.

    ## Up to first five removed samples are:

    ## 10UB.1   

    ## ...

    ## 50OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
ps4.rare <- rarefy_even_depth(ps4, sample.size = 3990, replace = TRUE, rngseed = 123 , trimOTUs = TRUE, verbose = TRUE) 
```

    ## `set.seed(123)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(123); .Random.seed` for the full vector

    ## ...

    ## 46OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
ps5.rare <- rarefy_even_depth(ps5, sample.size = 3990, replace = TRUE, rngseed = 123 , trimOTUs = TRUE, verbose = TRUE) 
```

    ## `set.seed(123)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(123); .Random.seed` for the full vector

    ## ...

    ## 83OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
ps6.rare <- rarefy_even_depth(ps6, sample.size = 3990, replace = TRUE, rngseed = 123 , trimOTUs = TRUE, verbose = TRUE) 
```

    ## `set.seed(123)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(123); .Random.seed` for the full vector

    ## ...

    ## 52OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
ps7.rare <- rarefy_even_depth(ps7, sample.size = 3990, replace = TRUE, rngseed = 123 , trimOTUs = TRUE, verbose = TRUE) 
```

    ## `set.seed(123)` was used to initialize repeatable random subsampling.

    ## Please record this for your records so others can reproduce.

    ## Try `set.seed(123); .Random.seed` for the full vector

    ## ...

    ## 61OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
ps1.rare #  27 taxa and 23 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 25 taxa and 23 samples ]
    ## sample_data() Sample Data:       [ 23 samples by 13 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 25 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 25 tips and 24 internal nodes ]

``` r
ps2.rare # 10 taxa and 6 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 10 taxa and 6 samples ]
    ## sample_data() Sample Data:       [ 6 samples by 13 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 10 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 10 tips and 9 internal nodes ]

``` r
ps3.rare # 55 taxa and 29 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 50 taxa and 29 samples ]
    ## sample_data() Sample Data:       [ 29 samples by 13 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 50 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 50 tips and 49 internal nodes ]

``` r
ps4.rare # 61 taxa and 15 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 54 taxa and 15 samples ]
    ## sample_data() Sample Data:       [ 15 samples by 13 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 54 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 54 tips and 53 internal nodes ]

``` r
ps5.rare # 18 taxa and 24 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 17 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 17 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 17 tips and 16 internal nodes ]

``` r
ps6.rare # 51 taxa and 18 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 48 taxa and 18 samples ]
    ## sample_data() Sample Data:       [ 18 samples by 12 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 48 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 48 tips and 47 internal nodes ]

``` r
ps7.rare # 39 taxa and 51 samples
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 39 taxa and 51 samples ]
    ## sample_data() Sample Data:       [ 51 samples by 15 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 39 taxa by 8 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 39 tips and 38 internal nodes ]

#### Beta diversity

``` r
# Define function to get distance matrix 
# I have to use distance() from vegan packages due to in compatibility 
# https://github.com/joey711/phyloseq/issues/918
DistTab <- function(ps, DistMethod){
  df <- psmelt(ps)
  df.cast <- dcast(df, For_AVE + SampleID ~ REC, value.var = "Abundance")
  df.cast <- df.cast[,-1]
  row.names(df.cast) <- df.cast[,1]
  df.cast <- df.cast[,-1]
  dist <- vegdist(df.cast, method = DistMethod) %>% as.matrix()
}

# create ps list with all phyloseq objects
psList.rare <- list(ps1.rare, ps2.rare, ps3.rare, ps4.rare, ps5.rare, ps6.rare, ps7.rare)
names(psList.rare) <- c(1:7)

# use the function to write bray distance list and write out as csv files
BrayList <- lapply(psList.rare, DistTab, DistMethod = "bray")
mapply(write.csv, BrayList, file = paste0(path, "distance/bray_", names(psList.rare), ".csv"))
```

    ## $`1`
    ## NULL
    ## 
    ## $`2`
    ## NULL
    ## 
    ## $`3`
    ## NULL
    ## 
    ## $`4`
    ## NULL
    ## 
    ## $`5`
    ## NULL
    ## 
    ## $`6`
    ## NULL
    ## 
    ## $`7`
    ## NULL

``` r
# write out weighted list as well
WeightList <- lapply(psList.rare, UniFrac, weighted = TRUE)
WeightList <- lapply(WeightList, as.matrix)
mapply(write.csv, WeightList, file = paste0(path, "distance/weighted_", names(psList.rare), ".csv"))
```

    ## $`1`
    ## NULL
    ## 
    ## $`2`
    ## NULL
    ## 
    ## $`3`
    ## NULL
    ## 
    ## $`4`
    ## NULL
    ## 
    ## $`5`
    ## NULL
    ## 
    ## $`6`
    ## NULL
    ## 
    ## $`7`
    ## NULL

``` r
# Import organized and subsetted distance csv file for making bar plots
```

BBBBBB USE CCA TRIPLOT TO SEE IF THERE'S ASSOCIATION BETWEEN METHOD AND TAXA

AAAAAA Alpha could be more useful when using actual milk samples?

Percent recovery of GBS DNA with meconium extractions generated from various extraction kits. Data are mean ± S
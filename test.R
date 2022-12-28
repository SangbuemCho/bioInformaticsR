
# Package required
## BiocManager package
### install.packages("BiocManager") 
## data2 package
## BiocManager::install("dada2", version = "3.16")

library(dada2)

# load data
path = "/Users/sangbuemcho/Documents/bioinformatics_R/Miyoung" # where the sequence files located
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE)) # get all forward sequence files, pattern should be altered
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE)) # get all reverse sequence files, pattern should be altered
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) # create list of file name

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # get all filtered forward sequence files, pattern should be altered
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) # get all filtered reverse sequence files, pattern should be altered
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# Filter and Trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn the Error rates

load("tempSave.rda") # data from line 24 ~ 32
# save(taxa, errF, errR, mergers, file = "testSave.rda")
errF <- learnErrors(filtFs, multithread=TRUE) # takes huge time and ram
errR <- learnErrors(filtRs, multithread=TRUE) # takes huge time and ram
plotErrors(errF, nominalQ=TRUE)


# Sample inference

dadaFs <- dada(filtFs, err=errF, multithread=TRUE) 
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Merges paired 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])


# Construct sequence Table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))


# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab)
rownames(seqtab.nochim)


# Assign Taxonomy
# taxa file "silva_nr99_v138.1_train_set.fa.gz"  source (below URL)
## https://zenodo.org/record/4587955#.Y6v9fezP1hE 

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE) # takes long times
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


# Alpha diversity (phyloseq)


# BiocManager::install("phyloseq") # package install should be advanced
library(phyloseq)
# BiocManager::install("Biostrings") # package install should be advanced
library(Biostrings)
# install.packages("ggplot2")
library(ggplot2)
theme_set(theme_bw())


## construct data frame file
samples.out <- rownames(seqtab.nochim) # seqtab.nochim was defined at line 54
samples.out

group = rep(c("CON", "TRT"), each = 8) # In this analysis, there are two groups and 8 times
day = rep(seq(2, 9, 1), times = 2)

samdf <- data.frame(Group=group, Day=day)
rownames(samdf) <- samples.out

# construct a phyloseq object directly from the dada2 outputs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps

# BiocManager::install("DECIPHER") # package install should be advanced

library(DECIPHER); packageVersion("DECIPHER") 

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Visualize alpha-diversity:
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="Group")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu)) # probability evaluation
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")


# Barplot:
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))
top20
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~Group)


# relative abundance file creation
dat = ps.top20@otu_table
# remove NA
na.omit(dat) %>% as.data.frame() -> dat

# Conduct PCoA
euclidean_dist <- vegan::vegdist(dat, method = "euclidean")
euclidean_pcoa <- ecodist::pco(euclidean_dist)
euclidean_pcoa

euclidean_pcoa_df <- data.frame(pcoa1 = euclidean_pcoa$vectors[,1], 
                                pcoa2 = euclidean_pcoa$vectors[,2])


euclidean_plot <- ggplot(data = euclidean_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
  geom_point() +
  labs(x = "PC1",
       y = "PC2",
       title = "Euclidean PCoA with CLR transformation") +
  theme(title = element_text(size = 12)) # makes titles smaller

euclidean_plot

# PERMANOVA
dat$group = row.names(dat)
dat$group <- as.factor(dat$group)
rownames(dat) <- NULL

PERMANOVA(dat[, -1], dat$group)
X = dat[, -ncol(dat)]
D = DistContinuous(X)
X
D
PERMANOVA(D, dat$group)





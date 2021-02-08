
#install DESeq2 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("tximportData")
install.packages("readr")
BiocManager::install("tximport")
BiocManager::install("rhdf5")

library(tximportData)
library( "DESeq2" )
library(readr)
library("tximport")
library("rhdf5")
#browseVignettes("tximport")

# kallisto abundance.tsv
samples <- read.table("samples.txt", header=TRUE)
samples$condition <- factor(rep(c("untreated","treated"),each=3))
rownames(samples) <- samples$run

# Next we specify the path to the files using the appropriate columns of samples, and we read in a table that links transcripts to genes for this dataset.

# kallisto with TSV files

files <- list.files("kallisto/",full.names = T)
files <- paste0(files, "/","abundance.tsv.gz")
names(files) <- samples$run
tx2gene <- read_csv("tx2gene.gencode.v27.csv")
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, txOut = FALSE)


#kallisto abundance.h5

files <- list.files("kallisto/",full.names = T)
files <- paste0(files, "/","abundance.h5")
names(files) <- samples$run
tx2gene <- read_csv("tx2gene.gencode.v27.csv")
txi.kallisto.h5 <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE,  txOut = FALSE)

# txOut = FALSE, Summarizes abundances, counts, lengths,  from transcript- to gene-level.


## construct a DESeqDataSet from the txi.kallisto.h5 object and sample information in samples.

dds <- DESeqDataSetFromTximport(txi.kallisto.h5,
                                   colData = samples,
                                   design = ~ condition)

## Pre-filtering 
##While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we increase the speed of the transformation and testing functions within DESeq2. Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total.


keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Note on factor levels
# By default, R will choose a reference level for factors based on alphabetical order. Then, if you never tell the DESeq2 functions which level you want to compare against (e.g. which level represents the control group), the comparisons will be based on the alphabetical order of the levels. There are two solutions: you can either explicitly tell results which comparison to make using the contrast argument (this will be shown later), or you can explicitly set the factors levels. In order to see the change of reference levels reflected in the results names, you need to either run DESeq or nbinomWaldTest/nbinomLRT after the re-leveling operation. Setting the factor levels can be done in two ways, either using factor:
#dds$condition <- factor(dds$condition, levels = c("untreated","treated"))

dds$condition <- relevel(dds$condition, ref = "untreated")



# Collapsing technical replicates

# ddseg <- makeExampleDESeqDataSet(m=12)
# ddseg$condition
# # make data with two technical replicates for three samples
# ddseg$sample <- factor(sample(paste0("sample",rep(1:9, c(2,1,1,2,1,1,2,1,1)))))
# ddseg$run <- paste0("run",1:12)
# 
# ddsColl <- collapseReplicates(ddseg, ddseg$sample, ddseg$run)
# ddsColl$runsCollapsed


#Differential expression analysis

dds <- DESeq(dds)
res <- results(dds, alpha=0.1)
res
#The text, condition treated vs untreated, tells you that the estimates are of the logarithmic fold change log2(treated/untreated).


# Log fold change shrinkage for visualization and ranking

# Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. To shrink the LFC, we pass the dds object to the function lfcShrink. Below we specify to use the apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.
# https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html

resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="normal")


# We can summarize some basic tallies using the summary function.


summary(res)

# How many adjusted p-values were less than 0.1?
  
sum(res$padj < 0.1, na.rm=TRUE)

# change the alpha
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res$padj < 0.05, na.rm=TRUE)

plotMA(res, ylim=c(-2,2))

# It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
plotMA(resLFC, ylim=c(-2,2))

# After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices:
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

# Data transformations

# rlog

# he running times are shorter when using blind=FALSE and if the function DESeq has already been run, because then it is not necessary to re-estimate the dispersion values. 
# The assay function is used to extract the matrix of normalized values.

rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)


#  log2 transformation + 1
ntd <- normTransform(dds,f = rlog)
head(assay(ntd), 3)

# write result to csv

write.csv(as.data.frame(res), 
          file="condition_treated_results.csv")


### heat map


library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
row.names(df) <- colnames(heat.df)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# PCA
plotPCA(ntd, intgroup=c("condition"))

plotPCA(rld, intgroup=c("condition"))


# Heatmap of the sample-to-sample distances

# Another use of the transformed data is sample clustering. Here, we apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.

sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition)
colors <- colorRampPalette( rev(brewer.pal(6, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples. We have to provide a hierarchical clustering hc to the heatmap function based on the sample distances, or else the heatmap function would calculate a clustering based on the distances between the rows/columns of the distance matrix.
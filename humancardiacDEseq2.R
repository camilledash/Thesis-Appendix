#Goal: use rat count data from Naqvi (2019) to calculate sex bias in gene expression using deSeq2

#install packages/libraries
##dplyr
library(dplyr)
##tidyverse
#install.packages("tidyverse")
library(tidyverse)
##GEOquery
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("GEOquery")
library(GEOquery)
##DESeq2
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
library(ggpubr)
packageVersion('ggpubr')

#STEP 1: DATA READING AND WRANGLING
#read in counts data
human_counts_raw <- as.matrix(read.table(file = "gtex.filt.salmon.tximport.unadj.counts.txt"))
human_counts_raw_df <- as.data.frame(human_counts_raw)

#read in metadata
metadata <- as.data.frame(read.table(file = "human.metadata.txt"))
##subset metadata: only relevant columns, only heart
metadata_heart <- metadata[metadata$V4 == "Heart", ]
metadata_heart <- select(metadata_heart, c(1, 3, 4))

#edit/pare down data (only heart columns)
colnames(human_counts_raw_df)
human_counts_heart <- human_counts_raw_df %>% select(c(metadata_heart$V1))

##rename rows in metadata_heart so they reflect column names (sample names) in human_counts_heart data
colnames(human_counts_heart)
rownames(metadata_heart) <- c(metadata_heart$V1)

#aligning counts and metadata
all(colnames(human_counts_heart) %in% rownames(metadata_heart)) ##TRUE
#check that order matches
all(colnames(human_counts_heart) == rownames(metadata_heart)) ##TRUE

#STEP 2: CREATE DESeqDataSet OBJECT
dds_human <- DESeqDataSetFromMatrix(countData = round(human_counts_heart), 
                                  colData = metadata_heart, 
                                  design = ~V3) 
dds_human
##added "round" into above code b/c got error indicating that some values in assay are not integers

#prefiltering: remove rows w/ low gene counts (less than 10)
keep <- rowSums(counts(dds_human)) >=10
dds_human <- dds_human[keep,] ##removed ~2000 low gene counts (from original ~20000)

#set factor level (using male as reference level)
dds_human$V3 <- relevel(dds_human$V3, ref = "Male")
dds_human$V3 ##Male is listed first in levels, indicating it is reference level


#STEP 3: RUN DEseq
deseq_dds_human <- DESeq(dds_human)
results_dds_human <- results(deseq_dds_human, alpha = 0.05)
results_dds_human

#STEP 4: Explore results
##general
summary(results_dds_human)
resultsNames(deseq_dds_human)

##plotting
plotMA(results_dds_human)

#STEP 5: Export results
##reorder so lowest pvalues first
results_dds_human_ordered <- results_dds_human[order(results_dds_human$padj),]
results_dds_human_ordered
##make csv
write.csv(as.data.frame(results_dds_human_ordered),
          file = "humancardiacDEseq2_results.csv")

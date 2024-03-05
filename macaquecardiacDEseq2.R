#Goal: use macaque count data from Naqvi (2019) to calculate sex bias in gene expression using deSeq2

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

#STEP 1: DATA READING AND WRANGLING
#read in counts data
macaque_counts_raw <- as.matrix(read.table(file = "GSE125483_cyno.salmon.tximport.counts.txt"))
macaque_counts_raw_df <- as.data.frame(macaque_counts_raw)
##edit/pare down data (only heart columns)
colnames(macaque_counts_raw_df)
macaque_counts_raw_df %>% select(contains("Heart"))
macaque_counts_heart <- select(macaque_counts_raw_df, c(10, 23, 27, 47, 52, 60))

#read in metadata
gse <- getGEO(GEO = "GSE125483", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[3]]))
head(metadata)
##subset metadata: only relevant columns, only heart
metadata_heart <- select(metadata, c(1,8))
metadata_heart <- metadata_heart[metadata_heart$source_name_ch1 == "Heart", ]

##create column for sex (since not in metadata): order of macaques in metadata_heart is 1496, 1438, 1083, 1148, 1394, 1255, corresponding with M, F, M, F, M, F (according to macaque_counts_heart)
sex_column <- c("M", "F", "M", "F", "M", "F")
##add sex column to metadata_heart
metadata_heart$sex <- sex_column

##rename rows in metadata_heart so they reflect column names (sample names) in macaque_counts_heart data
colnames(macaque_counts_heart)
metadata_heart$title
rownames(metadata_heart) <- c("Cyno_Male_Heart_L16_1496", "Cyno_Female_Heart_L16_1438", "Cyno_Male_Heart_L16_1083", "Cyno_Female_Heart_L16_1148", "Cyno_Male_Heart_L16_1394", "Cyno_Female_Heart_L16_1255")

#aligning counts and metadata
all(colnames(macaque_counts_heart) %in% rownames(metadata_heart)) ##TRUE
##reorder columns in macaque_counts_heart so its same order as rows in metadata_heart
rownames(metadata_heart)
colnames(macaque_counts_heart)

macaque_counts_heart <- macaque_counts_heart[, c(5, 6, 1, 2, 4, 3)]

#check that order matches
all(colnames(macaque_counts_heart) == rownames(metadata_heart)) ##TRUE

#STEP 2: CREATE DESeqDataSet OBJECT
dds_macaque <- DESeqDataSetFromMatrix(countData = round(macaque_counts_heart), 
                                  colData = metadata_heart, 
                                  design = ~sex) 
dds_macaque
##added "round" into above code b/c got error indicating that some values in assay are not integers
#prefiltering: remove rows w/ low gene counts (less than 10)
keep <- rowSums(counts(dds_macaque)) >=10
dds_macaque <- dds_macaque[keep,] ##removed ~20000 low gene counts (from original ~45500)

#set factor level (using male as reference level)
dds_macaque$sex <- relevel(dds_macaque$sex, ref = "M")
dds_macaque$sex ##M is listed first in levels, indicating it is reference level

#STEP 3: RUN DEseq
deseq_dds_macaque <- DESeq(dds_macaque)
results_dds_macaque <- results(deseq_dds_macaque, alpha = 0.05)
results_dds_macaque

#STEP 4: Explore results
##general
summary(results_dds_macaque)
resultsNames(deseq_dds_macaque)

##plotting
plotMA(results_dds_macaque)

#STEP 5: EXPORT DATA
##reoder so lowest padj values are first
results_dds_macaque_reodered <- results_dds_macaque[order(results_dds_macaque$padj),]
write.csv(as.data.frame(results_dds_macaque_reodered), 
          file = "macaquecardiacDEseq_results.csv")

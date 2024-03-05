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


#STEP 1: DATA READING AND WRANGLING
#read in counts data
rat_counts_raw <- as.matrix(read.table(file = "GSE125483_rat.salmon.tximport.counts.txt"))
rat_counts_raw_df <- as.data.frame(rat_counts_raw)
##edit/pare down data (only heart columns)
colnames(rat_counts_raw_df)
rat_counts_raw_df %>% select(contains("Heart"))
rat_counts_heart <- select(rat_counts_raw_df, c(9, 26, 34, 50, 57, 65))

#read in metadata
gse <- getGEO(GEO = "GSE125483", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[2]]))
head(metadata)
##subset metadata: only relevant columns, only heart
metadata_heart <- select(metadata, c(1,8))
metadata_heart <- metadata_heart[metadata_heart$source_name_ch1 == "Heart", ]

##create column for sex (since not in metadata): order of rats in metadata_heart is 1421, 1505, 1408, 1080, 1281, 1180 corresponding with M, M, F, M, F, F (according to rat_counts_heart)
colnames(rat_counts_heart)
sex_column <- c("M", "M", "F", "M", "F", "F")
##add sex column to metadata_heart
metadata_heart$sex <- sex_column

##rename rows in metadata_heart so they reflect column names (sample names) in macaque_counts_heart data
colnames(rat_counts_heart)
metadata_heart$title
rownames(metadata_heart) <- c("Rat_Male_Heart_L16_1421", "Rat_Male_Heart_L16_1505", "Rat_Female_Heart_L16_1408", "Rat_Male_Heart_L16_1080", "Rat_Female_Heart_L16_1281", "Rat_Female_Heart_L16_1180")

#aligning counts and metadata
all(colnames(rat_counts_heart) %in% rownames(metadata_heart)) ##TRUE
##reorder columns in rat_counts_heart so its same order as rows in metadata_heart
rownames(metadata_heart)
colnames(rat_counts_heart)

rat_counts_heart <- rat_counts_heart[, c(5, 6, 4, 1, 3, 2)]

#check that order matches
all(colnames(rat_counts_heart) == rownames(metadata_heart)) ##TRUE


#STEP 2: CREATE DESeqDataSet OBJECT
dds_rat <- DESeqDataSetFromMatrix(countData = round(rat_counts_heart), 
                                      colData = metadata_heart, 
                                      design = ~sex) 
dds_rat
##added "round" into above code b/c got error indicating that some values in assay are not integers
#prefiltering: remove rows w/ low gene counts (less than 10)
keep <- rowSums(counts(dds_rat)) >=10
dds_rat <- dds_rat[keep,] ##removed ~15000 low gene counts (from original ~41000)

#set factor level (using male as reference level)
dds_rat$sex <- relevel(dds_rat$sex, ref = "M")
dds_rat$sex ##M is listed first in levels, indicating it is reference level

#STEP 3: RUN DEseq
deseq_dds_rat <- DESeq(dds_rat)
results_dds_rat <- results(deseq_dds_rat, alpha = 0.05)
results_dds_rat

#STEP 4: Explore results
##general
summary(results_dds_rat)
resultsNames(deseq_dds_rat)

##plotting
plotMA(results_dds_rat)

#STEP 5: EXPORT DATA
##reoder so lowest padj values are first
results_dds_rat_reodered <- results_dds_rat[order(results_dds_rat$padj),]
write.csv(as.data.frame(results_dds_rat_reodered), 
          file = "ratcardiacDEseq_results.csv")

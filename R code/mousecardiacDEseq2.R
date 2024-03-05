#Goal: use mouse count data from Naqvi (2019) to calculate sex bias in gene expression using deSeq2

#install packages/libraries
##dplyr
library(dplyr)
##tidyverse
#install.packages("tidyverse")
library(tidyverse)
##GEOquery
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("GEOquery")
library(GEOquery)
##DESeq2
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("DESeq2")
library(DESeq2)

#STEP 1: DATA READING AND WRANGLING
#read in counts data
mouse_counts_raw <- as.matrix(read.table(file = "GSE125483_mouse.salmon.tximport.counts.txt"))
mouse_counts_raw_df <- as.data.frame(mouse_counts_raw)
##edit/pare down data (only heart columns)
colnames(mouse_counts_raw_df)
mouse_counts_heart <- select(mouse_counts_raw_df, c(1, 9, 10, 48, 66, 76))

#read in metadata
gse <- getGEO(GEO = "GSE125483", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)
##subset metadata: only relevant columns, only heart
metadata_heart <- select(metadata, c(1,8))
metadata_heart <- metadata_heart[metadata_heart$source_name_ch1 == "Heart", ]
##create column for sex (since not in metadata): order of mice in metadata_heart is 1279, 1272, 1600, 1428, 1230, 1379, corresponding with F, M, F, M, F, M (according to mouse_counts_heart)
sex_column <- c("F", "M", "F", "M", "F", "M")
##add sex column to metadata_heart
metadata_heart$sex <- sex_column
##rename rows in metadata_heart so they reflect column names (sample names) in mouse_counts_heart data
colnames(mouse_counts_heart)
metadata_heart$title

rownames(metadata_heart) <- c("Mouse_Female_Heart_L16_1279", "Mouse_Male_Heart_L16_1272", "Mouse_Female_Heart_L16_1600", "Mouse_Male_Heart_L16_1428", "Mouse_Female_Heart_L16_1230", "Mouse_Male_Heart_L16_1379")

#aligning counts and metadata
all(colnames(mouse_counts_heart) %in% rownames(metadata_heart)) ##TRUE
##reorder columns in mouse_counts_heart so its same order as rows in metadata_heart
rownames(metadata_heart)
colnames(mouse_counts_heart)

mouse_counts_heart_reordered <- mouse_counts_heart[, c(3, 2, 5, 4, 1, 6)]

#check that order matches
all(colnames(mouse_counts_heart_reordered) == rownames(metadata_heart)) ##TRUE

#STEP 2: CREATE DESeqDataSet OBJECT
dds_mouse <- DESeqDataSetFromMatrix(countData = round(mouse_counts_heart_reordered), 
                        colData = metadata_heart, 
                       design = ~sex) 
dds_mouse
##added "round" into above code b/c got error indicating that some values in assay are not integers
#prefiltering: remove rows w/ low gene counts (less than 10)
keep <- rowSums(counts(dds_mouse)) >=10
dds_mouse <- dds_mouse[keep,] ##removed ~25000 low gene counts (~halved it)

#set factor level (using male as reference level)
dds_mouse$sex <- relevel(dds_mouse$sex, ref = "M")
dds_mouse$sex ##M is listed first in levels, indicating it is reference level

#STEP 3: RUN DEseq
deseq_dds_mouse <- DESeq(dds_mouse)
results_dds_mouse <- results(deseq_dds_mouse, alpha = 0.05)
results_dds_mouse

#STEP 4: Explore results
##general
summary(results_dds_mouse)
resultsNames(deseq_dds_mouse)

##plotting
plotMA(results_dds_mouse)

#STEP 5: EXPORT DATA
##reoder so lowest padj values are first
results_dds_mouse_reodered <- results_dds_mouse[order(results_dds_mouse$padj),]
write.csv(as.data.frame(results_dds_mouse_reodered), 
          file = "mousecardiacDEseq_results.csv")

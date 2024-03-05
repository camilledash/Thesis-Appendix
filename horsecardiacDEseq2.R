#Goal: use horse count data from Peng (2023) to calculate sex bias in gene expression using deSeq2

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
#tximport to read the .quant.sf files
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("tximport")
library(tximport)

#STEP 1: DATA READING AND WRANGLING
#read in counts data
##AH1
file_pathAH1 <- "AH1_Heart.quant.sf"
sample_AH1 <- data.frame(sampleName = "AH1_Heart_F", fileName = file_pathAH1)
txi_AH1 <- tximport(files = sample_AH1$fileName, type = "salmon", txOut = TRUE)
gene_counts_AH1 <- txi_AH1$counts
##AH2
file_pathAH2 <- "AH2_Heart.quant.sf"
sample_AH2 <- data.frame(sampleName = "AH2_Heart_F", fileName = file_pathAH2)
txi_AH2 <- tximport(files = sample_AH2$fileName, type = "salmon", txOut = TRUE)
gene_counts_AH2 <- txi_AH2$counts
##AH3
file_pathAH3 <- "AH3_Heart.quant.sf"
sample_AH3 <- data.frame(sampleName = "AH3_Heart_M", fileName = file_pathAH3)
txi_AH3 <- tximport(files = sample_AH3$fileName, type = "salmon", txOut = TRUE)
gene_counts_AH3 <- txi_AH3$counts
##AH4
file_pathAH4 <- "AH4_Heart.quant.sf"
sample_AH4 <- data.frame(sampleName = "AH4_Heart_M", fileName = file_pathAH4)
txi_AH4 <- tximport(files = sample_AH4$fileName, type = "salmon", txOut = TRUE)
gene_counts_AH4 <- txi_AH4$counts

#combine each individual into one big table
##convert to dataframe
gene_counts_AH1 <- as.data.frame(gene_counts_AH1)
gene_counts_AH2 <- as.data.frame(gene_counts_AH2)
gene_counts_AH3 <- as.data.frame(gene_counts_AH3)
gene_counts_AH4 <- as.data.frame(gene_counts_AH4)
##create column for row names (gene names)
gene_counts_AH1 <- rownames_to_column(gene_counts_AH1, var = "Gene")
gene_counts_AH2 <- rownames_to_column(gene_counts_AH2, var = "Gene")
gene_counts_AH3 <- rownames_to_column(gene_counts_AH3, var = "Gene")
gene_counts_AH4 <- rownames_to_column(gene_counts_AH4, var = "Gene")
#rename counts columns
colnames(gene_counts_AH1)[colnames(gene_counts_AH1) == "V1"] <- "AH1_Female"
colnames(gene_counts_AH2)[colnames(gene_counts_AH2) == "V1"] <- "AH2_Female"
colnames(gene_counts_AH3)[colnames(gene_counts_AH3) == "V1"] <- "AH3_Male"
colnames(gene_counts_AH4)[colnames(gene_counts_AH4) == "V1"] <- "AH4_Male"
#merge
horse_counts_heart <- left_join(gene_counts_AH1, gene_counts_AH2, by = "Gene")
horse_counts_heart <- left_join(horse_counts_heart, gene_counts_AH3, by = "Gene")
horse_counts_heart <- left_join(horse_counts_heart, gene_counts_AH4, by = "Gene")
rownames(horse_counts_heart) <- horse_counts_heart$Gene
horse_counts_heart <- horse_counts_heart[, -1]


#create metadata
metadata_heart <- data.frame(
  title = c("AH1_Female", "AH2_Female", "AH3_Male", "AH4_Male"), 
  source = c("heart", "heart", "heart", "heart"),
  sex = c("F", "F", "M", "M")
)

rownames(metadata_heart) <- c("AH1_Female", "AH2_Female", "AH3_Male", "AH4_Male")

#aligning counts and metadata
all(colnames(horse_counts_heart) %in% rownames(metadata_heart)) ##TRUE
#check ordering
all(colnames(horse_counts_heart) == rownames(metadata_heart)) ##TRUE

#STEP 2: CREATE DESeqDataSet OBJECT
dds_horse <- DESeqDataSetFromMatrix(countData = round(horse_counts_heart), 
                                  colData = metadata_heart, 
                                  design = ~sex) 
dds_horse
##added "round" into above code b/c got error indicating that some values in assay are not integers
#prefiltering: remove rows w/ low gene counts (less than 10)
keep <- rowSums(counts(dds_horse)) >=10
dds_horse <- dds_horse[keep,] ##removed ~12000 low gene counts (from original ~44000)

#set factor level (using male as reference level)
dds_horse$sex <- relevel(dds_horse$sex, ref = "M")
dds_horse$sex ##M is listed first in levels, indicating it is reference level

#STEP 3: RUN DEseq
deseq_dds_horse <- DESeq(dds_horse)
results_dds_horse <- results(deseq_dds_horse, alpha = 0.05)
results_dds_horse

#STEP 4: Explore results
##general
summary(results_dds_horse)
resultsNames(deseq_dds_horse)

##plotting
plotMA(results_dds_horse)

#STEP 5: EXPORT DATA
##reoder so lowest padj values are first
results_dds_horse_reodered <- results_dds_horse[order(results_dds_horse$padj),]
write.csv(as.data.frame(results_dds_horse_reodered), 
          file = "horsecardiacDEseq_results.csv")

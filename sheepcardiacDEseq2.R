#Goal: use sheep count data from BLANK to calculate sex bias in gene expression using deSeq2

#install packages/libraries
##dplyr
library(dplyr)
##tidyverse
#install.packages("tidyverse")
library(tidyverse)
##DESeq2
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)


#STEP 1: DATA READING AND WRANGLING
#read in counts data
file_path <- 'E-MTAB-3838-raw-counts.tsv.undecorated'
sheep_counts_raw <- read.table(file_path, header = TRUE, sep = "\t")

#read in study data
study_data <- read.table(file = "E-MTAB-3838-experiment-design.tsv", header = TRUE, sep = "\t", fill = TRUE)
colnames(study_data)
study_data <- study_data %>% filter(Sample.Characteristic.developmental.stage. == "adult")
study_data <- study_data[grepl("ventricle", study_data$Sample.Characteristic.organism.part., ignore.case = TRUE), ]
##from this, know we need runs: ERR489258, ERR489259, ERR489268, ERR489269

##edit/pare down data (only heart columns)
sheep_counts_heart <- select(sheep_counts_raw, c("Gene.ID", "ERR489258", "ERR489259", "ERR489268", "ERR489269"))
rownames(sheep_counts_heart) <- sheep_counts_heart$Gene.ID
sheep_counts_heart <- sheep_counts_heart[, -1]

##subset metadata: only relevant columns, only heart
metadata_heart <- select(study_data, c(1,8,10))

##create column for sex 
sex_column <- c("M", "M", "F", "F")
##add sex column to metadata_heart
metadata_heart$sex <- sex_column

##rename rows in metadata_heart so they reflect column names (sample names) in macaque_counts_heart data
metadata_heart$Run
rownames(metadata_heart) <- c("ERR489258", "ERR489259", "ERR489268", "ERR489269")

#aligning counts and metadata
all(colnames(sheep_counts_heart) %in% rownames(metadata_heart)) ##TRUE
#check that order matches
all(colnames(sheep_counts_heart) == rownames(metadata_heart)) ##TRUE


#STEP 2: CREATE DESeqDataSet OBJECT
dds_sheep <- DESeqDataSetFromMatrix(countData = round(sheep_counts_heart), 
                                  colData = metadata_heart, 
                                  design = ~sex) 
dds_sheep
##added "round" into above code b/c got error indicating that some values in assay are not integers
#prefiltering: remove rows w/ low gene counts (less than 10)
keep <- rowSums(counts(dds_sheep)) >=10
dds_sheep <- dds_sheep[keep,] ##removed ~12000 low gene counts (from original ~25000)

#set factor level (using male as reference level)
dds_sheep$sex <- relevel(dds_sheep$sex, ref = "M")
dds_sheep$sex ##M is listed first in levels, indicating it is reference level

#STEP 3: RUN DEseq
deseq_dds_sheep <- DESeq(dds_sheep)
results_dds_sheep <- results(deseq_dds_sheep, alpha = 0.05)
results_dds_sheep

#STEP 4: Explore results
##general
summary(results_dds_sheep)
resultsNames(deseq_dds_sheep)

##plotting
plotMA(results_dds_sheep)

#STEP 5: EXPORT DATA
##reoder so lowest padj values are first
results_dds_sheep_reodered <- results_dds_sheep[order(results_dds_sheep$padj),]
write.csv(as.data.frame(results_dds_sheep_reodered), 
          file = "sheepcardiacDEseq_results.csv")


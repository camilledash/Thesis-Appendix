#Goal: use dog count data from Naqvi (2019) to calculate sex bias in gene expression using deSeq2

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
dog_counts_raw <- as.matrix(read.table(file = "GSE125483_dog.salmon.tximport.counts.txt"))
dog_counts_raw_df <- as.data.frame(dog_counts_raw)
##edit/pare down data (only heart columns)
colnames(dog_counts_raw_df)
dog_counts_raw_df %>% select(contains("Heart"))
dog_counts_heart <- select(dog_counts_raw_df, c(7, 20, 32, 49, 62))

#read in metadata
gse <- getGEO(GEO = "GSE125483", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[4]]))
head(metadata)
##subset metadata: only relevant columns, only heart
metadata_heart <- select(metadata, c(1,8))
metadata_heart <- metadata_heart[metadata_heart$source_name_ch1 == "Heart", ]

##create column for sex (since not in metadata): order of dogs in metadata_heart is 1088, 1510, 1170, 1593, 1324, corresponding with F, F, M, F, M (according to dog_counts_heart)
sex_column <- c("F", "F", "M", "F", "M")
##add sex column to metadata_heart
metadata_heart$sex <- sex_column
##rename rows in metadata_heart so they reflect column names (sample names) in dog_counts_heart data
colnames(dog_counts_heart)
metadata_heart$title

rownames(metadata_heart) <- c("Dog_Female_Heart_L16_1088", "Dog_Female_Heart_L16_1510", "Dog_Male_Heart_L16_1170", "Dog_Female_Heart_L16_1593", "Dog_Male_Heart_L16_1324")

#aligning counts and metadata
all(colnames(dog_counts_heart) %in% rownames(metadata_heart)) ##TRUE
##reorder columns in dog_counts_heart so its same order as rows in metadata_heart
rownames(metadata_heart)
colnames(dog_counts_heart)

dog_counts_heart_reordered <- dog_counts_heart[, c(1, 4, 2, 5, 3)]

#check that order matches
all(colnames(dog_counts_heart_reordered) == rownames(metadata_heart)) ##TRUE

#STEP 2: CREATE DESeqDataSet OBJECT
dds_dog <- DESeqDataSetFromMatrix(countData = round(dog_counts_heart_reordered), 
                                    colData = metadata_heart, 
                                    design = ~sex) 
dds_dog
##added "round" into above code b/c got error indicating that some values in assay are not integers
#prefiltering: remove rows w/ low gene counts (less than 10)
keep <- rowSums(counts(dds_dog)) >=10
dds_dog <- dds_dog[keep,] ##removed ~15000 low gene counts (from original ~40000)

#set factor level (using male as reference level)
dds_dog$sex <- relevel(dds_dog$sex, ref = "M")
dds_dog$sex ##M is listed first in levels, indicating it is reference level

#STEP 3: RUN DEseq
deseq_dds_dog <- DESeq(dds_dog)
results_dds_dog <- results(deseq_dds_dog, alpha = 0.05)
results_dds_dog

#STEP 4: Explore results
##general
summary(results_dds_dog)
resultsNames(deseq_dds_dog)

##plotting
plotMA(results_dds_dog)

#STEP 5: EXPORT DATA
##reoder so lowest padj values are first
results_dds_dog_reodered <- results_dds_dog[order(results_dds_dog$padj),]
write.csv(as.data.frame(results_dds_dog_reodered), 
          file = "dogcardiacDEseq_results.csv")



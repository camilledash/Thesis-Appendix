##Goal: determine if there are intergroup differences in differential cardiac gene expression by placenta type

#load libraries
library(tidyverse)
library(dplyr)
library(pheatmap)
#install.packages("viridis")
library(viridis)
#load data
all_log2FC_padj_cardiac <- read_csv("allspecies_log2FC_padj.csv")
rows_to_remove <- c(2, 4, 6, 8, 10, 12, 14)
all_log2FC_cardiac <- all_log2FC_padj_cardiac[-rows_to_remove, ]

final_log2FC_cardiac <- all_log2FC_cardiac[, colSums(is.na(all_log2FC_cardiac)) == 0]
final_log2FC_cardiac <- as.data.frame(final_log2FC_cardiac)

rownames(final_log2FC_cardiac) <- final_log2FC_cardiac$...1
log2FC_cardiac_hist <- final_log2FC_cardiac[, -c(1,3)]
log2FC_cardiac_interdigitation <- final_log2FC_cardiac[, -c(1,2)]
#FIX INTERDIG
log2FC_cardiac_interdigitation[3, "Placenta_Interdig"] <- "lamel"
#w/o dog???
log_interdigitation_nodog <- log2FC_cardiac_interdigitation[-3, ]


#FIND SIGNIF GENES
padj_updown <- all_log2FC_padj_cardiac[, colSums(is.na(all_log2FC_padj_cardiac)) == 0]
padj_updown <- as.data.frame(padj_updown)
rownames(padj_updown) <- padj_updown$...1
padj_updown <- select(padj_updown, -c(1:3))

deg_df <- padj_updown[, colSums(padj_updown[grep("_padj$", rownames(padj_updown)),] < 0.05) > 0]

#invasiveness
#data adjust
deg_df <- deg_df[-c(2, 4, 6, 8, 10, 12, 14), ]
deg_df_hist <- deg_df %>%
  mutate(Placenta_Hist = c("hemo", "hemo", "endo", "hemo", "hemo", "epi", "epi"))
deg_df_hist <- deg_df_hist %>%
  select(Placenta_Hist, everything())

#anova
anova_results_df <- data.frame(Gene = character(), F_value = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each column (excluding the grouping variable)
for (col in colnames(deg_df_hist)[-which(colnames(deg_df_hist) == "Placenta_Hist")]) {
  # Perform ANOVA
  anova_result <- aov(deg_df_hist[[col]] ~ Placenta_Hist, data = deg_df_hist)
  
  # Extract F-statistic and p-value
  anova_summary <- summary(anova_result)
  F_value <- anova_summary[[1]][1, 4]  # Assuming the factor variable is the first element in the summary output
  p_value <- anova_summary[[1]][1, 5]   # Assuming the factor variable is the first element in the summary output
  
  # Store results in the data frame
  result_row <- data.frame(Gene = col, F_value = F_value, p_value = p_value, stringsAsFactors = FALSE)
  anova_results_df <- rbind(anova_results_df, result_row)
}

#number total signif genes
barrier_signif <- anova_results_df[anova_results_df$p_value < 0.05, ]

#sorted
sorted_anova_results <- anova_results_df[order(anova_results_df$p_value), ]
top_ten_anova <- sorted_anova_results[1:10, ]

#tukey on the ten
# Create an empty dataframe to store Tukey's HSD results
tukey_results_df <- data.frame(Gene = character(), Group1 = character(), Group2 = character(), Estimate = numeric(), p_adj = numeric(), stringsAsFactors = FALSE)

# Loop through each row in the top_ten_anova_results dataframe
for (i in 1:nrow(top_ten_anova)) {
  # Extract information for the current gene
  gene <- top_ten_anova$Gene[i]
  F_value <- top_ten_anova$F_value[i]
  p_value <- top_ten_anova$p_value[i]
  
  # If the p-value is significant (you may adjust the threshold)
  if (p_value < 0.05) {
    # Perform Tukey's HSD
    anova_result <- aov(log2FC_cardiac_hist[[gene]] ~ Placenta_Hist, data = log2FC_cardiac_hist)
    tukey_result <- TukeyHSD(anova_result)
    
    # Extract relevant information from Tukey's HSD
    tukey_info <- as.data.frame(tukey_result$`Placenta_Hist`)
    
    # Add Gene column
    tukey_info$Gene <- gene
    
    # Store results in the tukey_results_df dataframe
    tukey_results_df <- rbind(tukey_results_df, tukey_info)
  }
}

#subset 10 most signif
top_ten_anova$Gene
subset_hist <- select(log2FC_cardiac_hist, c("PON3", "USP4", "WNK1", "KDR", "MPP7", "DES", "VWA7", "NET1", "NFKB2", "WLS"))
subset_hist_padj <- select(all_log2FC_padj_cardiac, c("...1","PON3", "USP4", "WNK1", "KDR", "MPP7", "DES", "VWA7", "NET1", "NFKB2", "WLS"))

#EXPORT
write.csv(as.data.frame(anova_results_df), file = "anovaFEB1.csv")
write.csv(as.data.frame(top_ten_anova), file = "toptenanova_FEB1.csv")
write.csv(as.data.frame(tukey_results_df), file = "toptentukey_FEB1.csv")

#interdigitation
#data adjust
deg_df_interdig <- deg_df %>%
  mutate(Placenta_Interdig = c("vil", "laby", "lamel", "vil", "laby", "vil", "vil"))
deg_df_interdig <- deg_df_interdig[-3, ]
deg_df_interdig <- deg_df_interdig %>%
  select(Placenta_Interdig, everything())

#mw
mann_whitney_results_df <- data.frame(Gene = character(), U_statistic = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Number of permutations
num_permutations <- 100000

set.seed(123)  # Set seed for reproducibility, change as needed

# Loop through each gene
for (col in colnames(deg_df_interdig)[-which(colnames(deg_df_interdig) == "Placenta_Interdig")]) {
  # Extract data for each group
  laby_data <- deg_df_interdig[deg_df_interdig$Placenta_Interdig == "laby", col]
  vil_data <- deg_df_interdig[deg_df_interdig$Placenta_Interdig == "vil", col]
  
  # Observed Mann-Whitney U statistic
  observed_statistic <- wilcox.test(laby_data, vil_data)$statistic
  
  # Initialize a vector to store permuted test statistics
  permuted_statistics <- numeric(num_permutations)
  
  # Permutation testing
  for (i in 1:num_permutations) {
    # Randomly permute the group labels
    permuted_labels <- sample(c(laby_data, vil_data))
    
    # Split the permuted data into groups
    permuted_group1 <- permuted_labels[1:length(laby_data)]
    permuted_group2 <- permuted_labels[(length(laby_data) + 1):length(permuted_labels)]
    
    # Compute the Mann-Whitney U statistic for the permuted data
    permuted_statistic <- wilcox.test(permuted_group1, permuted_group2)$statistic
    
    # Store the permuted test statistic
    permuted_statistics[i] <- permuted_statistic
  }
  
  # Calculate the p-value based on the permutation distribution
  p_value_permutation <- (sum(permuted_statistics >= observed_statistic) + 1) / (num_permutations + 1)
  
  # Store results in the data frame
  result_row <- data.frame(Gene = col, U_statistic = observed_statistic, p_value = p_value_permutation, stringsAsFactors = FALSE)
  mann_whitney_results_df <- rbind(mann_whitney_results_df, result_row)
}

#No signif results



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


#FIND GENES UPREG IN ONE GROUP AND NOT THE OTHERS
deg_df_transposed <- as.data.frame(t(deg_df))
#EPI
#upreg in horse
horseup_condition <- deg_df_transposed %>%
  filter(`Horse_padj` < 0.05 & `Horse_log2FC` > 0)
#upreg in sheep
sheepup_condition <- deg_df_transposed %>%
  filter(`Sheep_padj` < 0.05 & `Sheep_log2FC` > 0)
#UPREG IN EPI
epi_up_rows <- intersect(rownames(horseup_condition), rownames(sheepup_condition))
print(epi_up_rows) # "INPP5F"  "CEP85"   "GIPC2"   "GPBP1L1" "TAPBPL"  "EEF2"    "CTSB"   
#downreg in horse
horsedown_condition <- deg_df_transposed %>%
  filter(`Horse_padj` < 0.05 & `Horse_log2FC` < 0)
#upreg in sheep
sheepdown_condition <- deg_df_transposed %>%
  filter(`Sheep_padj` < 0.05 & `Sheep_log2FC` < 0)
#UPREG IN EPI
epi_down_rows <- intersect(rownames(horsedown_condition), rownames(sheepdown_condition))
print(epi_down_rows) # "NNT"      "XIRP2"    "CRIM1"    "NDUFS1"   "LRRC39"   "SLC25A12" "CMYA5"    "CFL2"     "SYNPO2L"  "MAP4"     "EIF3A"   "MAP1B"    "ALDH6A1"  "GLG1"     "MYH7"     "SVIL"     "TCAIM"    "EHBP1"    "GOLGA4"   "SORBS2"   "PLEC"     "RDX"     "ZNF106"   "PRUNE2"   "NR4A1"    "NQO1"     "TMEM164"  "AGL"      "TOB2"     "DSP"      "SREBF1"   "NRAP"

#HEMO
#upreg in human
humanup_condition <- deg_df_transposed %>%
  filter(`Human_padj` < 0.05 & `Human_log2FC` > 0)
#upreg in sheep
mouseup_condition <- deg_df_transposed %>%
  filter(`Mouse_padj` < 0.05 & `Mouse_log2FC` > 0)
#upreg in macaque
macaqueup_condition <- deg_df_transposed %>%
  filter(`Macaque_padj` < 0.05 & `Macaque_log2FC` > 0)
#upreg in rat
ratup_condition <- deg_df_transposed %>%
  filter(`Rat_padj` < 0.05 & `Rat_log2FC` > 0)
#UPREG IN HEMO
hemo_up_rows <- Reduce(intersect, list(rownames(humanup_condition), rownames(mouseup_condition), rownames(macaqueup_condition), rownames(ratup_condition))) #none

#downreg in human
humandown_condition <- deg_df_transposed %>%
  filter(`Human_padj` < 0.05 & `Human_log2FC` < 0)
#downreg in sheep
mousedown_condition <- deg_df_transposed %>%
  filter(`Mouse_padj` < 0.05 & `Mouse_log2FC` < 0)
#downreg in macaque
macaquedown_condition <- deg_df_transposed %>%
  filter(`Macaque_padj` < 0.05 & `Macaque_log2FC` < 0)
#downreg in rat
ratdown_condition <- deg_df_transposed %>%
  filter(`Rat_padj` < 0.05 & `Rat_log2FC` < 0)
#DOWNREG IN HEMO
hemo_down_rows <- Reduce(intersect, list(rownames(humandown_condition), rownames(mousedown_condition), rownames(macaquedown_condition), rownames(ratdown_condition))) # just KDM5D

#ENDO (DOG)
#upreg in dog
dogup_condition <- deg_df_transposed %>%
  filter(`Dog_padj` < 0.05 & `Dog_log2FC` > 0) #just PLIN3
#downreg in dog
dogdown_condition <- deg_df_transposed %>%
  filter(`Dog_padj` < 0.05 & `Dog_log2FC` < 0) #KDM5D, VWA7

#ARE THERE ANY EPI UP THAT ARE ALSO UP IN HEMOS
#human
if (any(epi_up_rows %in% rownames(humanup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} #NO
#mouse
if (any(epi_up_rows %in% rownames(mouseup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} #NO
#macaque
if (any(epi_up_rows %in% rownames(macaqueup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} #NO
#rat
if (any(epi_up_rows %in% rownames(ratup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
}
#dog
if (any(epi_up_rows %in% rownames(dogup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
}#no

#ANY DOWN EPI ALSO IN DOWN HEMO/ENDO
if (any(epi_down_rows %in% rownames(humandown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} #YES
human_epi_down_matching <- intersect(epi_down_rows, rownames(humandown_condition))
#mouse
if (any(epi_down_rows %in% rownames(mousedown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} #YES
mouse_epi_down_matching <- intersect(epi_down_rows, rownames(mousedown_condition))
#macaque
if (any(epi_down_rows %in% rownames(macaquedown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} #NO
#rat
if (any(epi_down_rows %in% rownames(ratdown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
}
#dog
if (any(epi_down_rows %in% rownames(dogdown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
}#no

#ENDO
#upregged
#human
if (any("PLIN3" %in% rownames(humanup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
}
#mouse
if (any("PLIN3" %in% rownames(mouseup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
}
#macaque
if (any("PLIN3" %in% rownames(macaqueup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
}
#rat
if (any("PLIN3" %in% rownames(humanup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
}
#horse
if (any("PLIN3" %in% rownames(horseup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
}
#sheep
if (any("PLIN3" %in% rownames(sheepup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
}
#downregged
#human
if (any(rownames(dogdown_condition) %in% rownames(humandown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} #yes
human_endo_down_matching <- intersect(rownames(dogdown_condition), rownames(mousedown_condition)) #KDM5D
#mouse
if (any("VWA7" %in% rownames(mousedown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} 
#rat
if (any("VWA7" %in% rownames(ratdown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} 
#macaque
if (any("VWA7" %in% rownames(macaquedown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} 
#horse
if (any("VWA7" %in% rownames(horsedown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} 
#sheep
if (any("VWA7" %in% rownames(sheepdown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} 



#INTERDIG
#LABY
#UPREG IN LABY
laby_up_rows <- intersect(rownames(mouseup_condition), rownames(ratup_condition))
print(laby_up_rows) #DDX3X  
#DOWNREG in LABY
laby_down_rows <- intersect(rownames(mousedown_condition), rownames(ratdown_condition))
print(laby_down_rows) #KDM5D

#UPREG IN VIL
vil_up_rows <- Reduce(intersect, list(rownames(humanup_condition), rownames(horseup_condition), rownames(macaqueup_condition), rownames(sheepup_condition))) #none
#DOWNREG IN VIL
vil_down_rows <- Reduce(intersect, list(rownames(humandown_condition), rownames(horsedown_condition), rownames(macaquedown_condition), rownames(sheepdown_condition))) #none

#uniquely upregulated in laby
#human
if (any("DDX3X" %in% rownames(humanup_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} #YES matches

#uniquely downregulated in laby
#human
if (any("KDM5D" %in% rownames(humandown_condition))) {
  print("There are matches.")
} else {
  print("No matches found.")
} #YES matches






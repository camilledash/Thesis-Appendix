##Goal: create matrix/df with genes and their log2FC for each species

#Installing packages/libraries
##dplyr
library(dplyr)
##tidyverse
library(tidyverse)
#rtracklayer (for horse data)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("rtracklayer")
library(rtracklayer)

#STEP 1: importing data
##deseq data
humanDEseq2 <- read_csv("humancardiacDEseq2_results.csv")
mouseDEseq2 <- read_csv("mousecardiacDEseq_results.csv")
dogDEseq2 <- read_csv("dogcardiacDEseq_results.csv")
macaqueDEseq2 <- read_csv("macaquecardiacDEseq_results.csv")
ratDEseq2 <- read_csv("ratcardiacDEseq_results.csv")
horseDEseq2 <- read_csv("horsecardiacDEseq_results.csv")
sheepDEseq2 <- read_csv("sheepcardiacDEseq_results.csv")

###special horse data
gff3_horse <- "all_samples.gff3"
horse_annotations <- import(gff3_horse)
horse_annotations <- as.data.frame(horse_annotations)
horse_annotations <- select(horse_annotations, c(10:12, 16))

##orthologs (naqvi)
orthologs <- as.data.frame(read.table(file = "one2oneorth_emblids.txt"))

##orthologs (horse)
horseorthologs <- as.data.frame(read.table(file = "horse_mart_export.txt", header = TRUE, sep = "\t"))
horseorthologs <- horseorthologs %>% filter(Human.homology.type == "ortholog_one2one")

##orthologs (sheep)
texel_ncbi_ortho <- as.data.frame(read.table(file = "texelncbi_mart_export.txt", header = TRUE, sep = "\t"))
ram_ncbi <- as.data.frame(read.table(file = "ramncbi_mart_export.txt", header = TRUE, sep = "\t"))

sheeporthologs <- as.data.frame(read.table(file = "ramsheep_mart_export.txt", header = TRUE, sep = "\t"))
sheeporthologs <- sheeporthologs %>% filter(Human.homology.type == "ortholog_one2one")

#STEP 2: renaming columns
humanDEseq2 <- humanDEseq2 %>%
  rename("...1" = "GeneNames", "log2FoldChange" = "Human_log2FC", "padj" = "Human_padj")
mouseDEseq2 <- mouseDEseq2 %>%
  rename("...1" = "GeneNames", "log2FoldChange" = "Mouse_log2FC", "padj" = "Mouse_padj")
dogDEseq2 <- dogDEseq2 %>%
  rename("...1" = "GeneNames", "log2FoldChange" = "Dog_log2FC", "padj" = "Dog_padj")
macaqueDEseq2 <- macaqueDEseq2 %>%
  rename("...1" = "GeneNames", "log2FoldChange" = "Macaque_log2FC", "padj" = "Macaque_padj")
ratDEseq2 <- ratDEseq2 %>%
  rename("...1" = "GeneNames", "log2FoldChange" = "Rat_log2FC", "padj" = "Rat_padj")
horseDEseq2 <- horseDEseq2 %>%
  rename("...1" = "GeneNames", "log2FoldChange" = "Horse_log2FC", "padj" = "Horse_padj")
sheepDEseq2 <- sheepDEseq2 %>%
  rename("...1" = "GeneNames", "log2FoldChange" = "Sheep_log2FC", "padj" = "Sheep_padj")

#STEP 3: removing genes not present in the orthologs dataset
##human
filtered_human <- humanDEseq2 %>%  
  filter(GeneNames %in% orthologs$Gene.name)
filtered_human$GeneNames %in% orthologs$Gene.name ##all TRUE
##mouse
filtered_mouse <- mouseDEseq2 %>%  
  filter(GeneNames %in% orthologs$Mouse.gene.stable.ID)
filtered_mouse$GeneNames %in% orthologs$Mouse.gene.stable.ID ##all TRUE
##dog
filtered_dog <- dogDEseq2 %>% 
  filter(GeneNames %in% orthologs$Dog.gene.stable.ID)
filtered_dog$GeneNames %in% orthologs$Dog.gene.stable.ID ##all TRUE
##macaque
filtered_macaque <- macaqueDEseq2 %>% 
  filter(GeneNames %in% orthologs$Crab.eating.macaque.gene.stable.ID)
filtered_macaque$GeneNames %in% orthologs$Crab.eating.macaque.gene.stable.ID ##all TRUE
##rat
filtered_rat <- ratDEseq2 %>% 
  filter(GeneNames %in% orthologs$Rat.gene.stable.ID)
filtered_rat$GeneNames %in% orthologs$Rat.gene.stable.ID ##all TRUE
##horse
###first, check w/ horse annotations to get as many real names as possible
horse_annotations$Parent <- as.character(horse_annotations$Parent)
pb_horse <- horseDEseq2[grepl("^PB", horseDEseq2$GeneNames), ] #only non ensembl values
pb_horse_parent <- left_join(pb_horse, horse_annotations, by = c("GeneNames" = "Parent"))
pb_horse_parent <- pb_horse_parent[complete.cases(pb_horse_parent$ID), ] #0 results
pb_horse_id <- left_join(pb_horse, horse_annotations, by = c("GeneNames" = "ID"))
pb_horse_id <- pb_horse_id[complete.cases(pb_horse_id$Parent), ] #worked
pb_horse_id <- pb_horse_id %>%
  mutate(GeneNames = Parent) #replace PBs with real names
pb_horse_id <- select(pb_horse_id, c(1:7)) 
#replace horse gene IDs w/ horse transcript IDs
pb_horse_id <- pb_horse_id %>% 
  filter(GeneNames %in% horseorthologs$Gene.stable.ID)
pb_horse_id <- full_join(pb_horse_id, horseorthologs, by = c("GeneNames" ="Gene.stable.ID"))
pb_horse_id <- pb_horse_id %>%
  mutate(GeneNames = Transcript.stable.ID)
pb_horse_id <- pb_horse_id[complete.cases(pb_horse_id$baseMean), ]
pb_horse_id <- select(pb_horse_id, c(1:7)) 
#add back w/ other genes
horseDEseq2 <- horseDEseq2[grepl("^EN", horseDEseq2$GeneNames), ]
horseDEseq2 <- bind_rows(horseDEseq2, pb_horse_id)
#now filter etc
filtered_horse <- horseDEseq2 %>% 
  filter(GeneNames %in% horseorthologs$Transcript.stable.ID)

##sheep
###first, join sheep deseq2 w/ oarv3.1 ncbi names
sheepncbi <- full_join(sheepDEseq2, texel_ncbi_ortho, by = c("GeneNames" = "Gene.stable.ID"))
sheepncbi <- sheepncbi[complete.cases(sheepncbi$baseMean), ]
sheepncbi <- sheepncbi[complete.cases(sheepncbi$NCBI.gene..formerly.Entrezgene..ID), ]
###now, join w/ ram ncbi and ensembl names
sheepncbi <- full_join(sheepncbi, ram_ncbi, by = "NCBI.gene..formerly.Entrezgene..ID")
sheepncbi <- sheepncbi[complete.cases(sheepncbi$baseMean), ]
sheepncbi <- sheepncbi[complete.cases(sheepncbi$Gene.stable.ID), ]
###finally, join w/ ram ensembl and human orthos
sheepncbi <- full_join(sheepncbi, sheeporthologs, by = "Gene.stable.ID")
sheepncbi <- sheepncbi[complete.cases(sheepncbi$baseMean), ]

#STEP 4: merge data with orthologs
##full join & remove unecessary rows
#human (dont need to join since already has standard ortholog name, just removing unnecessary columns)
human_stand <- dplyr::select(humanDEseq2, c(1, 3, 7))
#mouse
mouse_stand <- full_join(filtered_mouse, orthologs, by = c("GeneNames" ="Mouse.gene.stable.ID"))
mouse_stand <- dplyr::select(mouse_stand, c(3, 7, 16))
#dog
dog_stand <- full_join(filtered_dog, orthologs, by = c("GeneNames" = "Dog.gene.stable.ID"))
dog_stand <- dplyr::select(dog_stand, c(3, 7, 16))
#macaque
macaque_stand <- full_join(filtered_macaque, orthologs, by = c("GeneNames" = "Crab.eating.macaque.gene.stable.ID"))
macaque_stand <- dplyr::select(macaque_stand, c(3, 7, 16))
#rat
rat_stand <- full_join(filtered_rat, orthologs, by = c("GeneNames" = "Rat.gene.stable.ID"))
rat_stand <- dplyr::select(rat_stand, c(3, 7, 16))
#horse
horse_stand <- full_join(filtered_horse, horseorthologs, by = c("GeneNames" = "Transcript.stable.ID"))
horse_stand <- horse_stand %>%
  distinct(Human.gene.stable.ID, .keep_all = TRUE)
horse_stand <- dplyr::select(horse_stand, c(3, 7, 11))
#sheep
sheep_stand <- sheepncbi %>% 
  distinct(Human.gene.stable.ID, .keep_all = TRUE)
sheep_stand <- dplyr::select(sheep_stand, c(3, 7, 14))
sheep_stand <- sheep_stand %>% filter(Human.gene.name != "")
#print(sheep_stand$Human.gene.name[duplicated(sheep_stand$Human.gene.name)])

#STEP 5: merging data across species
allspecies_log2FC_raw <- human_stand %>%
  full_join(mouse_stand, by = c("GeneNames" = "Gene.name")) %>%
  full_join(dog_stand, by = c("GeneNames" = "Gene.name")) %>% 
  full_join(macaque_stand, by = c("GeneNames" = "Gene.name")) %>% 
  full_join(rat_stand, by = c("GeneNames" = "Gene.name")) %>% 
  full_join(horse_stand, by = c("GeneNames" = "Human.gene.name")) %>% 
  full_join(sheep_stand, by = c("GeneNames" = "Human.gene.name"))
allspecies_log2FC_raw <- allspecies_log2FC_raw %>% filter(GeneNames != "")
allspecies_log2FC_raw <- allspecies_log2FC_raw %>%
  distinct(GeneNames, .keep_all = TRUE)

View(allspecies_log2FC_raw)

#STEP 6: shape the table
##make Gene Names the row name (not just another column)
allspecies_log2FC_raw <- allspecies_log2FC_raw %>% 
  column_to_rownames(var = "GeneNames")
##make rows columns and vice versa
allspecies_log2FC_final <- data.frame(t(allspecies_log2FC_raw))
##make it a tibble
#allspecies_log2FC_final <- as_tibble(allspecies_log2FC_final)
##make columns for placenta type
placenta_hist <- c("hemo", "hemo", "hemo", "hemo", "endo", "endo", "hemo", "hemo", "hemo", "hemo", "epi", "epi", "epi", "epi")
placenta_interdig <- c("vil", "vil", "laby", "laby", "laby", "laby", "vil", "vil", "laby", "laby", "vil", "vil", "vil", "vil")
##insert columns
allspecies_log2FC_final <- allspecies_log2FC_final %>%
  mutate(Placenta_Hist = placenta_hist) %>% 
  mutate(Placenta_Interdig = placenta_interdig)
##reorder columns
allspecies_log2FC_final <- allspecies_log2FC_final %>% 
  select(Placenta_Hist, Placenta_Interdig, everything())

#STEP 7: export data
write.csv(as.data.frame(allspecies_log2FC_final),
          file = "allspecies_log2FC_padj.csv")

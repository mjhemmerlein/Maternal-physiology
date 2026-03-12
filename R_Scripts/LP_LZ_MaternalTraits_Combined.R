library(readxl)
library(dplyr)
library(edgeR)
library(limma)
library(variancePartition)
library(BiocParallel)
library(tidyverse)


# Combined data -----------
# Generates keep variable from all data

LP_LZ_Sample_Info <- read_csv("Raw_Data/MetaData_LP_Maternal_resid_scaled.csv")
LP_LZ_Sample_Info <- as.data.frame(LP_LZ_Sample_Info)
LP_LZ_Sample_Info = LP_LZ_Sample_Info[-c(80,81),]
rownames(LP_LZ_Sample_Info) <- LP_LZ_Sample_Info$Sample_ID_LZ
LP_LZ_Sample_Info <- LP_LZ_Sample_Info %>% filter(Strain == "BW" | Strain == "ME")
LP_LZ_Sample_Info <- LP_LZ_Sample_Info[LP_LZ_Sample_Info$Sample_ID_LZ != "LZ089",]
LP_LZ_Sample_Info$Group <- as.factor(LP_LZ_Sample_Info$Group)


Pman_rawreads_LP_LZ <- read_xlsx("Raw_Data/LP_LZ_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads_LP_LZ <- as.data.frame(Pman_rawreads_LP_LZ)
Pman_rawreads_LP_LZ = `row.names<-`(Pman_rawreads_LP_LZ, Pman_rawreads_LP_LZ$Geneid)
Pman_rawreads_LP_LZ <- Pman_rawreads_LP_LZ %>% filter(!is.na(Geneid))
Pman_rawreads_LP_LZ <- Pman_rawreads_LP_LZ[,-c(1:6)]
Pman_rawreads_LP_LZ <- subset(Pman_rawreads_LP_LZ, select = -c(RNA201216ZC_LZ089_S16_L001_fastp_pman_Halign_liberal.bam))

# Filter data
Pman_readcounts_LP_LZ <- as.matrix(Pman_rawreads_LP_LZ)
dPman_0_LP_LZ <- DGEList(Pman_readcounts_LP_LZ)
dPman_0_LP_LZ <- calcNormFactors(dPman_0_LP_LZ)
keep_LP_LZ <- rowSums(cpm(dPman_0_LP_LZ) > 0.5) >= 60
dPman_LP_LZ <- dPman_0_LP_LZ[keep_LP_LZ,]
dim(dPman_LP_LZ)

plotMDS(dPman_LP_LZ, col = as.numeric(LP_LZ_Sample_Info$Strain), labels = LP_LZ_Sample_Info$Strain)

# Check and match columns
Check_LP <- LP_LZ_Sample_Info$LZ_Seq_Name
colnames(Pman_rawreads_LP_LZ) == Check_LP
colnames(Pman_rawreads_LP_LZ) <- rownames(LP_LZ_Sample_Info)

# IMPORTANT: Also update column names in dPman to match
colnames(dPman_LP_LZ) <- rownames(LP_LZ_Sample_Info)

# Define your traits
traits <- c("resid_F.Mass_Gain",
            "resid_F.Mass_GD1",
            "resid_F.Mass_GD18.5",
            "resid_FoodCon",
            "resid_Lung",
            "resid_Lung_drymass",
            "resid_LV_Sep",
            "resid_RV",
            "scaled_Gluc",
            "scaled_Hb",
            "scaled_HctAv",
            "scaled_Start_AvgLean")


# Initialize storage
results_list <- list()
sig_counts <- data.frame(
  trait = character(),
  interaction = integer(),
  n_samples = integer(),
  stringsAsFactors = FALSE
)

# Set up parallel processing
param <- SnowParam(8, "SOCK", progressbar = FALSE)

# Loop through traits
for (trait in traits) {
  
  # Remove samples with NA for this trait
  na_mask <- !is.na(LP_LZ_Sample_Info[[trait]])
  LP_LZ_Sample_Info_sub <- LP_LZ_Sample_Info[na_mask, ]
  dPman_sub <- dPman_LP_LZ[, na_mask]
  
  # Ensure column names match row names
  colnames(dPman_sub) <- rownames(LP_LZ_Sample_Info_sub)
  
  # Define formula using reformulate
  form <- reformulate(c(paste0(trait, "*O2"), "(1|Mom)"))
  
  # Run voom with dream weights
  vobjDream <- voomWithDreamWeights(
    dPman_sub,
    form, 
    LP_LZ_Sample_Info_sub, 
    BPPARAM = param, 
    plot = FALSE
  )
  
  # Fit model
  fitmm <- dream(vobjDream, form, LP_LZ_Sample_Info_sub, ddf = "Kenward-Roger")
  fitmm <- eBayes(fitmm)
  
  # Extract results
  # For main trait effect
  DE_trait <- topTable(fitmm, coef = trait, sort.by = "P", n = Inf)
  
  # For hypoxia
  DE_o2 <- topTable(fitmm, coef='O22H', sort.by = "P", n = Inf)
  
  # Count significant genes
  n_sig_trait <- sum(DE_trait$adj.P.Val < 0.05)
  n_sig_o2 <- sum(DE_o2$adj.P.Val < 0.05)
  
  # Store results
  results_list[[trait]] <- list(
    DE_trait = DE_trait,
    DE_o2 = DE_o2
  )
  
  # Store counts
  sig_counts <- rbind(sig_counts, data.frame(
    trait = trait,
    traitDE = n_sig_trait,
    treatment = n_sig_o2
  ))
  
  # Save individual trait results
  # write.csv(DE_trait, file = paste0("LP_LZ_Combined_DE_trait_", trait, ".csv"), row.names = TRUE)
  # write.csv(DE_o2, file = paste0("LP_LZ_Combined_DE_O2_", trait, ".csv"), row.names = TRUE)
}

# Save summary table
write.csv(sig_counts, file = "LP_LZ_Combined_SIG_summary.csv", row.names = FALSE)




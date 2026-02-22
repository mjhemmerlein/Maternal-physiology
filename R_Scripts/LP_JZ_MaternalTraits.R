library(readxl)
library(dplyr)
library(edgeR)
library(limma)
library(variancePartition)
library(BiocParallel)
library(tidyverse)


# Combined data -----------
# Generates keep variable from all data

LP_JZ_Sample_Info <- read_csv("Raw_Data/MetaData_LP_Maternal_resid_scaled.csv")
LP_JZ_Sample_Info <- as.data.frame(LP_JZ_Sample_Info)
LP_JZ_Sample_Info = LP_JZ_Sample_Info[-c(80,81),]
rownames(LP_JZ_Sample_Info) <- LP_JZ_Sample_Info$Sample_ID_JZ
LP_JZ_Sample_Info <- LP_JZ_Sample_Info %>% filter(Strain == "BW" | Strain == "ME")
LP_JZ_Sample_Info$Group <- as.factor(LP_JZ_Sample_Info$Group)


Pman_rawreads_LP_JZ <- read_xlsx("Raw_Data/LP_JZ_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads_LP_JZ <- as.data.frame(Pman_rawreads_LP_JZ)
Pman_rawreads_LP_JZ = `row.names<-`(Pman_rawreads_LP_JZ, Pman_rawreads_LP_JZ$Geneid)
Pman_rawreads_LP_JZ <- Pman_rawreads_LP_JZ %>% filter(!is.na(Geneid))
Pman_rawreads_LP_JZ <- Pman_rawreads_LP_JZ[,-c(1:6)]

# Check and match columns
Check_LP <- LP_JZ_Sample_Info$JZ_Seq_Name
colnames(Pman_rawreads_LP_JZ) == Check_LP
colnames(Pman_rawreads_LP_JZ) <- rownames(LP_JZ_Sample_Info)


# Filter data
Pman_readcounts_LP_JZ <- as.matrix(Pman_rawreads_LP_JZ)
dPman_0_LP_JZ <- DGEList(Pman_readcounts_LP_JZ)
dPman_0_LP_JZ <- calcNormFactors(dPman_0_LP_JZ)
keep_LP_JZ <- rowSums(cpm(dPman_0_LP_JZ) > 0.5) >= 60
dPman_LP_JZ <- dPman_0_LP_JZ[keep_LP_JZ,]
dim(dPman_LP_JZ)


# Filter to BW -------
LP_JZ_Sample_Info = LP_JZ_Sample_Info %>% filter(Strain == "BW")
BW = LP_JZ_Sample_Info$JZ_Seq_Name


Pman_rawreads_LP_JZ <- read_xlsx("Raw_Data/LP_JZ_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads_LP_JZ <- as.data.frame(Pman_rawreads_LP_JZ)
Pman_rawreads_LP_JZ = `row.names<-`(Pman_rawreads_LP_JZ, Pman_rawreads_LP_JZ$Geneid)
Pman_rawreads_LP_JZ <- Pman_rawreads_LP_JZ %>% filter(!is.na(Geneid))
Pman_rawreads_LP_JZ <- Pman_rawreads_LP_JZ[,-c(1:6)]
Pman_rawreads_LP_JZ <- Pman_rawreads_LP_JZ[,names(Pman_rawreads_LP_JZ) %in% BW]
Pman_readcounts_LP_JZ <- as.matrix(Pman_rawreads_LP_JZ)

dPman_0 <- DGEList(Pman_readcounts_LP_JZ)
dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)

dPman_BW <- dPman_0[keep_LP_JZ,]
dim(dPman_BW)
plotMDS(dPman_BW, col = as.numeric(LP_JZ_Sample_Info$Strain), labels = LP_JZ_Sample_Info$Strain)

# Check and match columns
Check_LP <- LP_JZ_Sample_Info$JZ_Seq_Name
colnames(Pman_rawreads_LP_JZ) == Check_LP
colnames(Pman_rawreads_LP_JZ) <- rownames(LP_JZ_Sample_Info)

colnames(dPman_BW) <- rownames(LP_JZ_Sample_Info)

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
  na_mask <- !is.na(LP_JZ_Sample_Info[[trait]])
  LP_JZ_Sample_Info_sub <- LP_JZ_Sample_Info[na_mask, ]
  dPman_sub <- dPman_BW[, na_mask]
  
  # Ensure column names match row names
  colnames(dPman_sub) <- rownames(LP_JZ_Sample_Info_sub)
  
  # Define formula using reformulate
  form <- reformulate(c(paste0(trait, "*O2"), "(1|Mom)"))
  
  # Run voom with dream weights
  vobjDream <- voomWithDreamWeights(
    dPman_sub,
    form, 
    LP_JZ_Sample_Info_sub, 
    BPPARAM = param, 
    plot = FALSE
  )
  
  # Fit model
  fitmm <- dream(vobjDream, form, LP_JZ_Sample_Info_sub, ddf = "Kenward-Roger")
  fitmm <- eBayes(fitmm)
  
  # Extract results
  # For main trait effect
  DE_trait <- topTable(fitmm, coef = trait, sort.by = "P", n = Inf)
  
  # For interaction
  interaction_coef <- paste0(trait, ":O22H")
  DE_ixn <- topTable(fitmm, coef = interaction_coef, sort.by = "P", n = Inf)
  
  # Count significant genes
  n_sig_trait <- sum(DE_trait$adj.P.Val < 0.05)
  n_sig_ixn <- sum(DE_ixn$adj.P.Val < 0.05)
  
  # Store results
  results_list[[trait]] <- list(
    DE_trait = DE_trait,
    DE_ixn = DE_ixn,
    n_samples = sum(na_mask)
  )
  
  # Store counts
  sig_counts <- rbind(sig_counts, data.frame(
    trait = trait,
    trait_DE = n_sig_trait,
    interaction_DE = n_sig_ixn,
    n_samples = sum(na_mask)
  ))
  
  # Save individual trait results
  write.csv(DE_trait, file = paste0("LP_JZ_BW_DE_trait_", trait, ".csv"), row.names = TRUE)
  write.csv(DE_ixn, file = paste0("LP_JZ_BW_DE_ixn_", trait, ".csv"), row.names = TRUE)
}

# Save summary table
write.csv(sig_counts, file = "LP_JZ_BW_SIG_summary.csv", row.names = FALSE)

###########################



# Combined data -----------
# Generates keep variable from all data

LP_JZ_Sample_Info <- read_csv("Raw_Data/MetaData_LP_Maternal_resid_scaled.csv")
LP_JZ_Sample_Info <- as.data.frame(LP_JZ_Sample_Info)
LP_JZ_Sample_Info = LP_JZ_Sample_Info[-c(80,81),]
rownames(LP_JZ_Sample_Info) <- LP_JZ_Sample_Info$Sample_ID_JZ
LP_JZ_Sample_Info <- LP_JZ_Sample_Info %>% filter(Strain == "BW" | Strain == "ME")
LP_JZ_Sample_Info$Group <- as.factor(LP_JZ_Sample_Info$Group)


Pman_rawreads_LP_JZ <- read_xlsx("Raw_Data/LP_JZ_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads_LP_JZ <- as.data.frame(Pman_rawreads_LP_JZ)
Pman_rawreads_LP_JZ = `row.names<-`(Pman_rawreads_LP_JZ, Pman_rawreads_LP_JZ$Geneid)
Pman_rawreads_LP_JZ <- Pman_rawreads_LP_JZ %>% filter(!is.na(Geneid))
Pman_rawreads_LP_JZ <- Pman_rawreads_LP_JZ[,-c(1:6)]

# Check and match columns
Check_LP <- LP_JZ_Sample_Info$JZ_Seq_Name
colnames(Pman_rawreads_LP_JZ) == Check_LP
colnames(Pman_rawreads_LP_JZ) <- rownames(LP_JZ_Sample_Info)

# Filter data
Pman_readcounts_LP_JZ <- as.matrix(Pman_rawreads_LP_JZ)
dPman_0_LP_JZ <- DGEList(Pman_readcounts_LP_JZ)
dPman_0_LP_JZ <- calcNormFactors(dPman_0_LP_JZ)
keep_LP_JZ <- rowSums(cpm(dPman_0_LP_JZ) > 0.5) >= 60
dPman_LP_JZ <- dPman_0_LP_JZ[keep_LP_JZ,]
dim(dPman_LP_JZ)

# Filter to ME -------
LP_JZ_Sample_Info = LP_JZ_Sample_Info %>% filter(Strain == "ME")
ME = LP_JZ_Sample_Info$JZ_Seq_Name

Pman_rawreads_LP_JZ <- read_xlsx("Raw_Data/LP_JZ_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads_LP_JZ <- as.data.frame(Pman_rawreads_LP_JZ)
Pman_rawreads_LP_JZ = `row.names<-`(Pman_rawreads_LP_JZ, Pman_rawreads_LP_JZ$Geneid)
Pman_rawreads_LP_JZ <- Pman_rawreads_LP_JZ %>% filter(!is.na(Geneid))
Pman_rawreads_LP_JZ <- Pman_rawreads_LP_JZ[,-c(1:6)]
Pman_rawreads_LP_JZ <- Pman_rawreads_LP_JZ[,names(Pman_rawreads_LP_JZ) %in% ME]
Pman_readcounts_LP_JZ <- as.matrix(Pman_rawreads_LP_JZ)

dPman_0 <- DGEList(Pman_readcounts_LP_JZ)
dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)

dPman_ME <- dPman_0[keep_LP_JZ,]
dim(dPman_ME)
plotMDS(dPman_ME, col = as.numeric(LP_JZ_Sample_Info$Strain), labels = LP_JZ_Sample_Info$Strain)

# Check and match columns
Check_LP <- LP_JZ_Sample_Info$JZ_Seq_Name
colnames(Pman_rawreads_LP_JZ) == Check_LP
colnames(Pman_rawreads_LP_JZ) <- rownames(LP_JZ_Sample_Info)

colnames(dPman_ME) <- rownames(LP_JZ_Sample_Info)

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
  na_mask <- !is.na(LP_JZ_Sample_Info[[trait]])
  LP_JZ_Sample_Info_sub <- LP_JZ_Sample_Info[na_mask, ]
  dPman_sub <- dPman_ME[, na_mask]
  
  # Ensure column names match row names
  colnames(dPman_sub) <- rownames(LP_JZ_Sample_Info_sub)
  
  # Define formula using reformulate
  form <- reformulate(c(paste0(trait, "*O2"), "(1|Mom)"))
  
  # Run voom with dream weights
  vobjDream <- voomWithDreamWeights(
    dPman_sub,
    form, 
    LP_JZ_Sample_Info_sub, 
    BPPARAM = param, 
    plot = FALSE
  )
  
  # Fit model
  fitmm <- dream(vobjDream, form, LP_JZ_Sample_Info_sub, ddf = "Kenward-Roger")
  fitmm <- eBayes(fitmm)
  
  # Extract results
  # For main trait effect
  DE_trait <- topTable(fitmm, coef = trait, sort.by = "P", n = Inf)
  
  # For interaction
  interaction_coef <- paste0(trait, ":O22H")
  DE_ixn <- topTable(fitmm, coef = interaction_coef, sort.by = "P", n = Inf)
  
  # Count significant genes
  n_sig_trait <- sum(DE_trait$adj.P.Val < 0.05)
  n_sig_ixn <- sum(DE_ixn$adj.P.Val < 0.05)
  
  # Store results
  results_list[[trait]] <- list(
    DE_trait = DE_trait,
    DE_ixn = DE_ixn,
    n_samples = sum(na_mask)
  )
  
  # Store counts
  sig_counts <- rbind(sig_counts, data.frame(
    trait = trait,
    trait_DE = n_sig_trait,
    interaction_DE = n_sig_ixn,
    n_samples = sum(na_mask)
  ))
  
  # Save individual trait results
  write.csv(DE_trait, file = paste0("LP_JZ_ME_DE_trait_", trait, ".csv"), row.names = TRUE)
  write.csv(DE_ixn, file = paste0("LP_JZ_ME_DE_ixn_", trait, ".csv"), row.names = TRUE)
}

# Save summary table
write.csv(sig_counts, file = "LP_JZ_ME_SIG_summary.csv", row.names = FALSE)


### Checking individual traits


library(readxl)
library(dplyr)
library(edgeR)
library(limma)
library(variancePartition)
library(BiocParallel)
library(tidyverse)


# Combined data -----------
# Generates keep variable from all data

EP_Sample_Info <- read_csv("Raw_Data/MetaData_EP_Maternal_resid_scaled.csv")
EP_Sample_Info <- as.data.frame(EP_Sample_Info)
rownames(EP_Sample_Info) <- EP_Sample_Info$Sample_ID
EP_Sample_Info <- EP_Sample_Info %>% filter(Strain == "BW" | Strain == "ME")
EP_Sample_Info$Group <- as.factor(EP_Sample_Info$Group)
EP_Sample_Info$Group_Timepoint = paste0("EP_", EP_Sample_Info$Group)

Pman_rawreads_EP <- read_xlsx("Raw_Data/EP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads_EP <- as.data.frame(Pman_rawreads_EP)
Pman_rawreads_EP <- Pman_rawreads_EP %>% filter(!is.na(Geneid))
row.names(Pman_rawreads_EP) <- Pman_rawreads_EP$Geneid
Pman_rawreads_EP <- Pman_rawreads_EP[,-c(1:6)]

# Check and match columns
Check_EP <- EP_Sample_Info$Seq_Name
colnames(Pman_rawreads_EP) == Check_EP
colnames(Pman_rawreads_EP) <- rownames(EP_Sample_Info)

# Filter data
Pman_readcounts_EP <- as.matrix(Pman_rawreads_EP)
dPman_0_EP <- DGEList(Pman_readcounts_EP)
dPman_0_EP <- calcNormFactors(dPman_0_EP)
keep_EP <- rowSums(cpm(dPman_0_EP) > 0.5) >= 60
dPman_EP <- dPman_0_EP[keep_EP,]
dim(dPman_EP)


# Filter to BW -------
EP_Sample_Info = EP_Sample_Info %>% filter(Strain == "BW")
BW = EP_Sample_Info$Seq_Name

Pman_rawreads_EP <- read_xlsx("Raw_Data/EP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads_EP <- as.data.frame(Pman_rawreads_EP)
Pman_rawreads_EP <- Pman_rawreads_EP %>% filter(!is.na(Geneid))
row.names(Pman_rawreads_EP) <- Pman_rawreads_EP$Geneid
Pman_rawreads_EP <- Pman_rawreads_EP[,-c(1:6)]
Pman_rawreads_EP <- Pman_rawreads_EP[,names(Pman_rawreads_EP) %in% BW]
Pman_readcounts_EP <- as.matrix(Pman_rawreads_EP)

dPman_0 <- DGEList(Pman_readcounts_EP)
dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)

dPman_BW <- dPman_0[keep_EP,]
dim(dPman_BW)
plotMDS(dPman_BW, col = as.numeric(EP_Sample_Info$Strain), labels = EP_Sample_Info$Strain)

# Check and match columns
Check_EP <- EP_Sample_Info$Seq_Name
colnames(Pman_rawreads_EP) == Check_EP
colnames(Pman_rawreads_EP) <- rownames(EP_Sample_Info)

# IMPORTANT: Also update column names in dPman_BW to match
colnames(dPman_BW) <- rownames(EP_Sample_Info)

# Define your traits
traits <- c("scaled_Gluc",
            "resid_Lung_drymass")

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
  na_mask <- !is.na(EP_Sample_Info[[trait]])
  EP_Sample_Info_sub <- EP_Sample_Info[na_mask, ]
  dPman_sub <- dPman_BW[, na_mask]
  
  # Ensure column names match row names
  colnames(dPman_sub) <- rownames(EP_Sample_Info_sub)
  
  # Define formula using reformulate
  form <- reformulate(c(paste0(trait, "*O2"), "(1|Mom)"))
  
  # Run voom with dream weights
  vobjDream <- voomWithDreamWeights(
    dPman_sub,
    form, 
    EP_Sample_Info_sub, 
    BPPARAM = param, 
    plot = FALSE
  )
  
  # Fit model
  fitmm <- dream(vobjDream, form, EP_Sample_Info_sub, ddf = "Kenward-Roger")
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
    traitDE = n_sig_trait,
    interaction = n_sig_ixn,
    n_samples = sum(na_mask)
  ))
  
  # Save individual trait results
  write.csv(DE_trait, file = paste0("byhand_BW_DE_trait_", trait, ".csv"), row.names = TRUE)
  write.csv(DE_ixn, file = paste0("byhand_BW_DE_ixn_", trait, ".csv"), row.names = TRUE)
}

# Save summary table
write.csv(sig_counts, file = "byhand_BW_SIG_summary.csv", row.names = FALSE)

###########################

# Combined data -----------
# Generates keep variable from all data

EP_Sample_Info <- read_csv("Raw_Data/MetaData_EP_Maternal_resid_scaled.csv")
EP_Sample_Info <- as.data.frame(EP_Sample_Info)
rownames(EP_Sample_Info) <- EP_Sample_Info$Sample_ID
EP_Sample_Info <- EP_Sample_Info %>% filter(Strain == "BW" | Strain == "ME")
EP_Sample_Info$Group <- as.factor(EP_Sample_Info$Group)
EP_Sample_Info$Group_Timepoint = paste0("EP_", EP_Sample_Info$Group)

Pman_rawreads_EP <- read_xlsx("Raw_Data/EP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads_EP <- as.data.frame(Pman_rawreads_EP)
Pman_rawreads_EP <- Pman_rawreads_EP %>% filter(!is.na(Geneid))
row.names(Pman_rawreads_EP) <- Pman_rawreads_EP$Geneid
Pman_rawreads_EP <- Pman_rawreads_EP[,-c(1:6)]

# Check and match columns
Check_EP <- EP_Sample_Info$Seq_Name
colnames(Pman_rawreads_EP) == Check_EP
colnames(Pman_rawreads_EP) <- rownames(EP_Sample_Info)

# Filter data
Pman_readcounts_EP <- as.matrix(Pman_rawreads_EP)
dPman_0_EP <- DGEList(Pman_readcounts_EP)
dPman_0_EP <- calcNormFactors(dPman_0_EP)
keep_EP <- rowSums(cpm(dPman_0_EP) > 0.5) >= 60
dPman_EP <- dPman_0_EP[keep_EP,]
dim(dPman_EP)

# Filter to ME -------
EP_Sample_Info = EP_Sample_Info %>% filter(Strain == "ME")
ME = EP_Sample_Info$Seq_Name
Pman_rawreads_EP <- read_xlsx("Raw_Data/EP_Pman_ExtMMFrac_readcounts_Exon.xlsx")
Pman_rawreads_EP <- as.data.frame(Pman_rawreads_EP)
Pman_rawreads_EP <- Pman_rawreads_EP %>% filter(!is.na(Geneid))
row.names(Pman_rawreads_EP) <- Pman_rawreads_EP$Geneid
Pman_rawreads_EP <- Pman_rawreads_EP[,-c(1:6)]
Pman_rawreads_EP <- Pman_rawreads_EP[,names(Pman_rawreads_EP) %in% ME]
Pman_readcounts_EP <- as.matrix(Pman_rawreads_EP)
dPman_0 <- DGEList(Pman_readcounts_EP)
dPman_0 <- calcNormFactors(dPman_0)
dim(dPman_0)
dPman_ME <- dPman_0[keep_EP,]
dim(dPman_ME)
plotMDS(dPman_ME, col = as.numeric(EP_Sample_Info$Strain), labels = EP_Sample_Info$Strain)

# Check and match columns
Check_EP <- EP_Sample_Info$Seq_Name
colnames(Pman_rawreads_EP) == Check_EP
colnames(Pman_rawreads_EP) <- rownames(EP_Sample_Info)

# IMPORTANT: Also update column names in dPman_ME to match
colnames(dPman_ME) <- rownames(EP_Sample_Info)

# Define your traits
traits <- c("resid_F.Mass_GD1",
            "scaled_HctAv")

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
  na_mask <- !is.na(EP_Sample_Info[[trait]])
  EP_Sample_Info_sub <- EP_Sample_Info[na_mask, ]
  dPman_sub <- dPman_ME[, na_mask]
  
  # Ensure column names match row names
  colnames(dPman_sub) <- rownames(EP_Sample_Info_sub)
  
  # Define formula using reformulate
  form <- reformulate(c(paste0(trait, "*O2"), "(1|Mom)"))
  
  # Run voom with dream weights
  vobjDream <- voomWithDreamWeights(
    dPman_sub,
    form, 
    EP_Sample_Info_sub, 
    BPPARAM = param, 
    plot = FALSE
  )
  
  # Fit model
  fitmm <- dream(vobjDream, form, EP_Sample_Info_sub, ddf = "Kenward-Roger")
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
    interaction = n_sig_ixn,
    n_samples = sum(na_mask)
  ))
  
  # Save individual trait results
  write.csv(DE_trait, file = paste0("byhand_ME_DE_trait_", trait, ".csv"), row.names = TRUE)
  write.csv(DE_ixn, file = paste0("byhand_ME_DE_ixn_", trait, ".csv"), row.names = TRUE)
}

# Save summary table
write.csv(sig_counts, file = "byhand_ME_SIG_summary.csv", row.names = FALSE)

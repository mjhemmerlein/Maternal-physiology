library('readxl')
library('dplyr')
library(edgeR)
library(limma)
library(variancePartition)
library(BiocParallel)
library(tidyverse)


### write fxns ------
extract_third_field <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  if (length(parts) >= 3) return(parts[3]) else return(name)
}

run_dream_analysis <- function(
    strain,
    variable_of_interest,
    metadata,
    count_matrix,
    output_prefix,
    n_threads = 8
) {
  # Subset metadata to target strain
  meta_sub <- metadata %>% 
    filter(Strain == strain)
  
  # Keep only samples with non-missing values in required columns
  meta_sub <- meta_sub %>% 
    filter(!is.na(.data[[variable_of_interest]]))
  
  # Subset count matrix to matching sample IDs
  sample_ids <- rownames(meta_sub)
  dge_sub <- dge[,sample_ids]
  
  # Construct formula dynamically
  interaction_formula <- as.formula(
    paste0("~ ", variable_of_interest, "*O2 + (1|Mom)")
  )
  
  # Run voom + dream
  param <- SnowParam(n_threads, "SOCK", progressbar = TRUE)
  
  vobj <- voomWithDreamWeights(dge_sub, interaction_formula, meta_sub, BPPARAM = param, plot = TRUE)
  fit <- dream(vobj, interaction_formula, meta_sub, ddf = "Kenward-Roger")
  fit <- eBayes(fit)
  
  # Extract coefficient names
  coef_O2 <- "O22H"  
  coef_main <- variable_of_interest
  coef_interaction <- paste0(variable_of_interest, ":O22H")
  
  # Extract results
  DE_O2 <- topTable(fit, coef = coef_O2, sort.by = "P", n = Inf)
  DE_var <- topTable(fit, coef = coef_main, sort.by = "P", n = Inf)
  DE_interact <- topTable(fit, coef = coef_interaction, sort.by = "P", n = Inf)
  
  # Print summary
  cat("Significant genes (adj.P.Val < 0.05):\n")
  cat("  O2:", sum(DE_O2$adj.P.Val < 0.05), "\n")
  cat("  ", variable_of_interest, ":", sum(DE_var$adj.P.Val < 0.05), "\n")
  cat("  Interaction:", sum(DE_interact$adj.P.Val < 0.05), "\n")
  
  # Step 11: Write results
  write.csv(DE_O2, paste0(output_prefix, "_DE_O2.csv"))
  write.csv(DE_var, paste0(output_prefix, "_DE_", variable_of_interest, ".csv"))
  write.csv(DE_interact, paste0(output_prefix, "_DE_", variable_of_interest, "_x_O2.csv"))
}



### Read in data for analysis + apply common processing of read counts -----
#metadata
EP_maternal = read.csv("Raw_Data/MetaData_EP_Maternal_resid_scaled.csv", header = TRUE)
rownames(EP_maternal) = EP_maternal$Sample_ID

#counts + renaming sample columns
raw_counts <- read_xlsx("Raw_Data/EP_Pman_ExtMMFrac_readcounts_Exon.xlsx") %>% as.data.frame()
rownames(raw_counts) = raw_counts$Geneid
new_colnames <- sapply(colnames(raw_counts), extract_third_field)
colnames(raw_counts) <- new_colnames
raw_counts_cl <- raw_counts[,-c(1:6)]

#basic read count processing for filtering
raw_counts_cl <- as.matrix(raw_counts_cl)
dge <- DGEList(raw_counts_cl)

dge <- calcNormFactors(dge)
dim(dge)
keep <- rowSums(cpm(dge) > 0.5 ) >= 60
dge <- dge[keep,]
dim(dge)

### Variables to Analyze ----
variables_to_analyze <- c("resid_FoodCon", "resid_F.Mass_Gain", "resid_Lung_drymass",
                          "resid_RV", "resid_LV_Sep", "resid_Lung_H20", "scaled_Hb",
                          "scaled_Start_AvgLean", "scaled_HctAv", "scaled_Gluc")  # Add more as needed


# Initialize a data frame to store DE summaries
de_summary <- data.frame(
  VariableOfInterest = character(),
  GenesDEByVariable = integer(),
  GenesDEByO2 = integer(),
  GenesDEByTraitxO2 = integer(),
  stringsAsFactors = FALSE
)

# Loop over each variable and run analysis
for (var in variables_to_analyze) {
  prefix <- paste0("EP_ME_", var)
  
  # Run analysis
  run_dream_analysis(
    strain = "BW",
    variable_of_interest = var,
    metadata = EP_maternal,
    count_matrix = raw_counts_cl,
    output_prefix = prefix
  )
  
  # Load results
  DE_var <- read.csv(paste0(prefix, "_DE_", var, ".csv"))
  DE_O2 <- read.csv(paste0(prefix, "_DE_O2.csv"))
  DE_interact <- read.csv(paste0(prefix, "_DE_", var, "_x_O2.csv"))
  
  # Count significant genes
  n_var <- sum(DE_var$adj.P.Val < 0.05, na.rm = TRUE)
  n_o2 <- sum(DE_O2$adj.P.Val < 0.05, na.rm = TRUE)
  n_interact <- sum(DE_interact$adj.P.Val < 0.05, na.rm = TRUE)
  
  # Append to summary
  de_summary <- rbind(de_summary, data.frame(
    VariableOfInterest = var,
    GenesDEByVariable = n_var,
    GenesDEByO2 = n_o2,
    GenesDEByTraitxO2 = n_interact
  ))
}

# Save summary
write.csv(de_summary, "EP_BW_DE_Summary.csv", row.names = FALSE)

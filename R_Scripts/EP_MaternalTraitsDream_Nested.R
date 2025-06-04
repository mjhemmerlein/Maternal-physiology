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




### Run analysis ----
run_dream_analysis(
          strain = "BW",
          variable_of_interest = "scaled_Start_AvgLean",
          metadata = EP_maternal,
          count_matrix = raw_counts_cl,
          output_prefix = "EP_BW_scaled_Start_AvgLean"
        )

run_dream_analysis(
  strain = "ME",
  variable_of_interest = "scaled_Start_AvgLean",
  metadata = EP_maternal,
  count_matrix = raw_counts_cl,
  output_prefix = "EP_ME_scaled_Start_AvgLean"
)

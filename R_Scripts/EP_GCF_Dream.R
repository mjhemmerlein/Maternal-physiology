
# Packages ------------------------------------------------------------------
library('BiocParallel')
library('dplyr')
library('ggplot2')
library('lme4')
library('readxl')
library('variancePartition')
library('edgeR')
library('Matrix')

######### Combined Analysis -----------------------------------------------------
Sample_Info <- read_xlsx("Dream_Raw_Data/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info = Sample_Info %>% filter(Strain == "BW" | Strain == "ME")

# Read in Files + QC
Pman_bairdii_rawreads <- read.csv(file = "Dream_Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]

# Check read table vs sample info
Check = Sample_Info$Seq_Name
colnames(Pman_bairdii_rawreads) == Check

colnames(Pman_bairdii_rawreads) = rownames(Sample_Info)

Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)
keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
dPman_bairdii <- dPman_bairdii_0[keep,]
dim(dPman_bairdii)
plotMDS(dPman_bairdii, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Strain)

rownames(Sample_Info) == colnames(dPman_bairdii)

# Interaction
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Strain*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_ixn <- topTable( fitmm, coef='StrainME:O22H', sort.by = "P", n = Inf)

length(DE_ixn$logFC[which(DE_ixn$adj.P.Val < 0.05)])
length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Summarize Reads
library(emmeans)

Strain <- Sample_Info$Strain
O2 <- Sample_Info$O2
summary <- data.frame()
all_genes <- row.names(DE_ixn)

for (p in all_genes) {

  gene_id <- p
  test_data <- vobjDream$E[gene_id,]
  test_model <- lmer(test_data ~ Strain*O2 + (1|(Sample_Info$Mom)))
  anova(test_model)
  output <- summary(pairs(emmeans(test_model, ~ Strain*O2), adjust = "BH"))
  output_line <- output$p.value
  output_line[7] <- gene_id

  summary <- rbind(summary, output_line)
}

output_colnames <- output$contrast
colnames(summary) <- output_colnames

corrCounts <- t(vobjDream$E)
corrCounts <- as.data.frame(corrCounts)
corrCounts$Treatment <- Sample_Info$Group
corrCounts$ID <- row.names(corrCounts)
corrCounts$Check <- Sample_Info$Sample_ID
corrCounts$Check == row.names(corrCounts)

MeanCounts <- corrCounts %>%
  group_by(corrCounts$Treatment) %>%
  summarise_at(colnames(corrCounts), funs(mean(., na.rm=TRUE)))
meanCounts <- t(MeanCounts)
meanCounts <- as.data.frame(meanCounts)
colnames(meanCounts) <- meanCounts[1,]
meanCounts <- meanCounts[-1,]
meanCounts$ID <- row.names(meanCounts)

final <- merge(meanCounts, summary, by.x = "ID", by.y = c(7))

# Write Out
write.csv(final, "Dream_Output/Dream_RawFiles/EP_dreamGCF_Contrasts.csv")
write.csv(DE_strain, "Dream_Output/Dream_RawFiles/EP_dreamGCF_strainDE.csv")
write.csv(DE_o2, "Dream_Output/Dream_RawFiles/EP_dreamGCF_o2DE.csv")
write.csv(DE_ixn, "Dream_Output/Dream_RawFiles/EP_dreamGCF_ixnDE.csv")


##### BW Only ---------------------------------------------------------
Sample_Info <- read_xlsx("Raw_Data/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info = Sample_Info %>% filter(Strain == "BW")
BW = Sample_Info$Seq_Name

# Read in Files + QC
Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,names(Pman_bairdii_rawreads) %in% BW]

# Check read table vs sample info
Check = Sample_Info$Seq_Name
colnames(Pman_bairdii_rawreads) == Check

colnames(Pman_bairdii_rawreads) = rownames(Sample_Info)

Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)

dPman_bairdii_BW <- dPman_bairdii_0[keep,]
dim(dPman_bairdii_BW)
plotMDS(dPman_bairdii_BW, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Strain)

rownames(Sample_Info) == colnames(dPman_bairdii_BW)

# Interaction  
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_BW, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream(vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "Dream_Output/Dream_RawFiles/EP_BW_dreamGCF_o2DE_background.csv")


##### ME Only ---------------------------------------------------------
Sample_Info <- read_xlsx("Raw_Data/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info = Sample_Info %>% filter(Strain == "ME")
ME = Sample_Info$Seq_Name

# Read in Files + QC
Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,names(Pman_bairdii_rawreads) %in% ME]

# Check read table vs sample info
Check = Sample_Info$Seq_Name
colnames(Pman_bairdii_rawreads) == Check

colnames(Pman_bairdii_rawreads) = rownames(Sample_Info)

Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)

dPman_bairdii_ME <- dPman_bairdii_0[keep,]
dim(dPman_bairdii_ME)
plotMDS(dPman_bairdii_ME, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Strain)

rownames(Sample_Info) == colnames(dPman_bairdii_ME)

# Interaction  
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_ME, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream(vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "Dream_Output/Dream_RawFiles/EP_ME_dreamGCF_o2DE_background.csv")

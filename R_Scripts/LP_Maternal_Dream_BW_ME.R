#Packages -------------------------------------------------------------------------
library('BiocParallel')
library('dplyr')
library('ggplot2')
library('lme4')
library('readxl')
library('variancePartition')
library('edgeR')
library('tidyr')
library('RColorBrewer')

# ### BW ONLY ----------------------------------------------------
# # Hb -----------------------------------------------------------
# #Read in data
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "BW")
# LP_maternal <- LP_maternal[-which(rownames(LP_maternal) == "LZ089"),]
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$Hb), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ####
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Hb*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_Hb <- topTable( fitmm, coef='Hb', sort.by = "P", n = Inf)
# DE_Hb_o2 <- topTable( fitmm, coef='Hb:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_Hb$logFC[which(DE_Hb$adj.P.Val < 0.05)])
# length(DE_Hb_o2$logFC[which(DE_Hb_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_BW_dreamHb_o2DE.csv")
# # write.csv(DE_Hb, "LP_BW_dreamHb_DE.csv")
# # write.csv(DE_Hb_o2, "LP_BW_dreamHb_Hb_o2DE.csv")
# 
# # Hct Avg -------------------------------------------------------
# #Read in data
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "BW")
# LP_maternal <- LP_maternal[-which(rownames(LP_maternal) == "LZ089"),]
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$HctAv), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ HctAv*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_HctAv <- topTable( fitmm, coef='HctAv', sort.by = "P", n = Inf)
# DE_HctAv_o2 <- topTable( fitmm, coef='HctAv:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_HctAv$logFC[which(DE_HctAv$adj.P.Val < 0.05)])
# length(DE_HctAv_o2$logFC[which(DE_HctAv_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_BW_dreamHct_o2DE.csv")
# # write.csv(DE_HctAv, "LP_BW_dreamHct_DE.csv")
# # write.csv(DE_HctAv_o2, "LP_BW_dreamHct_Hct_o2DE.csv")
# 
# RV -------------------------------------------------------
#Read in data
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
LP_maternal = LP_maternal %>% filter(Strain == "BW")
LP_maternal <- LP_maternal[-which(rownames(LP_maternal) == "LZ089"),]

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)

Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)
keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
dPman_bairdii <- dPman_bairdii_0[keep,]
dim(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$RV), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
rownames(LP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ RV*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_RV <- topTable( fitmm, coef='RV', sort.by = "P", n = Inf)
DE_RV_o2 <- topTable( fitmm, coef='RV:O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_RV$logFC[which(DE_RV$adj.P.Val < 0.05)])
length(DE_RV_o2$logFC[which(DE_RV_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "LP_BW_dreamRV_o2DE.csv")
# write.csv(DE_RV, "LP_BW_dreamRV_DE.csv")
write.csv(DE_RV_o2, "LP_BW_dreamRV_RV_o2DE.csv")

# LV + Sept  -------------------------------------------------------
#Read in data
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
LP_maternal = LP_maternal %>% filter(Strain == "BW")
LP_maternal <- LP_maternal[-which(rownames(LP_maternal) == "LZ089"),]

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)

Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)
keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
dPman_bairdii <- dPman_bairdii_0[keep,]
dim(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$LV_Sep), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
rownames(LP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ LV_Sep*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_LV_Sep <- topTable( fitmm, coef='LV_Sep', sort.by = "P", n = Inf)
DE_LV_Sep_o2 <- topTable( fitmm, coef='LV_Sep:O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_LV_Sep$logFC[which(DE_LV_Sep$adj.P.Val < 0.05)])
length(DE_LV_Sep_o2$logFC[which(DE_LV_Sep_o2$adj.P.Val < 0.05)])

# Write Out
# write.csv(DE_o2, "LP_BW_dreamLV_o2DE.csv")
# write.csv(DE_LV_Sep, "LP_BW_dreamLV_DE.csv")
write.csv(DE_LV_Sep_o2, "LP_BW_dreamLV_LV_o2DE.csv")

# # Female mass gain ------------------------------------------------
# #Read in data
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "BW")
# LP_maternal <- LP_maternal[-which(rownames(LP_maternal) == "LZ089"),]
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$F.Mass_Gain), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ F.Mass_Gain*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_F.Mass_Gain <- topTable( fitmm, coef='F.Mass_Gain', sort.by = "P", n = Inf)
# DE_F.Mass_Gain_o2 <- topTable( fitmm, coef='F.Mass_Gain:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_F.Mass_Gain$logFC[which(DE_F.Mass_Gain$adj.P.Val < 0.05)])
# length(DE_F.Mass_Gain_o2$logFC[which(DE_F.Mass_Gain_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_BW_dreamFMassGain_o2DE.csv")
# # write.csv(DE_F.Mass_Gain, "LP_BW_dreamFMassGain_DE.csv")
# # write.csv(DE_F.Mass_Gain_o2, "LP_BW_dreamFMassGain_FMassGain_o2DE.csv")
# 
# # Glucose  -------------------------------------------------------
# #Read in data
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "BW")
# LP_maternal <- LP_maternal[-which(rownames(LP_maternal) == "LZ089"),]
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$Gluc), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Gluc*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_Gluc <- topTable( fitmm, coef='Gluc', sort.by = "P", n = Inf)
# DE_Gluc_o2 <- topTable( fitmm, coef='Gluc:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_Gluc$logFC[which(DE_Gluc$adj.P.Val < 0.05)])
# length(DE_Gluc_o2$logFC[which(DE_Gluc_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_BW_dreamGluc_o2DE.csv")
# # write.csv(DE_Gluc, "LP_BW_dreamGluc_DE.csv")
# # write.csv(DE_Gluc_o2, "LP_BW_dreamGluc_Gluc_o2DE.csv")
# 
# # Lung  -------------------------------------------------------
# #Read in data
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "BW")
# LP_maternal <- LP_maternal[-which(rownames(LP_maternal) == "LZ089"),]
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$Lung), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Lung*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_Lung <- topTable( fitmm, coef='Lung', sort.by = "P", n = Inf)
# DE_Lung_o2 <- topTable( fitmm, coef='Lung:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_Lung$logFC[which(DE_Lung$adj.P.Val < 0.05)])
# length(DE_Lung_o2$logFC[which(DE_Lung_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_BW_dreamLung_o2DE.csv")
# # write.csv(DE_Lung, "LP_BW_dreamLung_DE.csv")
# # write.csv(DE_Lung_o2, "LP_BW_dreamLung_Lung_o2DE.csv")
# 
# # Lung H20  -------------------------------------------------------
# #Read in data
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "BW")
# LP_maternal <- LP_maternal[-which(rownames(LP_maternal) == "LZ089"),]
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$Lung_H20), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Lung_H20*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_Lung_H20 <- topTable( fitmm, coef='Lung_H20', sort.by = "P", n = Inf)
# DE_Lung_H20_o2 <- topTable( fitmm, coef='Lung_H20:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_Lung_H20$logFC[which(DE_Lung_H20$adj.P.Val < 0.05)])
# length(DE_Lung_H20_o2$logFC[which(DE_Lung_H20_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_BW_dreamLungH20_o2DE.csv")
# # write.csv(DE_Lung_H20, "LP_BW_dreamLungH20_DE.csv")
# # write.csv(DE_Lung_H20_o2, "LP_BW_dreamLungH20_LungH20_o2DE.csv")
# 
# # Lung Dry Mass  -------------------------------------------------------
# #Read in data
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "BW")
# LP_maternal <- LP_maternal[-which(rownames(LP_maternal) == "LZ089"),]
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$Lung_drymass), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Lung_drymass*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_Lung_drymass <- topTable( fitmm, coef='Lung_drymass', sort.by = "P", n = Inf)
# DE_Lung_drymass_o2 <- topTable( fitmm, coef='Lung_drymass:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_Lung_drymass$logFC[which(DE_Lung_drymass$adj.P.Val < 0.05)])
# length(DE_Lung_drymass_o2$logFC[which(DE_Lung_drymass_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_BW_dreamLungDM_o2DE.csv")
# # write.csv(DE_Lung_drymass, "LP_BW_dreamLungDM_DE.csv")
# # write.csv(DE_Lung_drymass_o2, "LP_BW_dreamLungDM_LungDM_o2DE.csv")
# 
# ### ME ONLY ----------------------------------------------------
# # Hb -----------------------------------------------------------
# #Read in data
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "ME")
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$Hb), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ####
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Hb*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_Hb <- topTable( fitmm, coef='Hb', sort.by = "P", n = Inf)
# DE_Hb_o2 <- topTable( fitmm, coef='Hb:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_Hb$logFC[which(DE_Hb$adj.P.Val < 0.05)])
# length(DE_Hb_o2$logFC[which(DE_Hb_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_ME_dreamHb_o2DE.csv")
# # write.csv(DE_Hb, "LP_ME_dreamHb_DE.csv")
# # write.csv(DE_Hb_o2, "LP_ME_dreamHb_Hb_o2DE.csv")
# 
# # Hct Avg -------------------------------------------------------
# #Read in data
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "ME")
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$HctAv), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ HctAv*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_HctAv <- topTable( fitmm, coef='HctAv', sort.by = "P", n = Inf)
# DE_HctAv_o2 <- topTable( fitmm, coef='HctAv:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_HctAv$logFC[which(DE_HctAv$adj.P.Val < 0.05)])
# length(DE_HctAv_o2$logFC[which(DE_HctAv_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_ME_dreamHct_o2DE.csv")
# # write.csv(DE_HctAv, "LP_ME_dreamHct_DE.csv")
# # write.csv(DE_HctAv_o2, "LP_ME_dreamHct_Hct_o2DE.csv")
# 
# # RV -------------------------------------------------------
# #Read in data
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "ME")
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$RV), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ RV*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_RV <- topTable( fitmm, coef='RV', sort.by = "P", n = Inf)
# DE_RV_o2 <- topTable( fitmm, coef='RV:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_RV$logFC[which(DE_RV$adj.P.Val < 0.05)])
# length(DE_RV_o2$logFC[which(DE_RV_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_ME_dreamRV_o2DE.csv")
# # write.csv(DE_RV, "LP_ME_dreamRV_DE.csv")
# # write.csv(DE_RV_o2, "LP_ME_dreamRV_RV_o2DE.csv")
# 
# LV + Sept  -------------------------------------------------------
#Read in data
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
LP_maternal = LP_maternal %>% filter(Strain == "ME")

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)

Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)
keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
dPman_bairdii <- dPman_bairdii_0[keep,]
dim(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$LV_Sep), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
rownames(LP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ LV_Sep*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_LV_Sep <- topTable( fitmm, coef='LV_Sep', sort.by = "P", n = Inf)
DE_LV_Sep_o2 <- topTable( fitmm, coef='LV_Sep:O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_LV_Sep$logFC[which(DE_LV_Sep$adj.P.Val < 0.05)])
length(DE_LV_Sep_o2$logFC[which(DE_LV_Sep_o2$adj.P.Val < 0.05)])

# Write Out
# write.csv(DE_o2, "LP_ME_dreamLV_o2DE.csv")
write.csv(DE_LV_Sep, "LP_ME_dreamLV_DE.csv")
# write.csv(DE_LV_Sep_o2, "LP_ME_dreamLV_LV_o2DE.csv")

# # Female mass gain ------------------------------------------------
# #Read in data
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "ME")
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$F.Mass_Gain), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ F.Mass_Gain*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_F.Mass_Gain <- topTable( fitmm, coef='F.Mass_Gain', sort.by = "P", n = Inf)
# DE_F.Mass_Gain_o2 <- topTable( fitmm, coef='F.Mass_Gain:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_F.Mass_Gain$logFC[which(DE_F.Mass_Gain$adj.P.Val < 0.05)])
# length(DE_F.Mass_Gain_o2$logFC[which(DE_F.Mass_Gain_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_ME_dreamFMassGain_o2DE.csv")
# # write.csv(DE_F.Mass_Gain, "LP_ME_dreamFMassGain_DE.csv")
# # write.csv(DE_F.Mass_Gain_o2, "LP_ME_dreamFMassGain_FMassGain_o2DE.csv")
# 
# # Glucose  -------------------------------------------------------
# #Read in data 
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "ME")
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$Gluc), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Gluc*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_Gluc <- topTable( fitmm, coef='Gluc', sort.by = "P", n = Inf)
# DE_Gluc_o2 <- topTable( fitmm, coef='Gluc:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_Gluc$logFC[which(DE_Gluc$adj.P.Val < 0.05)])
# length(DE_Gluc_o2$logFC[which(DE_Gluc_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_ME_dreamGluc_o2DE.csv")
# # write.csv(DE_Gluc, "LP_ME_dreamGluc_DE.csv")
# # write.csv(DE_Gluc_o2, "LP_ME_dreamGluc_Gluc_o2DE.csv")
# 
# # Lung  -------------------------------------------------------
# #Read in data 
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "ME")
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$Lung), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Lung*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_Lung <- topTable( fitmm, coef='Lung', sort.by = "P", n = Inf)
# DE_Lung_o2 <- topTable( fitmm, coef='Lung:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_Lung$logFC[which(DE_Lung$adj.P.Val < 0.05)])
# length(DE_Lung_o2$logFC[which(DE_Lung_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_ME_dreamLung_o2DE.csv")
# # write.csv(DE_Lung, "LP_ME_dreamLung_DE.csv")
# # write.csv(DE_Lung_o2, "LP_ME_dreamLung_Lung_o2DE.csv")
# 
# # Lung H20  -------------------------------------------------------
# #Read in data 
# LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
# LP_maternal = as.data.frame(LP_maternal)
# LP_maternal = LP_maternal[-c(80,81),]
# rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
# LP_maternal = LP_maternal %>% filter(Strain == "ME")
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
# Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
# 
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
# colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# LP_maternal = LP_maternal[!is.na(LP_maternal$Lung_H20), ]
# dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
# rownames(LP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Lung_H20*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_Lung_H20 <- topTable( fitmm, coef='Lung_H20', sort.by = "P", n = Inf)
# DE_Lung_H20_o2 <- topTable( fitmm, coef='Lung_H20:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_Lung_H20$logFC[which(DE_Lung_H20$adj.P.Val < 0.05)])
# length(DE_Lung_H20_o2$logFC[which(DE_Lung_H20_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "LP_ME_dreamLungH20_o2DE.csv")
# # write.csv(DE_Lung_H20, "LP_ME_dreamLungH20_DE.csv")
# # write.csv(DE_Lung_H20_o2, "LP_ME_dreamLungH20_LungH20_o2DE.csv")
# 
# Lung Dry Mass  -------------------------------------------------------
#Read in data
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
LP_maternal = LP_maternal %>% filter(Strain == "ME")

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)

Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
colnames(Pman_bairdii_rawreads) = rownames(LP_maternal)
Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)
keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
dPman_bairdii <- dPman_bairdii_0[keep,]
dim(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$Lung_drymass), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)
rownames(LP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Lung_drymass*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_Lung_drymass <- topTable( fitmm, coef='Lung_drymass', sort.by = "P", n = Inf)
DE_Lung_drymass_o2 <- topTable( fitmm, coef='Lung_drymass:O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_Lung_drymass$logFC[which(DE_Lung_drymass$adj.P.Val < 0.05)])
length(DE_Lung_drymass_o2$logFC[which(DE_Lung_drymass_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "LP_ME_dreamLungDM_o2DE.csv")
write.csv(DE_Lung_drymass, "LP_ME_dreamLungDM_DE.csv")
# write.csv(DE_Lung_drymass_o2, "LP_ME_dreamLungDM_LungDM_o2DE.csv")


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

# ### BW ONLY ---------------------------------------------------
# # Hb -----------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "BW")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$Hb), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ####
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Hb*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamHb_o2DE.csv")
# # write.csv(DE_Hb, "dreamHb_DE.csv")
# # write.csv(DE_Hb_o2, "dreamHb_Hb_o2DE.csv")
# 
# 
# # Hct Avg -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "BW")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$HctAv), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ HctAv*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamHct_o2DE.csv")
# # write.csv(DE_HctAv, "dreamHct_DE.csv")
# # write.csv(DE_HctAv_o2, "dreamHct_Hct_o2DE.csv")
# 
# 
# # RV -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "BW")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$RV), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ RV*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamRV_o2DE.csv")
# # write.csv(DE_RV, "dreamRV_DE.csv")
# # write.csv(DE_RV_o2, "dreamRV_RV_o2DE.csv")
# 
# 
# # LV + Sept  -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "BW")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$LV_Sep), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ LV_Sep*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_LV_Sep <- topTable( fitmm, coef='LV_Sep', sort.by = "P", n = Inf)
# DE_LV_Sep_o2 <- topTable( fitmm, coef='LV_Sep:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_LV_Sep$logFC[which(DE_LV_Sep$adj.P.Val < 0.05)])
# length(DE_LV_Sep_o2$logFC[which(DE_LV_Sep_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "dreamLV_o2DE.csv")
# # write.csv(DE_LV_Sep, "dreamLV_DE.csv")
# # write.csv(DE_LV_Sep_o2, "dreamLV_LV_o2DE.csv")
# 
# 
# # Female mass gain ------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "BW")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$F.Mass_Gain), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ F.Mass_Gain*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamFMassGain_o2DE.csv")
# # write.csv(DE_F.Mass_Gain, "dreamFMassGain_DE.csv")
# # write.csv(DE_F.Mass_Gain_o2, "dreamFMassGain_FMassGain_o2DE.csv")
# 
# # Glucose  -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "BW")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$Gluc), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Gluc*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamGluc_o2DE.csv")
# # write.csv(DE_Gluc, "dreamGluc_DE.csv")
# # write.csv(DE_Gluc_o2, "dreamGluc_Gluc_o2DE.csv")
# 
# # Lung  -------------------------------------------------------
# #Read in data 
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "BW")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC 
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$Lung), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Lung*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# write.csv(DE_o2, "EP_BW_dreamLung_o2DE.csv")
# write.csv(DE_Lung, "EP_BW_dreamLung_DE.csv")
# write.csv(DE_Lung_o2, "EP_BW_dreamLung_Lung_o2DE.csv")
# 
# Lung H20  -------------------------------------------------------
#Read in data
EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
EP_maternal = as.data.frame(EP_maternal)
rownames(EP_maternal) = EP_maternal$Sample_ID
EP_maternal = EP_maternal %>% filter(Strain == "BW")
EP_maternal$Group <- as.factor(EP_maternal$Group)

# Read in Files + QC
Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]

colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)
keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
dPman_bairdii <- dPman_bairdii_0[keep,]
dim(dPman_bairdii)

EP_maternal = EP_maternal[!is.na(EP_maternal$Lung_H20), ]
dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]

dim(dPman_bairdii)
plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
rownames(EP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Lung_H20*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_Lung_H20 <- topTable( fitmm, coef='Lung_H20', sort.by = "P", n = Inf)
DE_Lung_H20_o2 <- topTable( fitmm, coef='Lung_H20:O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_Lung_H20$logFC[which(DE_Lung_H20$adj.P.Val < 0.05)])
length(DE_Lung_H20_o2$logFC[which(DE_Lung_H20_o2$adj.P.Val < 0.05)])

# Write Out
# write.csv(DE_o2, "EP_BW_dreamLungH20_o2DE.csv")
write.csv(DE_Lung_H20, "EP_BW_dreamLungH20_DE.csv")
# write.csv(DE_Lung_H20_o2, "EP_BW_dreamLungH20_LungH20_o2DE.csv")

# # Lung Dry Mass  -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "BW")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$Lung_drymass), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Lung_drymass*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamLungDM_drymass_o2DE.csv")
# # write.csv(DE_Lung_drymass, "dreamLungDM_DE.csv")
# # write.csv(DE_Lung_drymass_o2, "dreamLungDM_LungDM_o2DE.csv")
# 
# ### ME ONLY ---------------------------------------------------
# # Hb -----------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "ME")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$Hb), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ####
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Hb*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamHb_o2DE.csv")
# # write.csv(DE_Hb, "dreamHb_DE.csv")
# # write.csv(DE_Hb_o2, "dreamHb_Hb_o2DE.csv")
# 
# 
# # Hct Avg -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "ME")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$HctAv), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ HctAv*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamHct_o2DE.csv")
# # write.csv(DE_HctAv, "dreamHct_DE.csv")
# # write.csv(DE_HctAv_o2, "dreamHct_Hct_o2DE.csv")
# 
# 
# # RV -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "ME")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$RV), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ RV*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamRV_o2DE.csv")
# # write.csv(DE_RV, "dreamRV_DE.csv")
# # write.csv(DE_RV_o2, "dreamRV_RV_o2DE.csv")
# 
# 
# # LV + Sept  -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "ME")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$LV_Sep), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ LV_Sep*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
# fitmm = eBayes(fitmm)
# 
# DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
# DE_LV_Sep <- topTable( fitmm, coef='LV_Sep', sort.by = "P", n = Inf)
# DE_LV_Sep_o2 <- topTable( fitmm, coef='LV_Sep:O22H', sort.by = "P", n = Inf)
# 
# length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
# length(DE_LV_Sep$logFC[which(DE_LV_Sep$adj.P.Val < 0.05)])
# length(DE_LV_Sep_o2$logFC[which(DE_LV_Sep_o2$adj.P.Val < 0.05)])
# 
# # Write Out
# # write.csv(DE_o2, "dreamLV_o2DE.csv")
# # write.csv(DE_LV_Sep, "dreamLV_DE.csv")
# # write.csv(DE_LV_Sep_o2, "dreamLV_LV_o2DE.csv")
# 
# 
#  # Female mass gain ------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "ME")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$F.Mass_Gain), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ F.Mass_Gain*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamFMassGain_o2DE.csv")
# # write.csv(DE_F.Mass_Gain, "dreamFMassGain_DE.csv")
# # write.csv(DE_F.Mass_Gain_o2, "dreamFMassGain_FMassGain_o2DE.csv")
# 
# # # Glucose  -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "ME")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$Gluc), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Gluc*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamGluc_o2DE.csv")
# # write.csv(DE_Gluc, "dreamGluc_DE.csv")
# # write.csv(DE_Gluc_o2, "dreamGluc_Gluc_o2DE.csv")
# 
# # Lung  -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "ME")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$Lung), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Lung*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamLung_o2DE.csv")
# # write.csv(DE_Lung, "dreamLung_DE.csv")
# # write.csv(DE_Lung_o2, "dreamLung_Lung_o2DE.csv")
# 
#  # Lung H20  -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "ME")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$Lung_H20), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Lung_H20*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamLungH20_o2DE.csv")
# # write.csv(DE_Lung_H20, "dreamLungH20_DE.csv")
# # write.csv(DE_Lung_H20_o2, "dreamLungH20_LungH20_o2DE.csv")
# 
# # Lung Dry Mass  -------------------------------------------------------
# #Read in data
# EP_maternal = read_xlsx("Raw_Data/MetaData_EP_Maternal.xlsx")
# EP_maternal = as.data.frame(EP_maternal)
# rownames(EP_maternal) = EP_maternal$Sample_ID
# EP_maternal = EP_maternal %>% filter(Strain == "ME")
# EP_maternal$Group <- as.factor(EP_maternal$Group)
# 
# # Read in Files + QC
# Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
# Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
# Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
# 
# colnames(Pman_bairdii_rawreads) = rownames(EP_maternal)
# Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
# dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)
# 
# dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
# dim(dPman_bairdii_0)
# keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
# dPman_bairdii <- dPman_bairdii_0[keep,]
# dim(dPman_bairdii)
# 
# EP_maternal = EP_maternal[!is.na(EP_maternal$Lung_drymass), ]
# dPman_bairdii = dPman_bairdii[,rownames(EP_maternal)]
# 
# dim(dPman_bairdii)
# plotMDS(dPman_bairdii, col = as.numeric(EP_maternal$Strain), labels = EP_maternal$Strain)
# rownames(EP_maternal) == colnames(dPman_bairdii)
# 
# ###
# param = SnowParam(8, "SOCK", progressbar=TRUE)
# form <- ~ Lung_drymass*O2 + (1|Mom)
# vobjDream = voomWithDreamWeights(dPman_bairdii, form, EP_maternal, BPPARAM=param, plot = T)
# fitmm = dream( vobjDream, form, EP_maternal, ddf = "Kenward-Roger")
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
# # write.csv(DE_o2, "dreamLungDM_drymass_o2DE.csv")
# # write.csv(DE_Lung_drymass, "dreamLungDM_DE.csv")
# # write.csv(DE_Lung_drymass_o2, "dreamLungDM_LungDM_o2DE.csv")

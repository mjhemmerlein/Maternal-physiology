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

# Hb ----------------------------------------------------------- 
#Read in data 
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
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
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)

rownames(LP_maternal) == colnames(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$Hb), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
rownames(LP_maternal) == colnames(dPman_bairdii)

####
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Hb*Strain + Hb*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_Hb <- topTable( fitmm, coef='Hb', sort.by = "P", n = Inf)
DE_Hb_o2 <- topTable( fitmm, coef='Hb:O22H', sort.by = "P", n = Inf)
DE_Hb_strain <- topTable( fitmm, coef='Hb:StrainME', sort.by = "P", n = Inf)

length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_Hb$logFC[which(DE_Hb$adj.P.Val < 0.05)])
length(DE_Hb_o2$logFC[which(DE_Hb_o2$adj.P.Val < 0.05)])
length(DE_Hb_strain$logFC[which(DE_Hb_strain$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_strain, "LP_dreamHb_strainDE.csv")
write.csv(DE_o2, "LP_dreamHb_o2DE.csv")
write.csv(DE_Hb, "LP_dreamHb_DE.csv")
write.csv(DE_Hb_o2, "LP_dreamHb_Hb_o2DE.csv")
write.csv(DE_Hb_strain, "LP_dreamHb_Hb_strainDE.csv")


# Hct Avg -------------------------------------------------------
#Read in data 
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
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
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)

rownames(LP_maternal) == colnames(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$HctAv), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
rownames(LP_maternal) == colnames(dPman_bairdii)

### 
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ HctAv*Strain + HctAv*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_HctAv <- topTable( fitmm, coef='HctAv', sort.by = "P", n = Inf)
DE_HctAv_o2 <- topTable( fitmm, coef='HctAv:O22H', sort.by = "P", n = Inf)
DE_HctAv_strain <- topTable( fitmm, coef='HctAv:StrainME', sort.by = "P", n = Inf)

length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_HctAv$logFC[which(DE_HctAv$adj.P.Val < 0.05)])
length(DE_HctAv_o2$logFC[which(DE_HctAv_o2$adj.P.Val < 0.05)])
length(DE_HctAv_strain$logFC[which(DE_HctAv_strain$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_strain, "LP_dreamHct_strainDE.csv")
write.csv(DE_o2, "LP_dreamHct_o2DE.csv")
write.csv(DE_HctAv, "LP_dreamHct_DE.csv")
write.csv(DE_HctAv_o2, "LP_dreamHct_Hct_o2DE.csv")
write.csv(DE_HctAv_strain, "LP_dreamHct_Hct_strainDE.csv")


# RV -------------------------------------------------------
#Read in data 
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
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
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)

rownames(LP_maternal) == colnames(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$RV), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
rownames(LP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ RV*Strain + RV*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_RV <- topTable( fitmm, coef='RV', sort.by = "P", n = Inf)
DE_RV_o2 <- topTable( fitmm, coef='RV:O22H', sort.by = "P", n = Inf)
DE_RV_strain <- topTable( fitmm, coef='RV:StrainME', sort.by = "P", n = Inf)

length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_RV$logFC[which(DE_RV$adj.P.Val < 0.05)])
length(DE_RV_o2$logFC[which(DE_RV_o2$adj.P.Val < 0.05)])
length(DE_RV_strain$logFC[which(DE_RV_strain$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_strain, "LP_dreamRV_strainDE.csv")
write.csv(DE_o2, "LP_dreamRV_o2DE.csv")
write.csv(DE_RV, "LP_dreamRV_DE.csv")
write.csv(DE_RV_o2, "LP_dreamRV_RV_o2DE.csv")
write.csv(DE_RV_strain, "LP_dreamHb_RV_strainDE.csv")


# LV + Sept  -------------------------------------------------------
#Read in data 
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
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
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)

rownames(LP_maternal) == colnames(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$LV_Sep), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
rownames(LP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ LV_Sep*Strain + LV_Sep*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_LV_Sep <- topTable( fitmm, coef='LV_Sep', sort.by = "P", n = Inf)
DE_LV_Sep_o2 <- topTable( fitmm, coef='LV_Sep:O22H', sort.by = "P", n = Inf)
DE_LV_Sep_strain <- topTable( fitmm, coef='LV_Sep:StrainME', sort.by = "P", n = Inf)

length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_LV_Sep$logFC[which(DE_LV_Sep$adj.P.Val < 0.05)])
length(DE_LV_Sep_o2$logFC[which(DE_LV_Sep_o2$adj.P.Val < 0.05)])
length(DE_LV_Sep_strain$logFC[which(DE_LV_Sep_strain$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_strain, "LP_dreamLV_strainDE.csv")
write.csv(DE_o2, "LP_dreamLV_o2DE.csv")
write.csv(DE_LV_Sep, "LP_dreamLV_DE.csv")
write.csv(DE_LV_Sep_o2, "LP_dreamLV_LV_o2DE.csv")
write.csv(DE_LV_Sep_strain, "LP_dreamHb_LV_strainDE.csv")


# Female mass gain ------------------------------------------------
#Read in data 
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
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
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)

rownames(LP_maternal) == colnames(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$F.Mass_Gain), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
rownames(LP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ F.Mass_Gain*Strain + F.Mass_Gain*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_F.Mass_Gain <- topTable( fitmm, coef='F.Mass_Gain', sort.by = "P", n = Inf)
DE_F.Mass_Gain_o2 <- topTable( fitmm, coef='F.Mass_Gain:O22H', sort.by = "P", n = Inf)
DE_F.Mass_Gain_strain <- topTable( fitmm, coef='F.Mass_Gain:StrainME', sort.by = "P", n = Inf)

length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_F.Mass_Gain$logFC[which(DE_F.Mass_Gain$adj.P.Val < 0.05)])
length(DE_F.Mass_Gain_o2$logFC[which(DE_F.Mass_Gain_o2$adj.P.Val < 0.05)])
length(DE_F.Mass_Gain_strain$logFC[which(DE_F.Mass_Gain_strain$adj.P.Val < 0.05)])

#Write Out
write.csv(DE_strain, "LP_dreamFMassGain_strainDE.csv")
write.csv(DE_o2, "LP_dreamFMassGain_o2DE.csv")
write.csv(DE_F.Mass_Gain, "LP_dreamFMassGain_DE.csv")
write.csv(DE_F.Mass_Gain_o2, "LP_dreamFMassGain_FMassGain_o2DE.csv")
write.csv(DE_F.Mass_Gain_strain, "LP_dreamFMassGain_FMassGain_strainDE.csv")


# Glucose  -------------------------------------------------------
#Read in data 
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
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
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)

rownames(LP_maternal) == colnames(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$Gluc), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
rownames(LP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Gluc*Strain + Gluc*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_Gluc <- topTable( fitmm, coef='Gluc', sort.by = "P", n = Inf)
DE_Gluc_o2 <- topTable( fitmm, coef='Gluc:O22H', sort.by = "P", n = Inf)
DE_Gluc_strain <- topTable( fitmm, coef='Gluc:StrainME', sort.by = "P", n = Inf)

length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_Gluc$logFC[which(DE_Gluc$adj.P.Val < 0.05)])
length(DE_Gluc_o2$logFC[which(DE_Gluc_o2$adj.P.Val < 0.05)])
length(DE_Gluc_strain$logFC[which(DE_Gluc_strain$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_strain, "LP_dreamGluc_strainDE.csv")
write.csv(DE_o2, "LP_dreamGluc_o2DE.csv")
write.csv(DE_Gluc, "LP_dreamGluc_DE.csv")
write.csv(DE_Gluc_o2, "LP_dreamGluc_Gluc_o2DE.csv")
write.csv(DE_Gluc_strain, "LP_dreamGluc_Gluc_strainDE.csv")


# Lung  -------------------------------------------------------
#Read in data 
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
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
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)

rownames(LP_maternal) == colnames(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$Lung), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
rownames(LP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Lung*Strain + Lung*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_Lung <- topTable( fitmm, coef='Lung', sort.by = "P", n = Inf)
DE_Lung_o2 <- topTable( fitmm, coef='Lung:O22H', sort.by = "P", n = Inf)
DE_Lung_strain <- topTable( fitmm, coef='Lung:StrainME', sort.by = "P", n = Inf)


length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_Lung$logFC[which(DE_Lung$adj.P.Val < 0.05)])
length(DE_Lung_o2$logFC[which(DE_Lung_o2$adj.P.Val < 0.05)])
length(DE_Lung_strain$logFC[which(DE_Lung_strain$adj.P.Val < 0.05)])

# Write Out ---------------------------------------------------
write.csv(DE_strain, "LP_dreamLung_strainDE.csv")
write.csv(DE_o2, "LP_dreamLung_o2DE.csv")
write.csv(DE_Lung, "LP_dreamLung_DE.csv")
write.csv(DE_Lung_o2, "LP_dreamLung_Lung_o2DE.csv")
write.csv(DE_Lung_strain, "LP_dreamLung_Lung_strainDE.csv")


# Lung H20  -------------------------------------------------------
#Read in data 
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
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
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)

rownames(LP_maternal) == colnames(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$Lung_H20), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
rownames(LP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Lung_H20*Strain + Lung_H20*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_Lung_H20 <- topTable( fitmm, coef='Lung_H20', sort.by = "P", n = Inf)
DE_Lung_H20_o2 <- topTable( fitmm, coef='Lung_H20:O22H', sort.by = "P", n = Inf)
DE_Lung_H20_strain <- topTable( fitmm, coef='Lung_H20:StrainME', sort.by = "P", n = Inf)

length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_Lung_H20$logFC[which(DE_Lung_H20$adj.P.Val < 0.05)])
length(DE_Lung_H20_o2$logFC[which(DE_Lung_H20_o2$adj.P.Val < 0.05)])
length(DE_Lung_H20_strain$logFC[which(DE_Lung_H20_strain$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_strain, "LP_dreamLungH20_strainDE.csv")
write.csv(DE_o2, "LP_dreamLungH20_o2DE.csv")
write.csv(DE_Lung_H20, "LP_dreamLungH20_DE.csv")
write.csv(DE_Lung_H20_o2, "LP_dreamLungH20_LungH20_o2DE.csv")
write.csv(DE_Lung_H20_strain, "LP_dreamLungH20_LungH20_strainDE.csv")


# Lung Dry Mass  -------------------------------------------------------
#Read in data 
LP_maternal <- read_xlsx("Raw_Data/MetaData_LP_Maternal.xlsx")
LP_maternal = as.data.frame(LP_maternal)
LP_maternal = LP_maternal[-c(80,81),]
rownames(LP_maternal) = LP_maternal$Sample_ID_LZ
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
plotMDS(dPman_bairdii, col = as.numeric(LP_maternal$Strain), labels = LP_maternal$Strain)

rownames(LP_maternal) == colnames(dPman_bairdii)

LP_maternal = LP_maternal[!is.na(LP_maternal$Lung_drymass), ]
dPman_bairdii = dPman_bairdii[,rownames(LP_maternal)]

dim(dPman_bairdii)
rownames(LP_maternal) == colnames(dPman_bairdii)

###
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Lung_drymass*Strain + Lung_drymass*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii, form, LP_maternal, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, LP_maternal, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_strain <- topTable( fitmm, coef='StrainME', sort.by = "P", n = Inf, )
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_Lung_drymass <- topTable( fitmm, coef='Lung_drymass', sort.by = "P", n = Inf)
DE_Lung_drymass_o2 <- topTable( fitmm, coef='Lung_drymass:O22H', sort.by = "P", n = Inf)
DE_Lung_drymass_strain <- topTable( fitmm, coef='Lung_drymass:StrainME', sort.by = "P", n = Inf)

length(DE_strain$logFC[which(DE_strain$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_Lung_drymass$logFC[which(DE_Lung_drymass$adj.P.Val < 0.05)])
length(DE_Lung_drymass_o2$logFC[which(DE_Lung_drymass_o2$adj.P.Val < 0.05)])
length(DE_Lung_drymass_strain$logFC[which(DE_Lung_drymass_strain$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_strain, "LP_dreamLungDM_strainDE.csv")
write.csv(DE_o2, "LP_dreamLungDM_drymass_o2DE.csv")
write.csv(DE_Lung_drymass, "LP_dreamLungDM_DE.csv")
write.csv(DE_Lung_drymass_o2, "LP_dreamLungDM_LungDM_o2DE.csv")
write.csv(DE_Lung_drymass_strain, "LP_dreamLungDM_LungDM_strainDE.csv")



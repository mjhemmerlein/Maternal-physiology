
# Packages ------------------------------------------------------------------
library('BiocParallel')
library('dplyr')
library('ggplot2')
library('lme4')
library('readxl')
library('variancePartition')
library('edgeR')
library('Matrix')

### EARLY PREGNANCY --------------------------------------------------
Pman_bairdii_rawreads <- read.csv(file = "Raw_Data/EP_GCF_readcounts.csv", header = TRUE)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]

Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)
keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 60
dPman_bairdii <- dPman_bairdii_0[keep,]
dim(dPman_bairdii)

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

colnames(Pman_bairdii_rawreads) = rownames(Sample_Info)

Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)

dPman_bairdii_BW <- dPman_bairdii_0[keep,]
dim(dPman_bairdii_BW)
plotMDS(dPman_bairdii_BW, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Strain)

Sample_Info = Sample_Info[!is.na(Sample_Info$Sex), ]
dPman_bairdii_BW = dPman_bairdii_BW[,rownames(Sample_Info)]

rownames(Sample_Info) == colnames(dPman_bairdii_BW)

# Interaction  
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Sex*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_BW, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream(vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_sex <- topTable( fitmm, coef='SexM', sort.by = "P", n = Inf)
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_ixn <- topTable( fitmm, coef='SexM:O22H', sort.by = "P", n = Inf)

length(DE_sex$logFC[which(DE_sex$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_ixn$logFC[which(DE_ixn$adj.P.Val < 0.05)])

# Write Out
# write.csv(DE_sex, "Dream_Output/Dream_RawFiles/EP_BW_Sex_sex.csv")
# write.csv(DE_o2, "Dream_Output/Dream_RawFiles/EP_BW_Sex_o2.csv")
# write.csv(DE_ixn, "Dream_Output/Dream_RawFiles/EP_BW_Sex_ixn.csv")

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

colnames(Pman_bairdii_rawreads) = rownames(Sample_Info)

Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)

dPman_bairdii_ME <- dPman_bairdii_0[keep,]
dim(dPman_bairdii_ME)
plotMDS(dPman_bairdii_ME, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Strain)

Sample_Info = Sample_Info[!is.na(Sample_Info$Sex), ]
dPman_bairdii_ME = dPman_bairdii_ME[,rownames(Sample_Info)]

rownames(Sample_Info) == colnames(dPman_bairdii_ME)

# Interaction  
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Sex*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_ME, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream(vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_sex <- topTable( fitmm, coef='SexM', sort.by = "P", n = Inf)
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_ixn <- topTable( fitmm, coef='SexM:O22H', sort.by = "P", n = Inf)

length(DE_sex$logFC[which(DE_sex$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_ixn$logFC[which(DE_ixn$adj.P.Val < 0.05)])

# Write Out
# write.csv(DE_sex, "Dream_Output/Dream_RawFiles/EP_ME_Sex_sex.csv")
# write.csv(DE_o2, "Dream_Output/Dream_RawFiles/EP_ME_Sex_o2.csv")
# write.csv(DE_ixn, "Dream_Output/Dream_RawFiles/EP_ME_Sex_ixn.csv")

##### LATE PREGNANCY ----------------------------------------------------
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)

Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
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

#### BW ONLY ---------------------------------
Sample_Info <- read_xlsx("Raw_Data/MetaData_LP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
Sample_Info = Sample_Info[-c(80,81),]
rownames(Sample_Info) = Sample_Info$Sample_ID_LZ
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "LZ089"),]
Sample_Info = Sample_Info %>% filter(Strain == "BW")
BW = Sample_Info$Seq_Name

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,names(Pman_bairdii_rawreads) %in% BW]

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
form <- ~ Sex*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_BW, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_sex <- topTable( fitmm, coef='SexM', sort.by = "P", n = Inf)
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_ixn <- topTable( fitmm, coef='SexM:O22H', sort.by = "P", n = Inf)

length(DE_sex$logFC[which(DE_sex$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_ixn$logFC[which(DE_ixn$adj.P.Val < 0.05)])

# Write Out
# write.csv(DE_sex, "Dream_Output/Dream_RawFiles/LP_BW_Sex_sex.csv")
# write.csv(DE_o2, "Dream_Output/Dream_RawFiles/LL_BW_Sex_o2.csv")
# write.csv(DE_ixn, "Dream_Output/Dream_RawFiles/LL_BW_Sex_ixn.csv")

##### ME ONLY --------------------------------------
Sample_Info <- read_xlsx("Raw_Data/MetaData_LP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
Sample_Info = Sample_Info[-c(80,81),]
rownames(Sample_Info) = Sample_Info$Sample_ID_LZ
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "LZ089"),]
Sample_Info = Sample_Info %>% filter(Strain == "ME")
ME = Sample_Info$Seq_Name

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/LP_LZ_GCF_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA201216ZC_LZ089_S16_L001_fastp.bb.trimmed_pman_Halign_liberal.bam))
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,names(Pman_bairdii_rawreads) %in% ME]

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
form <- ~ Sex*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_ME, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_sex <- topTable( fitmm, coef='SexM', sort.by = "P", n = Inf)
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_ixn <- topTable( fitmm, coef='SexM:O22H', sort.by = "P", n = Inf)

length(DE_sex$logFC[which(DE_sex$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_ixn$logFC[which(DE_ixn$adj.P.Val < 0.05)])

# Write Out
# write.csv(DE_sex, "Dream_Output/Dream_RawFiles/LP_ME_Sex_sex.csv")
# write.csv(DE_o2, "Dream_Output/Dream_RawFiles/LL_ME_Sex_o2.csv")
# write.csv(DE_ixn, "Dream_Output/Dream_RawFiles/LL_ME_Sex_ixn.csv")

######### LEUCOPUS GCF -----------------------------------------------------
Sample_Info <- read_xlsx("Raw_Data/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W237"),]
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W311"),]
Sample_Info = Sample_Info %>% filter(Strain == "LL")

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/Pman_GCF_LL_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W237_S2_L001_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W311_S33_L002_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))

colnames(Pman_bairdii_rawreads) = rownames(Sample_Info)
Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)
keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 23.25
dPman_bairdii_LL <- dPman_bairdii_0[keep,]
dim(dPman_bairdii_LL)
plotMDS(dPman_bairdii_LL, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Sample_ID)

rownames(Sample_Info) == colnames(dPman_bairdii_LL)

# Interaction
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Sex*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_LL, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_sex <- topTable( fitmm, coef='SexM', sort.by = "P", n = Inf)
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_ixn <- topTable( fitmm, coef='SexM:O22H', sort.by = "P", n = Inf)

length(DE_sex$logFC[which(DE_sex$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_ixn$logFC[which(DE_ixn$adj.P.Val < 0.05)])

# Write Out
# write.csv(DE_sex, "Dream_Output/Dream_RawFiles/LL_GCF_Sex_sex.csv")
# write.csv(DE_o2, "Dream_Output/Dream_RawFiles/LL_GCF_Sex_o2.csv")
# write.csv(DE_ixn, "Dream_Output/Dream_RawFiles/LL_GCF_Sex_ixn.csv")

######### LEUCOPUS LIFTOFF -----------------------------------------------------

Sample_Info <- read_xlsx("Raw_Data/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W237"),]
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W311"),]
Sample_Info = Sample_Info %>% filter(Strain == "LL")

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/Pman_LiftOff_LL_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W237_S2_L001_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W311_S33_L002_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))

colnames(Pman_bairdii_rawreads) = rownames(Sample_Info)
Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)
keep <- rowSums(cpm(dPman_bairdii_0) > 0.5 ) >= 23.25
dPman_bairdii_LL <- dPman_bairdii_0[keep,]
dim(dPman_bairdii_LL)
plotMDS(dPman_bairdii_LL, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Sample_ID)

rownames(Sample_Info) == colnames(dPman_bairdii_LL)

# Interaction
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ Sex*O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_LL, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_sex <- topTable( fitmm, coef='SexM', sort.by = "P", n = Inf)
DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)
DE_ixn <- topTable( fitmm, coef='SexM:O22H', sort.by = "P", n = Inf)

length(DE_sex$logFC[which(DE_sex$adj.P.Val < 0.05)])
length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])
length(DE_ixn$logFC[which(DE_ixn$adj.P.Val < 0.05)])

# Write Out
# write.csv(DE_sex, "Dream_Output/Dream_RawFiles/LL_Lift_Sex_sex.csv")
# write.csv(DE_o2, "Dream_Output/Dream_RawFiles/LL_Lift_Sex_o2.csv")
# write.csv(DE_ixn, "Dream_Output/Dream_RawFiles/LL_Lift_Sex_ixn.csv")

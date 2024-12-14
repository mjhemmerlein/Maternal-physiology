
# Packages ------------------------------------------------------------------
library('BiocParallel')
library('dplyr')
library('ggplot2')
library('lme4')
library('readxl')
library('variancePartition')
library('edgeR')
library('Matrix')

######### LEUCOPUS GCF -----------------------------------------------------

# Combined
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
form <- ~ O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_LL, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "Dream_Output/Dream_RawFiles/LL_dreamGCF_o2DE.csv")

##### Females Only ---------------------------------------------------------

Sample_Info <- read_xlsx("Raw_Data/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W237"),]
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W311"),]
Sample_Info = Sample_Info %>% filter(Strain == "LL")
Sample_Info = Sample_Info %>% filter(Sex == "F")
FEM = Sample_Info$Seq_Name

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/Pman_GCF_LL_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W237_S2_L001_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W311_S33_L002_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,names(Pman_bairdii_rawreads) %in% FEM]

colnames(Pman_bairdii_rawreads) = rownames(Sample_Info)
Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)

dPman_bairdii_LL <- dPman_bairdii_0[keep,]
dim(dPman_bairdii_LL)
plotMDS(dPman_bairdii_LL, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Sample_ID)

rownames(Sample_Info) == colnames(dPman_bairdii_LL)

# Interaction
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_LL, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "Dream_Output/Dream_RawFiles/LL_FEM_dreamGCF_o2DE.csv")

##### Males Only ---------------------------------------------------------

Sample_Info <- read_xlsx("Raw_Data/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W237"),]
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W311"),]
Sample_Info = Sample_Info %>% filter(Strain == "LL")
Sample_Info = Sample_Info %>% filter(Sex == "M")
MALE = Sample_Info$Seq_Name

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/Pman_GCF_LL_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W237_S2_L001_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W311_S33_L002_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,names(Pman_bairdii_rawreads) %in% MALE]

colnames(Pman_bairdii_rawreads) = rownames(Sample_Info)
Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)

dPman_bairdii_LL <- dPman_bairdii_0[keep,]
dim(dPman_bairdii_LL)
plotMDS(dPman_bairdii_LL, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Sample_ID)

rownames(Sample_Info) == colnames(dPman_bairdii_LL)

# Interaction
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_LL, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "Dream_Output/Dream_RawFiles/LL_MALE_dreamGCF_o2DE.csv")

######### LEUCOPUS LIFTOFF -----------------------------------------------------

# Combined
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
form <- ~ O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_LL, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "Dream_Output/Dream_RawFiles/LL_dreamLift_o2DE.csv")

##### Females Only ---------------------------------------------------------

Sample_Info <- read_xlsx("Raw_Data/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W237"),]
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W311"),]
Sample_Info = Sample_Info %>% filter(Strain == "LL")
Sample_Info = Sample_Info %>% filter(Sex == "F")
FEM = Sample_Info$Seq_Name

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/Pman_LiftOff_LL_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W237_S2_L001_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W311_S33_L002_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,names(Pman_bairdii_rawreads) %in% FEM]

colnames(Pman_bairdii_rawreads) = rownames(Sample_Info)
Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)

dPman_bairdii_LL <- dPman_bairdii_0[keep,]
dim(dPman_bairdii_LL)
plotMDS(dPman_bairdii_LL, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Sample_ID)

rownames(Sample_Info) == colnames(dPman_bairdii_LL)

# Interaction
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_LL, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "Dream_Output/Dream_RawFiles/LL_FEM_dreamLift_o2DE.csv")

##### Males Only ---------------------------------------------------------

Sample_Info <- read_xlsx("Raw_Data/MetaData_EP.xlsx")
Sample_Info = as.data.frame(Sample_Info)
rownames(Sample_Info) = Sample_Info$Sample_ID
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W237"),]
Sample_Info <- Sample_Info[-which(rownames(Sample_Info) == "W311"),]
Sample_Info = Sample_Info %>% filter(Strain == "LL")
Sample_Info = Sample_Info %>% filter(Sex == "M")
MALE = Sample_Info$Seq_Name

# Read in Files + QC
Pman_bairdii_rawreads <- read_xlsx("Raw_Data/Pman_LiftOff_LL_readcounts.xlsx")
Pman_bairdii_rawreads = as.data.frame(Pman_bairdii_rawreads)
Pman_bairdii_rawreads <- `row.names<-`(Pman_bairdii_rawreads, Pman_bairdii_rawreads$Geneid)
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,-c(1:6)]
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W237_S2_L001_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))
Pman_bairdii_rawreads <- subset(Pman_bairdii_rawreads, select = -c(RNA220328ZC_W311_S33_L002_fastp.bb.trimmed_pmanLL_Halign_liberal.bam))
Pman_bairdii_rawreads <- Pman_bairdii_rawreads[,names(Pman_bairdii_rawreads) %in% MALE]

colnames(Pman_bairdii_rawreads) = rownames(Sample_Info)
Pman_bairdii_readcounts <- as.matrix(Pman_bairdii_rawreads)
dPman_bairdii_0 <- DGEList(Pman_bairdii_readcounts)

dPman_bairdii_0 <- calcNormFactors(dPman_bairdii_0)
dim(dPman_bairdii_0)

dPman_bairdii_LL <- dPman_bairdii_0[keep,]
dim(dPman_bairdii_LL)
plotMDS(dPman_bairdii_LL, col = as.numeric(Sample_Info$Strain), labels = Sample_Info$Sample_ID)

rownames(Sample_Info) == colnames(dPman_bairdii_LL)

# Interaction
param = SnowParam(8, "SOCK", progressbar=TRUE)
form <- ~ O2 + (1|Mom)
vobjDream = voomWithDreamWeights(dPman_bairdii_LL, form, Sample_Info, BPPARAM=param, plot = T)
fitmm = dream( vobjDream, form, Sample_Info, ddf = "Kenward-Roger")
fitmm = eBayes(fitmm)

DE_o2 <- topTable( fitmm, coef='O22H', sort.by = "P", n = Inf)

length(DE_o2$logFC[which(DE_o2$adj.P.Val < 0.05)])

# Write Out
write.csv(DE_o2, "Dream_Output/Dream_RawFiles/LL_MALE_dreamLift_o2DE.csv")


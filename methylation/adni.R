# Edits by Brian Lee (leebn@sas)
# Source code: https://github.com/TransBioInfoLab/ad-meta-analysis
# ---
#   title: "ADNI dataset"
# author: "Lanyu Zhang, Tiago C. Silva, Lissette Gomez, Lily Wang, (Edits: Brian Lee -- leebn@sas)"
# date: "`r Sys.Date()`"
# output:
#   rmarkdown::html_document:
#   theme: lumen
# toc: true
# number_sections: true
# df_print: paged
# code_download: false
# toc_float:
#   collapsed: yes
# toc_depth: 3
# editor_options:
#   chunk_output_type: inline    
# ---
  
# knitr::opts_chunk$set(echo = TRUE,warning = FALSE)

# Data retrival
library(RCurl)
library(ggpubr)
library(minfi)
library(dplyr)
cohort <- "ADNI"

# TODO: Edit this later and make these directories. 
data.dir <- "~/data/ROSMAPmethNew"
data.dir.table <- "~/data/ROSMAPmethNew/datatable" 
data.dir.raw <- file.path(data.dir,"step1_download/") 
data.dir.raw.idat <- "~/data/ROSMAPmeth/"
data.dir.raw.metadata <- file.path(data.dir.raw, "Metadata")
data.dir.read <- file.path(data.dir,"step2_read_minfi/") 
data.dir.bsfilter <- file.path(data.dir,"step3_bsConvFilter/") 
data.dir.clinical.filter <- file.path(data.dir,"step4_clinical_available_filtering/") 
data.dir.probes.qc <- file.path(data.dir,"step5_probesQC_filtering/") 
data.dir.probes.normalization <- file.path(data.dir,"step6_normalization/") 
data.dir.pca <- file.path(data.dir,"step7_pca_filtering/") 
data.dir.neuron <- file.path(data.dir,"step8_neuron_comp/") 
data.dir.single.cpg.pval <- file.path(data.dir,"step9_single_cpg_pval/") 
data.dir.residuals <- file.path(data.dir,"step10_residuals/") 
data.dir.median <- file.path(data.dir,"step11_median/") 
for(p in grep("data",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

## Download data (Done and in CUBIC)

# Data Pre-processing

# Description: 
#   
# - Read in idat files
# - Remove duplicated samples
# 
# Input: 
#   
#   - idat files
# 
# Output: 
#   
#   - RGSet.RDS
# - phenoData
# - BetaMatrixRaw

RGSet <- read.metharray.exp(
  base = data.dir.raw.idat,
  recursive = TRUE,
  verbose = TRUE
)
saveRDS(RGSet, paste0(data.dir.read, "RGSet.RDS"))

##### 5. BMIQ ##################################################################
library(wateRmelon)
library(RPMM)
library(sesame)
library(sesameData)

RGSet <- readRDS(file = paste0(data.dir.read, "RGSet.RDS"))
bs <- data.frame(bisulfiteConversion = bscon(RGSet))
bsFilteredOut <- row.names(bs)[bs$bisulfiteConversion < 88]
nb.samples <- ncol(RGSet)
RGSet <- RGSet[,!colnames(RGSet) %in% bsFilteredOut]
nb.samples.bc.filtered <-  ncol(RGSet)

save(RGSet,
     nb.samples,
     bs,
     nb.samples.bc.filtered, 
     file = paste0(data.dir.bsfilter, "RGSet_bsfiltered.rda"))

load(file = paste0(data.dir.bsfilter, "RGSet_bsfiltered.rda"))

ggpubr::gghistogram(bs$bisulfiteConversion,xlab = "bisulfite Conversion" )

load(file = paste0(data.dir.bsfilter, "RGSet_bsfiltered.rda"))
RGSet <- RGSet[,colnames(RGSet) %in% phenoData$Sentrix]
phenoData <- phenoData[phenoData$Sentrix %in% colnames(RGSet),]
colnames(RGSet) %in% paste0(phenoData$Sentrix_ID,"_",phenoData$Sentrix_Position) %>% table
phenoData <- phenoData[match(colnames(RGSet),phenoData$Sentrix),]

betaSet <- getBeta(RGSet)
identical(colnames(betaSet),phenoData$Sentrix) ##TRUE
colnames(betaSet) <- phenoData$Sample
save(betaSet,
     RGSet,
     phenoData,
     file =  paste0(data.dir.clinical.filter, "/ADNIbetaMatrixraw737ind.rda"))

# Sample subset

load(paste0(data.dir.clinical.filter, "/ADNIbetaMatrixraw737ind.rda"))
nb.samples.with.slide <- ncol(betaSet)
nb.probes <- nrow(betaSet)
#--------------------------------------------------
# Read metadata
#--------------------------------------------------
clinical <- readr::read_csv(
  file.path(data.dir.raw.metadata,"ROSMAP_Clinical_2019-05_v3.csv"), # TODO: Check if this exists. 
  col_types = readr::cols())
# change 90+ to 90
clinical$age_death <- as.numeric(gsub("\\+","",clinical$age_death))
biospecimen <- readr::read_csv(
  file.path(data.dir.raw.metadata,"ROSMAP_biospecimen_metadata.csv"), # TODO: Find this file's name/path?
  col_types = readr::cols()
)
idkey <- readr::read_csv(
  file.path(data.dir.raw.metadata,"ROSMAP_IDkey.csv"), # TODO: What's this?
  col_types = readr::cols()
)
#--------------------------------------------------
# keep only one entries with DNA methylation
idkey <- idkey[idkey$mwas_id %in% colnames(betaSet),]
# keep only one entry for each DNA methylation file 
idkey <- idkey[!duplicated(idkey$mwas_id),]
# keep samples with only Methylation, braaksc 
clinical <- unique(merge (clinical, idkey, by = "projid"))
clinical <- clinical[!is.na(clinical$mwas_id),]
clinical <- clinical[!is.na(clinical$braaksc),]
# merge phenoData with clinical data
dim(phenoData) #743  10
# data in both clinical & pheno data
phenoData <- merge(phenoData, clinical, by.x = "Sample", by.y = "mwas_id")
phenoData <- unique(subset(phenoData, select = Sample:cogdx))
betaSet <- betaSet[,colnames(betaSet) %in% phenoData$Sample]
phenoData <- phenoData[match(colnames(betaSet),phenoData$Sample),]
RGSet <- RGSet[,colnames(RGSet) %in% phenoData$Sentrix]
RGSet <- RGSet[,match(phenoData$Sentrix,colnames(RGSet))]
identical(colnames(betaSet),phenoData$Sample)
identical(colnames(RGSet),phenoData$Sentrix)

save(betaSet,
     phenoData,
     RGSet,
     file =  paste0(data.dir.clinical.filter, "/ROSMAPbetaMatrixraw734ind.rda"))

load(paste0(data.dir.clinical.filter, "/ROSMAPbetaMatrixraw734ind.rda"))
nb.samples.with.clinical <- ncol(betaSet)
dim(phenoData)
dim(betaSet)

detP <- detectionP(RGSet, type = "m+u")
failed.01 <- detP > 0.01
passedProbes <- rownames(failed.01)[rowMeans(failed.01) == 0]
sum(rowMeans(failed.01) == 0)  
betaSet <- betaSet[passedProbes, ]
dim(betaSet)
nb.probes.detectP <- nrow(betaSet)
###### (b) keep only probes that start with "cg"
betaSet <- subset (betaSet, substr(row.names(betaSet),1,2) == "cg")
dim(betaSet)
nb.probes.detectP.cg <- nrow(betaSet)

##### 2. drop probes that are on X/Y ###########################################
##### 3. drop probes where SNP with MAF >= 0.01 in the last 5 bp of the probe ##
library(DMRcate)

betaSet <- rmSNPandCH(object = betaSet,
                      dist = 5, 
                      mafcut = 0.01, 
                      and = TRUE,
                      rmcrosshyb = FALSE,
                      rmXY = TRUE)
nb.probes.cg.dmrcate <- nrow(betaSet)
dim(betaSet)

####### Save File
save(betaSet,
     nb.probes.detectP.cg,
     nb.probes.detectP,
     nb.probes.cg.dmrcate,
     phenoData,
     file =  paste0(data.dir.probes.qc, "/betas_CG_XY_SNPfiltered.rda"))

load(paste0(data.dir.probes.qc, "/betas_CG_XY_SNPfiltered.rda"))
dim(phenoData)
dim(betaSet)

# Normalization

# - Quantile normalization and BMIQ normalization

# Input: 
#   
#   - beta_CG_XY_SNPfiltered_mat.RDS
# - RGSet.RDS
# - pheno_df.RDS
# - full.annot.RDS
# 
# Output: 
#   
#   - bs.csv
# - pheno_df.RDS
# - QNBMIQ.RDS


##### 5. BMIQ ##################################################################
library(wateRmelon)
library(RPMM)
library(sesame)
library(sesameData)

# ```
## Quantile normalization

library(lumi)
betaQN <- lumiN(x.lumi = betaSet, method = "quantile")

### Order annotation in the same order as beta matrix
annotType <- sesameDataGet("EPIC.hg19.manifest") # TODO: Figure this out. 
annotType$designTypeNumeric <- ifelse(annotType$designType == "I",1,2)
### Density plot for type I and type II probes
library(sm)
betaQNCompleteCol1 <- betaQN[complete.cases(betaQN[,1]), ]
annotTypeCompleteCol1 <- annotType[row.names(betaQNCompleteCol1), ]
sm.density.compare(
  betaQNCompleteCol1[,1],
  annotTypeCompleteCol1$designTypeNumeric
)
type12 <- annotType$designTypeNumeric[match(rownames(betaQN),names(annotType))]

### BMIQ
set.seed (946)
betaQN_BMIQ <- apply(
  betaQN, 2,
  function(x){
    norm_ls <- BMIQ(x, design.v = type12, plots = FALSE)
    return (norm_ls$nbeta)
  }
)
saveRDS(betaQN_BMIQ, file.path(data.dir.probes.normalization,"ROSMAP_QNBMIQ.rds"))
saveRDS(phenoData, file.path(data.dir.probes.normalization,"pheno.rds"))

betaQN_BMIQ <- readRDS( file.path(data.dir.probes.normalization,"ROSMAP_QNBMIQ.rds"))
dim(betaQN_BMIQ)
phenoData <- readRDS(file.path(data.dir.probes.normalization,"pheno.rds"))
dim(phenoData)

# Outliers detection - PCA analysis


# Description: 
#   
#   1. estimate standard deviation for each probe
# 2. select most variable probes (e.g. n = 50,000)
# 3. pca analysis
# 4. Filter outliers
# 
# Input: 
#   
#   - QNBMIQ.rds
# - pheno_df.RDS
# 
# Output: 
#   
#   - PCs_usingBetas.csv, 
# - PCA plots
# - QNBMIQ_PCfiltered.RDS
# - pheno_df.RDS
# 
# ```{R, eval = TRUE}
devtools::source_gist("https://gist.github.com/tiagochst/d3a7b1639acf603916c315d23b1efb3e")
beta_mat <- readRDS( file.path(data.dir.probes.normalization,"ROSMAP_QNBMIQ.rds"))
dim(beta_mat)
phenoData <- readRDS(file.path(data.dir.probes.normalization,"pheno.rds"))
dim(phenoData)
identical(colnames(beta_mat), as.character(phenoData$Sample))
# transform to m values # TODO: Check if these manual transformations are OK? Also, ADNI doesn't have Braak stage. 
mvalue_mat <- log2(beta_mat / (1 - beta_mat))
phenoData$braaksc <- as.numeric(as.character(phenoData$braaksc))
phenoData$stage3 <- phenoData$braaksc

phenoData$stage3[phenoData$braaksc <= 2] <- '0-2'
phenoData$stage3[phenoData$braaksc > 2 & phenoData$braaksc < 5] <- '3-4'
phenoData$stage3[phenoData$braaksc >= 5] <- '5-6'

phenoData$sex <- phenoData$msex # Why?????
phenoData$sex[phenoData$msex == 1] <- 'Male'
phenoData$sex[phenoData$msex == 0] <- 'Female'
phenoData$batch <- as.factor(phenoData$batch)
phenoData$apoe_genotype <- as.factor(phenoData$apoe_genotype)
phenoData$race <- as.factor(phenoData$race)
phenoData$spanish <- as.factor(phenoData$spanish)
## 2. first order matrix by most variable probes on top
betaOrd_mat <- OrderDataBySd(beta_mat)
mOrd_mat <- OrderDataBySd(mvalue_mat)
betaOrd_matPlot <- betaOrd_mat[, as.character(phenoData$Sample)]
mOrd_matPlot <- mOrd_mat[, as.character(phenoData$Sample)]
identical(as.character(phenoData$Sample), colnames(betaOrd_matPlot))
identical(as.character(phenoData$Sample), colnames(mOrd_matPlot))
expSorted_mat = betaOrd_mat
pca <- prcomp(
  t(expSorted_mat[1:50000,]),
  scale = TRUE
)
d <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
meanPC1 <- mean(d$PC1)
sdPC1   <- sd(d$PC1)
meanPC2 <- mean(d$PC2)
sdPC2   <- sd(d$PC2)
out3sdPC1_1 <- meanPC1 - 3 * sdPC1
out3sdPC1_2 <- meanPC1 + 3 * sdPC1
out3sdPC2_1 <- meanPC2 - 3 * sdPC2
out3sdPC2_2 <- meanPC2 + 3 * sdPC2
d$outlier_PC1[d$PC1 >= out3sdPC1_1 & d$PC1 <= out3sdPC1_2] <- 0
d$outlier_PC1[d$PC1 < out3sdPC1_1 | d$PC1 > out3sdPC1_2] <- 1
d$outlier_PC2[d$PC2 >= out3sdPC2_1 & d$PC2 <= out3sdPC2_2] <- 0
d$outlier_PC2[d$PC2 < out3sdPC2_1 | d$PC2 > out3sdPC2_2] <- 1
write.csv(d, file.path(data.dir.pca,"ROSMAP_PCs_usingBetas.csv"))
# ```
# 
# ```{R, message = FALSE, results = 'hide'}
### 2. pca plots
library(ggplot2)
library(ggrepel)
phenoData$sample <- phenoData$Sample
# ```
# 
# ```{R, eval = FALSE, include = FALSE}
# beta values
byStage <- plotPCA(
  dataset = "ROSMAP: beta values",
  expSorted_mat = betaOrd_mat,
  pheno = phenoData,
  group_char = "stage3",
  ntop = 50000,
  center = TRUE,
  scale = TRUE
)
bySex <- plotPCA(
  dataset = "ROSMAP: beta values",
  expSorted_mat = betaOrd_mat,
  pheno = phenoData,
  group_char = "msex",
  ntop = 50000,
  center = TRUE,
  scale = TRUE
)
byBatch <- plotPCA(
  dataset = "ROSMAP: beta values",
  expSorted_mat = betaOrd_mat,
  pheno = phenoData,
  group_char = "batch",
  ntop = 50000,
  center = TRUE,
  scale = TRUE
)
byRace <- plotPCA(
  dataset = "ROSMAP: beta values",
  expSorted_mat = betaOrd_mat,
  pheno = phenoData,
  group_char = "race",
  ntop = 50000,
  center = TRUE,
  scale = TRUE
)
byApoe <- plotPCA(
  dataset = "ROSMAP: beta values",
  expSorted_mat = betaOrd_mat,
  pheno = phenoData,
  group_char = "apoe_genotype",
  ntop = 50000,
  center = TRUE,
  scale = TRUE
)
# M values
byStage <- plotPCA(
  dataset = "ROSMAP: M values",
  expSorted_mat = mOrd_mat,
  pheno = phenoData,
  group_char = "stage3",
  ntop = 50000,
  center = TRUE,
  scale = TRUE
)
bySex <- plotPCA(
  dataset = "ROSMAP: M values",
  expSorted_mat = mOrd_mat,
  pheno = phenoData,
  group_char = "msex",
  ntop = 50000,
  center = TRUE,
  scale = TRUE
)
byBatch <- plotPCA(
  dataset = "ROSMAP: M values",
  expSorted_mat = mOrd_mat,
  pheno = phenoData,
  group_char = "batch",
  ntop = 50000,
  center = TRUE,
  scale = TRUE
)
byRace <- plotPCA(
  dataset = "ROSMAP: M values",
  expSorted_mat = mOrd_mat,
  pheno = phenoData,
  group_char = "race",
  ntop = 50000,
  center = TRUE,
  scale = TRUE
)
byApoe <- plotPCA(
  dataset = "ROSMAP: M values",
  expSorted_mat = mOrd_mat,
  pheno = phenoData,
  group_char = "apoe_genotype",
  ntop = 50000,
  center = TRUE,
  scale = TRUE
)

## Filter samples by PCA, SAVE files
noOutliers <- d[which(d$outlier_PC1 == 0 & d$outlier_PC2 == 0), ]
betaQN_BMIQ_PCfiltered <- beta_mat[, rownames(noOutliers)]
dim(betaQN_BMIQ_PCfiltered)
saveRDS(betaQN_BMIQ_PCfiltered, file.path(data.dir.pca,"ROSMAP_QNBMIQ_PCfiltered.RDS"))
phenoData <- phenoData [phenoData$Sample %in% rownames(noOutliers) ,]
dim(phenoData)
saveRDS(phenoData, file.path(data.dir.pca,"pheno_PCfiltered.RDS"))

# Summary after QC steps

## Data and metadata

betaQN_BMIQ_PCfiltered <- readRDS(file.path(data.dir.pca,"ROSMAP_QNBMIQ_PCfiltered.RDS")) 
nb.samples.with.clinical.after.pca <- ncol(betaQN_BMIQ_PCfiltered)
dim(betaQN_BMIQ_PCfiltered)
phenoData <- readRDS(file.path(data.dir.pca,"pheno_PCfiltered.RDS")) 
dim(phenoData)

## Numbers of samples and probes removed in each step

df.samples <- data.frame("Number of samples" =  c(nb.samples, 
                                                  nb.samples.bc.filtered,
                                                  nb.samples.with.slide,
                                                  nb.samples.with.clinical, 
                                                  nb.samples.with.clinical.after.pca),
                         "Description" = c("total number of samples",
                                           "samples with bisulfate conversion > 88",
                                           "samples with slide information",
                                           "samples with clinical data",
                                           "Samples after PCA"),
                         "Difference" = c("-",
                                          nb.samples.bc.filtered - nb.samples ,
                                          nb.samples.with.slide - nb.samples.bc.filtered ,
                                          nb.samples.with.clinical - nb.samples.with.slide ,
                                          nb.samples.with.clinical.after.pca - nb.samples.with.clinical)
)    
df.samples                     
# Create summary table
df.probes <- data.frame("Number of probes" = c(nb.probes,
                                               nb.probes.detectP, 
                                               nb.probes.detectP.cg,
                                               nb.probes.cg.dmrcate),
                        "Description" = c("total number of probes in raw data",
                                          "detection P < 0.01",
                                          "only probes that start with cg",
                                          "DMRcate"),
                        "Difference" = c("-",
                                         nb.probes.detectP - nb.probes ,
                                         nb.probes.detectP.cg - nb.probes.detectP,
                                         nb.probes.cg.dmrcate - nb.probes.detectP.cg)
)
df.probes
save(df.samples,df.probes,file = file.path(data.dir.table, "ROSMAP_table.rda"))

# Compute neuron proportion

# 
# Data from  https://www.tandfonline.com/doi/full/10.4161/epi.23924
# 
# - Input: London_PFC_QNBMIQ_PCfiltered.RDS, pheno107_PFC_df.RDS
# - Output: pheno107_PFC_withNeuronProp_df.RDS

objects <- load("~/data/CETS/cets.RData") # TODO: Check if there's either an EPIC array version of this file or if I need to translate the signals??? Or just delete this portion?
objects

## Get reference profile from Caucasions + controls 
idx <- list(
  controlNeuron = pdBrain$celltype == "N" & pdBrain$diag == "Control" & pdBrain$ethnicity == "Caucasian",
  controlGlia   = pdBrain$celltype == "G" & pdBrain$diag == "Control" & pdBrain$ethnicity == "Caucasian"
)
refProfile <- getReference(brain, idx)
head(refProfile)
##### 2. Estimate proportions of neurons in PFC samples ########################
### Limit to 10,000 cpgs in the refProfile dataset
pfc <- readRDS(file.path(data.dir.pca,"ROSMAP_QNBMIQ_PCfiltered.RDS")) 
selected <- rownames(pfc) %in% rownames(refProfile)
table(selected)
pfc.refcpgs <- pfc[selected, ] 
### Estimate proportion of neurons
prop <- data.frame(estProportion(pfc.refcpgs, profile = refProfile))
colnames(prop) <- "prop.neuron"
##### 3. Merge pfc.refcpgs with phenotype file #################################
pheno <- readRDS(file.path(data.dir.pca,"pheno_PCfiltered.RDS"))
pheno_final <- merge(
  pheno,
  prop,
  by.x = "Sample",
  by.y = "row.names"
)
saveRDS(pheno_final, paste0(data.dir.neuron, "pheno_PFC_withNeuronProp_df.RDS"))
# ```
# 
# ```{R, include = FALSE, eval = FALSE}
pheno_final <- readRDS(paste0(data.dir.neuron, "pheno_PFC_withNeuronProp_df.RDS"))
# ```


# Linear regression by cpgs Methylation 

## Import datasets

# ```{R}
beta_mat <- readRDS(file.path(data.dir.pca,"ROSMAP_QNBMIQ_PCfiltered.RDS")) 
pheno_df <- readRDS(paste0(data.dir.neuron, "pheno_PFC_withNeuronProp_df.RDS")) 
# ```

## Test all regions

# ```{R, eval = TRUE}
### Compute M values
mval_mat <- log2(beta_mat / (1 - beta_mat))
pheno_df <- pheno_df[match(colnames(mval_mat),pheno_df$Sample),]
identical(pheno_df$Sample, colnames(mval_mat))
pheno_df$age.brain <- pheno_df$age_death
pheno_df$sex <- as.factor(pheno_df$sex)
pheno_df$slide <- as.factor(pheno_df$Sentrix_ID) ### ???
pheno_df$batch <- as.factor(pheno_df$batch)
# If rosmap cohort, don't forget batch effect
is(pheno_df$braaksc,"numeric")
is(pheno_df$age.brain,"numeric")
is(pheno_df$prop.neuron,"numeric")
#str(pheno_df)
# ```

# ```{R,  eval = FALSE}
predictors_char <- "braaksc" # Change to Diagnosis
covariates_char <- c("age.brain", "sex", "prop.neuron", "slide", "batch")
doParallel::registerDoParallel(cores = 20)
devtools::source_gist("https://gist.github.com/tiagochst/d3a7b1639acf603916c315d23b1efb3e")
results_ordered_df <- plyr::adply(mval_mat,1, function(row){
  
  sumOneRegion_df <- as.data.frame(t(row))
  
  result <- TestSingleRegion(
    predictors_char = predictors_char,
    covariates_char = covariates_char,
    pheno_df = pheno_df,
    sumOneRegion_df = sumOneRegion_df
  )
  result
}, .progress = "time",.parallel = TRUE,.id = "cpg")
colnames(results_ordered_df)[1] <- "cpg"
identical(row.names(mval_mat), results_ordered_df$cpg %>% as.character())
results_ordered_df$fdr <- p.adjust(
  results_ordered_df$pValue,
  method = "fdr"
)
write.csv(
  results_ordered_df,
  paste0(data.dir.single.cpg.pval, "ROSMAP_PFC_single_cpg_pVal_df.csv"),
  row.names = FALSE
)
# ```

# ```{R}
results_ordered_df <- readr::read_csv(paste0(data.dir.single.cpg.pval, "ROSMAP_PFC_single_cpg_pVal_df.csv"))
results_ordered_df
# ```


# Linear regression by regions median Methylation 

## Residuals control and coMethylated Regions

# 1. Take residuals
# 2. Find co-methylated regions
# 
# Input: 
#   
#   - QNBMIQ_PCfiltered
# - pheno_withNeuronProp_df
# 
# Output: 
#   
#   - QNBMIQ_PCfiltered_mvalResiduals
# - residuals_cometh_ls

###  Take residuals

# ```{R, eval = FALSE}
##### 1. Import datasets #######################################################
beta_mat <- readRDS(grep("QNBMIQ",dir(data.dir.pca,full.names = T),ignore.case = T,value = T)) 
pheno_df <- readRDS(dir(data.dir.neuron,full.names = T)) 
### Compute M values
mvalue_mat <- log2(beta_mat / (1 - beta_mat))
### Reorder samples based on pheno_df
mvalue_mat <- mvalue_mat[, pheno_df$Sample]
identical(colnames(mvalue_mat),  pheno_df$Sample)
### Take residuals
lmF <- function(mval){
  fitE <- lm(
    as.numeric(mval) ~ age_death + msex + prop.neuron + as.character(Slide) + batch, #add batch if rosmap
    data = pheno_df,
    na.action = na.exclude
  )
  residuals (fitE)
}
doParallel::registerDoParallel(cores = 4)
resid <- plyr::adply(mvalue_mat,1,.fun = lmF,.progress = "time",.parallel = TRUE)
rownames(resid) <- resid[,1]
resid[,1] <- NULL
colnames(resid) <- colnames(mvalue_mat)
dim(resid)
dim(mvalue_mat)
### Save dataset
saveRDS(
  resid,
  paste0(data.dir.residuals, 
         "ROSMAP_QNBMIQ_PCfiltered_mvalResiduals.RDS"
  )
)
# ```

### Find co-methylated regions

# ```{R, eval = FALSE}
### Import datasets
mvalue_residuals_mat <- readRDS(
  paste0(data.dir.residuals, 
         "ROSMAP_QNBMIQ_PCfiltered_mvalResiduals.RDS"
  )
)
### Call in functions
library(coMethDMR)
library(BiocParallel)
probes.cluster.all <- coMethDMR::getPredefinedCluster(arrayType = "450k",
                                                      clusterType = "regions")
### Find co-methylated clusters
ncores <- 6
coMeth_ls <- CoMethAllRegions(
  dnam = mvalue_residuals_mat,      
  betaToM = FALSE,                   
  CpGs_ls = probes.cluster.all,
  arrayType = "EPIC",
  rDropThresh_num = 0.4,
  minPairwiseCorr = NULL,
  method = "spearman",             
  returnAllCpGs = TRUE,              
  output = "all",
  nCores_int = ncores,
  progressbar = FALSE
)
saveRDS(
  coMeth_ls,
  paste0(data.dir.residuals, 
         "ROSMAP_residuals_cometh_input_ls.RDS")
)
# ```


## Linear regression by regions median Methylation 

# 1. Calculate medians by cluster and sample
# 2. linear regression
# 
# Input: 
#   
#   - QNBMIQ_PCfiltered,
# - pheno_withNeuronProp_df
# - residuals_cometh_input_ls
# 
# Output: 
#   
#   - info_df
# - mediansMval_df
# - linear_df

### Calculate medians by cluster and sample

# ```{R, eval = FALSE}
### Import datasets
beta_mat <- beta_mat <- readRDS(dir(data.dir.pca, pattern = "QNBMIQ", full.names = TRUE))
pheno_df <- readRDS(dir(data.dir.neuron,full.names = T)) 
mval_mat <- log2(beta_mat / (1 - beta_mat)) %>% as.matrix()
coMeth_ls <- readRDS(
  paste0(data.dir.residuals, 
         "ROSMAP_residuals_cometh_input_ls.RDS"
  )
)
library(tidyverse)
cohort <- "ROSMAP"
data.dir <- "~/data/ROSMAPmethNew"
data.dir.table <- "~/data/ROSMAPmethNew/datatable" 
data.dir.raw <- file.path(data.dir,"step1_download/") 
data.dir.raw.idat <- "~/data/ROSMAPmeth/"
data.dir.raw.metadata <- file.path(data.dir.raw, "Metadata")
data.dir.read <- file.path(data.dir,"step2_read_minfi/") 
data.dir.bsfilter <- file.path(data.dir,"step3_bsConvFilter/") 
data.dir.clinical.filter <- file.path(data.dir,"step4_clinical_available_filtering/") 
data.dir.probes.qc <- file.path(data.dir,"step5_probesQC_filtering/") 
data.dir.probes.normalization <- file.path(data.dir,"step6_normalization/") 
data.dir.pca <- file.path(data.dir,"step7_pca_filtering/") 
data.dir.neuron <- file.path(data.dir,"step8_neuron_comp/") 
data.dir.single.cpg.pval <- file.path(data.dir,"step9_single_cpg_pval/") 
data.dir.residuals <- file.path(data.dir,"step10_residuals/") 
data.dir.median <- file.path(data.dir,"step11_median/") 
for(p in grep("data",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

beta_mat <- beta_mat <- readRDS(dir(data.dir.pca, pattern = "QNBMIQ", full.names = TRUE))
pheno_df <- readRDS(dir(data.dir.neuron,full.names = T)) 
mval_mat <- log2(beta_mat / (1 - beta_mat)) %>% as.matrix()
coMeth_ls <- readRDS(
  paste0(data.dir.residuals, 
         "ROSMAP_residuals_cometh_input_ls.RDS"
  )
)

print("Imported fine")
### Create info dataset

# print(coMeth_ls)

#joined <- as.data.frame(rbindlist(coMeth_ls))

# input_cometh <- data.frame(
#   inputRegion = coMeth_ls$inputRegion_chr,
#   nCoMethRegion = coMeth_ls$nCoMethRegions_num,
#   coMethRegion = names(coMeth_ls$coMeth_ls),
#   nCpGs = unlist(lapply(coMeth_ls$coMeth_ls, length), use.names = FALSE),
#   stringsAsFactors = FALSE
# )
input_cometh <- NULL

for (i in 1:31199){
  temp = coMeth_ls[[i]]
  row <- data.frame(input_region=(temp$Region), nCoMethRegion=length(temp$Region), coMethRegion=(temp$CpG), nCpGs = length(temp$CpG), stringsAsFactors = F)
  rbind(input_cometh, row) -> input_cometh
  print(i)
}

print("Input cometh made")

#print(input_cometh)

input_cometh_nodup <- input_cometh[
  !duplicated(input_cometh$coMethRegion),
  ]
#print(input_cometh_nodup)


colnames(input_cometh_nodup) <- c(
  paste0(cohort, "_inputRegion"),
  paste0(cohort, "_nCoMethRegion"),
  paste0(cohort, "_coMethRegion"),
  paste0(cohort, "_nCpGs")
)

print("Formated weird")
saveRDS(
  input_cometh_nodup,
  paste0(data.dir.median, cohort, "_info_df.rds")
)
print("Saved")
### Take median of probes in each cluster for each sample
filename <-  paste0(data.dir.median, cohort, "_mediansMval_df.rds")
library(robustbase)

print("Start part 2")
## This is giving you some trouble. 
# test <- do.call(rbind, coMeth_ls)

mval_mat <- mval_mat[rownames(mval_mat) %in% unlist(input_cometh$CpG),]

medianMval.df <- NULL

for (i in 1:31199){
  temp = coMeth_ls[[i]]
  name = names(coMeth_ls)[[i]]
  cpgList = c(temp$CpG)
  row <- as.matrix(data.frame(colMedians(as.matrix(mval_mat[as.character(cpgList),]),na.rm=T)))
  colnames(row) <- name
  row <- t(row)
  rbind(medianMval.df, row) -> medianMval.df
  print(i)
}

# if(!file.exists(filename)){
#   medianMval.df <- plyr::ldply(
#     coMeth_ls$CpG[!duplicated((coMeth_ls$CpG))],
#     function(probes){
#       colMedians(as.matrix(mval_mat[as.character(probes),]),na.rm=T)
#     },
#     .progress = "time"
#   )
#   medianMval.df$.id <- NULL
#   colnames(medianMval.df) <- colnames(mval_mat)
#   saveRDS(medianMval.df, file = filename)
# } else {
#   print("File exists")
#   medianMval.df <- readRDS(filename)
# }
# ```
saveRDS(medianMval.df, file = filename)
### Test all regions -- linear regressions

# ```{R, eval = TRUE}
### Import datasets
info_df <- readRDS(
  paste0(data.dir.median, cohort, "_info_df.rds")
)
mediansMval_df <- readRDS(
  paste0(data.dir.median, cohort, "_mediansMval_df.rds")
)
pheno_df <- readRDS(dir(data.dir.neuron,full.names = T)) 
### Check variables before fitting model
pheno_df$Sample <- pheno_df$Sample_Group
identical(pheno_df$Sample, colnames(mediansMval_df))
#str(pheno_df)
pheno_df$msex <- as.factor(pheno_df$msex)
pheno_df$Slide <- as.factor(pheno_df$Slide)
# If rosmap cohort, don't forget batch effect
pheno_df$age.brain <- pheno_df$age_death
pheno_df$sex <- as.factor(pheno_df$sex)
pheno_df$slide <- as.factor(pheno_df$Sentrix_ID) ### ???
pheno_df$batch <- as.factor(pheno_df$batch)
# If rosmap cohort, don't forget batch effect
is(pheno_df$braaksc,"numeric")
is(pheno_df$age.brain,"numeric")
is(pheno_df$prop.neuron,"numeric")
#str(pheno_df)
# ```

# ```{R, eval = FALSE}
predictors_char <- "braaksc" # change to AD diagnosis
covariates_char <- c("age.brain", "msex", "prop.neuron", "Slide", "batch") # no prop.neuron potentially?
devtools::source_gist("https://gist.github.com/tiagochst/d3a7b1639acf603916c315d23b1efb3e")
res_df <- TestAllRegions_noInfo(
  predictors_char = predictors_char,
  covariates_char = covariates_char,
  pheno_df = pheno_df,
  summarizedRegions_df = mediansMval_df,
  cores = 1 
)
colnames(res_df) <- c(
  paste0(cohort, "_estimate"),
  paste0(cohort, "_se"),
  paste0(cohort, "_pVal"),
  paste0(cohort, "_fdr")
)
res_withInfo_df <- cbind(info_df, res_df)
saveRDS(
  res_withInfo_df,
  paste0(data.dir.median, cohort, "_linear_df.rds")
)
# ```

# ```{R}
#file <- dir(data.dir.median,pattern = paste0(".*linear_df"),
            # recursive = T,
            # full.names = TRUE,
            # ignore.case = T)
#file
#res_withInfo_df <- readRDS(file)
#dim(res_withInfo_df)
#res_withInfo_df

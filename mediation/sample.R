# Brian Lee (leebn@sas)
# 10 Feb 2021
# PLEASE do not distribute this copy outside of the Shen Lab/AD Omics Project!

# Sample implementation of a mediation analysis implemented in R via the mediation library!

library(mediation)

# If not installed:
# install.packages("mediation")

# These lookup tables available upon request; data taken directly from ADNI QT-PAD. 
fileList = list("/Volumes/ShenLab/11Feb/bl_lookuptable.csv", 
                "/Volumes/ShenLab/11Feb/m06_lookuptable.csv", 
                "/Volumes/ShenLab/11Feb/m12_lookuptable.csv")

# NOTE: Predictor (X) is SNP allelic dosage, Outcome (Y) is Cognitive Measurement value, Mediators (M) are imaging measures.
snps <- list("APOC1", "APOE", "TOMM40", "Chr12")
cog <- list("ADAS13", "CDRSB", "RAVLTlearning", "MMSE", "FAQ")
img <- list("FDG", "AV45", "WholeBrain", "Hippocampus", "Entorhinal", "Fusiform", 
            "MidTemp", "Ventricles")
            
 
for (file in fileList) {
  fileCounter = 1
  
  # DF to hold effect sizes + p values for significant results; see here for headers.
  output <- data.frame(Time=character(), Imaging=character(), Cognitive=character(), Indirect=numeric(), 
                       Direct=numeric(), Prop_Mediating=numeric(), Indirect_P=numeric(), Direct_P=numeric(),
                       Prop_Mediating_P=numeric(), stringsAsFactors=F)
  df <- read.csv(file)
  rowNum = 1
  for (snp in snps) {
    for (c in cog) {
      for (i in img) {
        temp <- df[c(snp, c, i, "AGE", "PTEDUCAT", "PTGENDER")] # data cleaning
        temp <- temp[complete.cases(temp),]
        
        medForm = paste0(i, " ~ ", snp, " + AGE + PTEDUCAT + factor(PTGENDER)", sep="") # Linear regression formula
        outForm = paste0(c, " ~ ", snp, " + ", i, " + AGE + PTEDUCAT + factor(PTGENDER)", sep="") # Linear regression formula
        
        med.fit <- lm(medForm, data=temp) # Relationship btwn SNP and Mediator (imaging measure)
        out.fit <- lm(outForm, data=temp) # Relationship btwn SNP and outcome (cog) with mediator (imaging) as covariate
        med.out <- mediate(med.fit, out.fit, treat=snp, mediator=c, robustSE=T, sims=1000) # Mediation analysis
        
        # Forgive the lazy code ...
        if (fileCounter == 1) {
          time = "bl"
        } else if (fileCounter == 2) {
          time = "m06"
        } else {
          time = "m12"
        }
        
        if (med.out$d.avg.p <= .05) {
          output[rowNum, ] <- list(time, i, c, med.out$d.avg, med.out$z.avg, med.out$n.avg, med.out$d.avg.p, 
                                 med.out$z.avg.p, med.out$n.avg.p) 
          rowNum = rowNum + 1
          print(rowNum)
        }
      }
    }
  }
  write.csv(output, paste0("/Volumes/ShenLab/11Feb/", time, "_mediation.csv", sep=""), quote=F)
  fileCounter = fileCounter + 1
  print(filCounter)
}

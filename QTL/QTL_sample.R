# Brian Lee (leebn@sas for questions)
# 27 Jan 2021
# MatrixEQTL: Shabalin, A.A. Matrix eQTL: Ultra fast eQTL analysis via large matrix operations. Bioinformatics 28, no. 10 (2012): 1353-1358.
# Code from website: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# *QTL (i.e. eQTL/mQTL) +  GWAS MatrixEQTL Pipeline
# Highly reccomend running using Microsoft R/using fast BLAS
# If on OS X: https://statistics.berkeley.edu/computing/faqs/linear-algebra-and-parallelized-linear-algebra-using-blas

library(RCurl) # Just to be safe when using CUBIC
library(MatrixEQTL)
useModel = modelLINEAR
SNP_file_name = ".../geno.txt"
expression_file_name = ".../pheno.txt"
covariates_file_name = ".../cov.txt"
output_file_name = ".../out.txt"
pvOutputThreshold = (0.05) / (52148 * 750174) # Note: Should use Bonferroni correction: (.05)/{(num(SNPs) * num(phenotypes)}
errorCovariance = numeric();

## Load SNP's

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 5000;      # read file in slices of 5,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 5000;      # read file in slices of 5,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

write.csv(me$all$eqtls, output_file_name, row.names=F, quote=F)

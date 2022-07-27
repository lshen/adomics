################################################################################
# 1. Download 1000 Genome genotypes.
#    Extract dosage from APOE region (1M bp padding) using bcftools.
################################################################################
# 1000 Genome genotype data downloaded from:
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/
# ALL.chr19.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz
#
# Extract APOE gene
# https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000130203;r=19:45409011-45412650
if (!file.exists("../data/APOE_region_dosage_matrix.rds")) {
  if (!file.exists("../data/APOE_region.vcf")) {
    system("bcftools view -r 19:44409011-46412650 ../data/ALL.chr19.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz > ../data/APOE_region.vcf")
  }
  if (!file.exists("../data/APOE_region_dosage.txt")) {
    system("bcftools +dosage -r 19:44409011-46412650 ../data/ALL.chr19.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz > ../data/APOE_region_dosage.txt")
  }
}

################################################################################
# 2. Create a dosage matrix with 1092 individuals (cols) and 6157 SNPs (rows).
################################################################################
library(tidyverse)
genotype = read_tsv("../data/APOE_region.vcf", comment = "#", col_names = FALSE)
dosage = read_tsv("../data/APOE_region_dosage.txt")
stopifnot(all.equal(genotype$X2, dosage$`[2]POS`))
dosage$rsid = genotype$X3  # rsid
# SNPs with NA rsids are removed:
dosage_matrix = as.matrix(dosage[dosage$rsid != ".", 5:1096])  # 1092 individuals
rownames(dosage_matrix) = dosage$rsid[dosage$rsid != "."]
dim(dosage_matrix)  # 27930 x 1092
# Some of the SNPs are missing in the LD reference panel and will be removed:
library(ieugwasr)
dosage_matrix = dosage_matrix[ld_reflookup(rownames(dosage_matrix), pop = "EUR"), ]
dim(dosage_matrix)  # 6157 x 1092
saveRDS(dosage_matrix, "../data/APOE_region_dosage_matrix.rds")

# WARNING: Old version. Do not use. 
# ################################################################################
# # 3. Simulate data that are specific to each setup
# ################################################################################
# set.seed(1234)
# # sample partition
# library(ieugwasr)
# dosage_matrix_path =  "../data/APOE_region_dosage_matrix.rds"
# stopifnot(file.exists(dosage_matrix_path))
# dosage_matrix = readRDS(dosage_matrix_path)
# n = sample_size = 500
# samples_all = sample(ncol(dosage_matrix), sample_size*2, replace=FALSE)
# GWAS_sample_indices = samples_all[seq_len(sample_size)]
# QTL_sample_indices  = samples_all[sample_size + seq_len(sample_size)]
# 
# EUR_path = "../data/EUR"
# if (!file.exists(paste0(EUR_path, ".bed"))) {
#   system("cd ../data/ && wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz &&  xzvf 1kg.v3.tgz EUR.bed EUR.bim EUR.fam")
# }
# 
# # plink_path = "/appl/plink-1.9-20210416/plink"
# # devtools::install_github("explodecomputer/plinkbinr")
# plink_path = plinkbinr::get_plink_exe()
# is_plink_available = file.exists(plink_path)
# if (!is_plink_available) {
#   plink_path = EUR_path = NULL
# }
# 
# # non-zero SNP-exposure effects: horizontal 1, 2, 3 and vertical
# # SNP-exposure effect
# IV_strength = 0.5
# effect = 0.1
# outcome = "continuous"
# xQTL_from_one_sample = TRUE
# overlap = 0
# number_of_SNPs_in_laplace_dist = 100
# SNP_exposure_effects_horizontal_1 = sample(c(
#   rep(0, nrow(dosage_matrix) - number_of_SNPs_in_laplace_dist),
#   rexp(number_of_SNPs_in_laplace_dist, rate = 1) *
#     sample(c(-1,1), number_of_SNPs_in_laplace_dist, replace=TRUE)
# ))
# SNP_exposure_effects_horizontal_2 = sample(c(
#   rep(0, nrow(dosage_matrix) - number_of_SNPs_in_laplace_dist),
#   rexp(number_of_SNPs_in_laplace_dist, rate = 1) *
#     sample(c(-1,1), number_of_SNPs_in_laplace_dist, replace=TRUE)
# ))
# SNP_exposure_effects_horizontal_3 = sample(c(
#   rep(0, nrow(dosage_matrix) - number_of_SNPs_in_laplace_dist),
#   rexp(number_of_SNPs_in_laplace_dist, rate = 1) *
#     sample(c(-1,1), number_of_SNPs_in_laplace_dist, replace=TRUE)
# ))
# SNP_exposure_effects_vertical = sample(c(
#   rep(0, nrow(dosage_matrix) - number_of_SNPs_in_laplace_dist),
#   rexp(number_of_SNPs_in_laplace_dist, rate = 1) *
#     sample(c(-1,1), number_of_SNPs_in_laplace_dist, replace=TRUE)
# ))
# 
# # select instrumental variables: horizontal and vertical
# for (pleiotropy in c("horizontal", "vertical")) {
#   # z_matrix_in_GWAS = t(dosage_matrix[, samples_in_GWAS])
#   # z_matrix_in_QTL  = t(dosage_matrix[, samples_in_QTL])
#   z_matrix = t(dosage_matrix[, c(GWAS_sample_indices, QTL_sample_indices)])
#   u = rnorm(2*n)
#   
#   # determine the number of the SNPs
#   p = ncol(z_matrix)
#   
#   # SNP-exposure effect
#   generate_SNP_effects = function(zmat, wt) {
#     zmat[] = t(t(zmat) * wt)
#     rowSums(zmat)
#   }
#   
#   if (pleiotropy == "horizontal") {
#     z = generate_SNP_effects(z_matrix, SNP_exposure_effects_horizontal_1)
#     z_first_half = generate_SNP_effects(z_matrix[,1:(floor(p/2))],
#                                         SNP_exposure_effects_horizontal_2[1:(floor(p/2))])
#     z_second_half = generate_SNP_effects(z_matrix[,(floor(p/2)+1):p],
#                                          SNP_exposure_effects_horizontal_3[(floor(p/2)+1):p])
#     x1 = 1/10 * IV_strength * z + u + rnorm(2*n) 
#     x2 = - 1/10 * IV_strength * z_first_half - u + rnorm(2*n)
#     x3 = 1/10 * IV_strength * z_second_half - u + rnorm(2*n)
#     v = effect * (x1 + x2 + x3) + u
#     odds = exp(-2 + v)
#   } else if (pleiotropy == "vertical") {
#     z = generate_SNP_effects(z_matrix, SNP_exposure_effects_vertical)
#     x1 = 1/10 * IV_strength * z + rnorm(2*n) + u
#     x2 = 2 * x1 - u + rnorm(2*n)
#     x3 = - x2 + u + rnorm(2*n)
#     v = effect * x3 - u + rnorm(2*n)
#     odds = exp(-2 + v)
#   } else stop("Unrecognized pleiotropy")
#   
#   
#   if (xQTL_from_one_sample) {
#     if (outcome == "continuous") {
#       y = v + rnorm(2*n)
#       samples_in_GWAS = 1:n
#       samples_in_xQTL = c(
#         seq_len(round(n*overlap)),
#         n + seq_len(n-round(n*overlap)))
#     }
#     samples_in_xQTL1 = samples_in_xQTL2 = samples_in_xQTL3 = samples_in_xQTL
#   } else {
#     if (outcome == "continuous") {
#       y = v + rnorm(2*n)
#       samples_in_GWAS = 1:n
#       samples_in_xQTL1 = n + seq_len(n)
#       samples_in_xQTL2 = n*2 + seq_len(n)
#       samples_in_xQTL3 = n*3 + seq_len(n)
#     }
#   }
#   
#   x1 = x1[samples_in_xQTL1]
#   x2 = x2[samples_in_xQTL2]
#   x3 = x3[samples_in_xQTL3]
#   y0 = y[samples_in_GWAS]
#   z1 = z_matrix[samples_in_xQTL1, ]
#   z2 = z_matrix[samples_in_xQTL2, ]
#   z3 = z_matrix[samples_in_xQTL3, ]
#   z0 = z_matrix[samples_in_GWAS, ]
#   
#   zx_est = zx_se = matrix(NA, nrow=3, ncol=p)
#   zy_est = zy_se = rep(NA, p)
#   for (k in 1:p) {
#     x1_lm_summary = summary(lm(x1~z1[,k]))
#     x2_lm_summary = summary(lm(x2~z2[,k]))
#     x3_lm_summary = summary(lm(x3~z3[,k]))
#     
#     x1_lm_coefficients = x1_lm_summary$coefficients
#     x2_lm_coefficients = x2_lm_summary$coefficients
#     x3_lm_coefficients = x3_lm_summary$coefficients
#     
#     zx_est[1, k] = x1_lm_coefficients[2, "Estimate"]
#     zx_se[1, k] = x1_lm_coefficients[2, "Std. Error"]
#     zx_est[2, k] = x2_lm_coefficients[2, "Estimate"]
#     zx_se[2, k] = x2_lm_coefficients[2, "Std. Error"]
#     zx_est[3, k] = x3_lm_coefficients[2, "Estimate"]
#     zx_se[3, k] = x3_lm_coefficients[2, "Std. Error"]
#     
#     if (outcome == "continuous") {
#       y0_lm_coefficients = summary(lm(y0~z0[,k]))$coefficients
#       
#       zy_est[k] = y0_lm_coefficients[2, "Estimate"]
#       zy_se[k] = y0_lm_coefficients[2, "Std. Error"]
#     }
#   }
#   
#   # select instrumental variables using ieugwasr
#   # we have plink_path
#   # get r2 matrix using ieugwasr::ld_matrix() and local references
#   min_p = colMins(2*pnorm(zx_est/zx_se, lower.tail = FALSE))
#   min_p[min_p > 1] = 1
#   # EUR_path = plink_path = NULL
#   ld_clump_result = ieugwasr::ld_clump(dat = data.frame(rsid = colnames(z_matrix), pval = min_p),
#                                        clump_kb = 100,
#                                        clump_r2 = 0.2,
#                                        clump_p = max(1e-4, min(min_p)),
#                                        bfile = EUR_path,
#                                        plink_bin = plink_path)
#   ld_matrix_signed = ieugwasr::ld_matrix(ld_clump_result$rsid,
#                                          with_allele = FALSE,
#                                          pop = "EUR",
#                                          bfile = EUR_path,
#                                          plink_bin = plink_path)
#   selected_SNPs = ld_clump_result$rsid
#   which_selected_SNPs = match(ld_clump_result$rsid, colnames(z_matrix))
#   if (pleiotropy == "horizontal") {
#     ld_matrix_signed_horizontal = ld_matrix_signed
#     which_selected_SNPs_horizontal = which_selected_SNPs
#   } else if (pleiotropy == "vertical") {
#     ld_matrix_signed_vertical = ld_matrix_signed
#     which_selected_SNPs_vertical = which_selected_SNPs
#   }
# }
# 
# save(GWAS_sample_indices, QTL_sample_indices,
#      SNP_exposure_effects_horizontal_1,
#      SNP_exposure_effects_horizontal_2,
#      SNP_exposure_effects_horizontal_3,
#      SNP_exposure_effects_vertical,
#      ld_matrix_signed_horizontal, which_selected_SNPs_horizontal,
#      ld_matrix_signed_vertical, which_selected_SNPs_vertical,
#      file="precomputed_IVs.RData")